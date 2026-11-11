from math import sqrt, cos, tan, radians, pi, exp
KG, FT, FT2, IN, G = 2.20462262, 3.2808399, 10.7639104, 39.3700787, 9.80665
AR, TAPER = 12.0, 0.292604704
SWEEP_C4, SWEEP_HT, SWEEP_VT = 30.77, 30.8, 34.5       # deg   
TC_ROOT, TC_VT = 0.121, 0.12                           # t/c   
CSW_RATIO = 35.57405256 / 169.4002503                  # S_csw / S_w
N_Z = 3.75                                             # ultimate load factor
W_S = 84363.16389 * G / 169.4002503                    # design wing loading [N/m2]
W0_GUESS = 84363.16389                                 # Class I MTOW [kg]
PAYLOAD = 18960.0                                      # [kg]
F_FUEL = 0.2907                                         # Class I f_f = W_fuel/W_TO: starting value of the iteration
r_0 = 1
E_OSWALD = 0.8                                        # Oswald factor, clean config (a higher value used, that corresponds to the performance of similar aircraft)  
ALPHA = 1.98     # alpha = 2 g R_0 / (eta_eng eta_p e_f sqrt(pi AR e))   
BETA = 0.01275      # beta  = g (1 + f_con) h_E / (0.7 eta_eng eta_p e_f)  
CD0_REF = {
    "wing":     0.0077,   # C_f FF IF S_wet / S_ref, wing                          -> scales r^-0.075 (PLACEHOLDER)
    "tails":    0.0016,   # horizontal + vertical tail                              -> scales r^0.39   (PLACEHOLDER)
    "fuselage": 0.00379,   # fixed hardware                                          -> scales 1/r      (PLACEHOLDER)
    "nacelles": 0.0016,   # nacelles + pylons, fixed hardware                       -> scales 1/r      (PLACEHOLDER)
    "gear":     0.000,   # dCD gear, retracted during cruise -> 1/r        (PLACEHOLDER)
    "upsweep":  0.001,   # (D/q)/S_ref, Eq. 4.14                                   -> scales 1/r      (PLACEHOLDER)
    "base":     0.000,   # (D/q)/S_ref, Eq. 4.15                                   -> scales 1/r      (PLACEHOLDER)
    "flap":     0.000,   # Eq. 4.10, flaps retracted in cruise -> 0                -> constant
    "wave":     0.001,   # Eq. 4.16/4.17 at M = 0.85 (needs M_crit, M_DD)          -> constant        (PLACEHOLDER)
}
K_EXCR = 0.97                                          # 3 % excrescence & leakage: C_D0 = subtotal / 0.97
WFUEL_LANDING = 0.943                                  # Eq. 4.26: W_L = W_TO - 0.943 W_fuel
FUS_LENGTH, FUS_DEPTH, FUS_WIDTH, FUS_SHELL = 44.92, 3.75, 3.85, 438.0   # m, m, m, m2 (shell = S_f)(PLACEHOLDERS)
V_H, V_V = 1.15, 0.086                                   # tail volume coefficients
L_H, L_V = 21.36, 21.36                                 # tail arms [m](PLACEHOLDERS)
A_H, A_V = 4, 1.6                                       # tail aspect ratios 
T_TAIL = False
GEAR_L_MAIN, GEAR_L_NOSE = 1.9, 1.65                   # strut lengths [m] (PLACEHOLDERS)
N_MW, N_MSS, N_NW = 4, 2, 2                            # main wheels, main struts, nose wheels (PLACEHOLDERS)
V_STALL = 63.0         
CL_CR = 0.415                                          
N_L = 4.5                                              # ultimate landing load factor (PLACEHOLDERS)
# landing weight is no longer a fixed fraction: W_L = W_TO - 0.943 W_fuel (Eq. 4.26), see WFUEL_LANDING
N_ENG, ENG_DRY_KG = 2, 3221   
                                         # design cruise C_L                     
# --- Raymer table 15.2 factors
K_INSTALL, K_ALLELSE, K_DOOR = 1.3, 0.17, 1.06


# ------------------------------------------------------------------ groups (kg)
def planform(S):
    b = sqrt(AR * S)
    cr = 2 * S / (b * (1 + TAPER))
    mac = 2 / 3 * cr * (1 + TAPER + TAPER ** 2) / (1 + TAPER)
    return b, mac

def cd0(r):
    """Eq. 5.4: cruise C_D0 at area ratio r = S_w/S_w,0 = W_TO/W_TO,0."""
    c = CD0_REF
    wing = c["wing"] * r ** -0.075
    tails = c["tails"] * r ** 0.39
    fixed = (c["fuselage"] + c["nacelles"] + c["gear"] + c["upsweep"] + c["base"]) / r
    const = c["flap"] + c["wave"]
    return (wing + tails + fixed + const) / K_EXCR

def fuel_fraction(CD0):
    """Eq. 5.3 (k) and Eq. 5.2 (f_f). Returns f_f, k."""
    x = CL_CR / sqrt(pi * AR * E_OSWALD * CD0)
    k = 2 * x / (1 + x ** 2)
    f_f = 1 - exp(-(ALPHA * sqrt(CD0) / k + BETA))
    return f_f, k

def weights(W0, f_f):
    S = W0 * G / W_S
    b, mac = planform(S)
    W = W0 * KG                                        # lb
    Sft = S * FT2
    # wing, 15.25
    wing = (0.0051 * (W * N_Z) ** 0.557 * Sft ** 0.649 * AR ** 0.5 * TC_ROOT ** -0.4
            * (1 + TAPER) ** 0.1 / cos(radians(SWEEP_C4)) * (CSW_RATIO * Sft) ** 0.1)
    Sh, Sv = V_H * S * mac / L_H, V_V * S * b / L_V
    bh = sqrt(A_H * Sh) * FT
    Sh, Sv, Lt = Sh * FT2, Sv * FT2, L_H * FT
    ht = (0.0379 * (1 + FUS_WIDTH * FT / bh) ** -0.25 * W ** 0.639 * N_Z ** 0.10 * Sh ** 0.75 / Lt
          * (0.3 * Lt) ** 0.704 / cos(radians(SWEEP_HT)) * A_H ** 0.166 * 1.3 ** 0.1)
    vt = (0.0026 * (2.0 if T_TAIL else 1.0) ** 0.225 * W ** 0.556 * N_Z ** 0.536 * (L_V * FT) ** -0.5
          * Sv ** 0.5 * (L_V * FT) ** 0.875 / cos(radians(SWEEP_VT)) * A_V ** 0.35 * TC_VT ** -0.5)
    L, D, Sf = FUS_LENGTH * FT, FUS_DEPTH * FT, FUS_SHELL * FT2
    Kws = 0.75 * (1 + 2 * TAPER) / (1 + TAPER) * (b * FT / L) * tan(radians(SWEEP_C4))
    fus = 0.3280 * K_DOOR * (W * N_Z) ** 0.5 * L ** 0.25 * Sf ** 0.302 * (1 + Kws) ** 0.04 * (L / D) ** 0.10
    Wl = W * (1 - WFUEL_LANDING * f_f)
    main = (0.0106 * Wl ** 0.888 * N_L ** 0.25 * (GEAR_L_MAIN * IN) ** 0.4 * N_MW ** 0.321
            * N_MSS ** -0.5 * (V_STALL * FT) ** 0.1)
    nose = 0.032 * Wl ** 0.646 * N_L ** 0.2 * (GEAR_L_NOSE * IN) ** 0.5 * N_NW ** 0.45
    out = dict(wing=wing, fuselage=fus, h_tail=ht, v_tail=vt, landing_gear=main + nose)
    out = {k: v / KG for k, v in out.items()}
    out["engines (installed)"] = K_INSTALL * N_ENG * ENG_DRY_KG
    out["all-else (0.17 W_TO)"] = K_ALLELSE * W0
    return out, S, Sh / FT2, Sv / FT2
def run(tol=0.005):
    W0 = W0_GUESS
    r = r_0
    print(f"{'it':>3} {'W_TO [kg]':>11} {'OEW [kg]':>10} {'S [m2]':>8} {'r':>6} "
          f"{'C_D0':>7} {'k':>6} {'f_f':>6} {'change':>8}")
    f_f = F_FUEL                                       # start from the Class I fuel fraction
    for it in range(1, 51):
        parts, S, Sh, Sv = weights(W0, f_f)            # gear uses W_L from the current f_f (Eq. 4.26)
        oew = sum(parts.values())
        CD0 = cd0(r)                                   # Eq. 5.4
        f_new, k = fuel_fraction(CD0)                  # Eq. 5.3, 5.2
        W_new = (oew + PAYLOAD) / (1 - f_new)          # closure (Eq. 4.39)
        err = abs(W_new - W0) / W0
        err_f = abs(f_new - f_f) / f_f
        print(f"{it:>3} {W0:>11,.0f} {oew:>10,.0f} {S:>8.1f} {r:>6.3f} "
              f"{CD0:>7.5f} {k:>6.3f} {f_f:>6.3f} {err:>8.2%}")
        W0, f_f = W_new, f_new
        r = W0/W0_GUESS                                # S_w / S_w,0 at fixed W/S
        if err < tol and err_f < tol:                  # converged in both W_TO and f_f
            break
    CD0 = cd0(r)
    f_f, k = fuel_fraction(CD0)
    parts, S, Sh, Sv = weights(W0, f_f)
    oew = sum(parts.values())
    LD = k * 0.5 * sqrt(pi * AR * E_OSWALD / CD0)
    print(f"\nW_TO = {W0:,.0f} kg   OEW = {oew:,.0f} kg ({oew / W0:.1%})   "
          f"S = {S:.1f} m2   S_h = {Sh:.1f} m2   S_v = {Sv:.1f} m2")
    print(f"C_D0 = {CD0:.5f}   k = {k:.3f}   (L/D)_CR = {LD:.1f}   "
          f"f_f = {f_f:.4f}  (Class I: {F_FUEL:.4f})   W_fuel = {f_f * W0:,.0f} kg\n")
    for k, v in parts.items():
        print(f"{k:<22}{v:>9,.0f} kg {v / W0:>7.1%}")
if __name__ == "__main__":
    W = 84363.16389 * KG
    S = 169.4002503 * FT2
    w = (0.0051 * (W * N_Z) ** 0.557 * S ** 0.649 * AR ** 0.5 * TC_ROOT ** -0.4 * (1 + TAPER) ** 0.1
         / cos(radians(SWEEP_C4)) * (CSW_RATIO * S) ** 0.1)
    print(f"self-test: wing = {w:,.0f} lb (sheet 20,919)  {'OK' if abs(w / 20919 - 1) < 0.005 else 'MISMATCH'}\n")
    run()
