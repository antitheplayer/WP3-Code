"""
AE2111-I Aircraft System Design - Work Package 3 (WP3.3a)
Class II Aircraft Drag Estimation Script (Interactive / Blank Template v2)

This script performs a complete Class II Component Drag Build-Up for subsonic/transonic 
transport aircraft based on Raymer and AE2111-I specifications.

USERS CAN PROVIDE INPUTS IN THREE WAYS:
1. Edit the USER_INPUTS dictionary at the top of this script.
2. Run interactively from terminal (python class2_drag_estimation-v2.py --interactive).
3. Pass a JSON file with parameters (python class2_drag_estimation-v2.py --json params.json).
"""

import math
import sys
import json
import argparse

# ==============================================================================
# USER INPUT CONFIGURATION (Fill in your aircraft parameters below)
# Set any value to None or a custom float/int.
# ==============================================================================
USER_INPUTS = {
    # --- Flight Condition ---
    'altitude_m': 9448.8,          # Altitude [m] (e.g. 10668 m / 35,000 ft)
    'cruise_mach': 0.85,          # Cruise Mach number [-] (e.g. 0.78)

    # --- Main Wing Geometry ---
    'wing_area': 169.4,           # Wing reference area S_ref [m^2] (e.g. 182.6)
    'wingspan': 38.72,             # Wingspan b [m] (e.g. 40.0)
    'aspect_ratio': 12,         # Aspect ratio AR [-] (e.g. 8.85)
    'sweep_LE_deg': 32.66,         # Leading-edge sweep angle [deg] (e.g. 28.0)
    'sweep_mid_deg': 25.24,        # Mid-chord sweep angle [deg] (e.g. 25.0)
    'wing_mac': 4.132,             # Mean Aerodynamic Chord MAC [m] (e.g. 4.98)
    'wing_tc': 0.121,              # Average thickness-to-chord ratio t/c [-] (e.g. 0.11)
    'wing_x_tmax': 0.30,          # Location of max thickness (x/c)_tmax [-] (default: 0.30)
    'wing_swet': None,            # Wing wetted area [m^2] (None = auto-calculated)

    # --- Fuselage Geometry --
    'fuse_tc_length': 11.216478,
    'fuse_nc_length': 8.4124,
    'fuse_length': 35.28,          # Fuselage overall length [m] (e.g. 37.5)
    'fuse_diameter': 4.01,        # Fuselage max diameter [m] (e.g. 4.0)
    'fuse_upsweep_deg': 6,      # Rear fuselage upsweep angle [deg] (default: 3.0)
    'fuse_base_area': 0.0,        # Fuselage unfaired base area [m^2] (default: 0.0)
    'fuse_swet': None,            # Fuselage wetted area [m^2] (None = auto-calculated)
    
    # --- Horizontal Tail Geometry ---
    'ht_area': 37.62,              # Horizontal tail area [m^2] (e.g. 45.0)
    'ht_mac': 3.26,               # Horizontal tail MAC [m] (e.g. 3.2)
    'ht_tc': 0.12,                # Thickness ratio t/c [-] (default: 0.10)
    'ht_sweep_mid_deg': 26.0,     # Mid-chord sweep [deg] (default: 28.0)
    'ht_swet': None,              # Horizontal tail wetted area [m^2] (None = auto-calc)

    # --- Vertical Tail Geometry ---
    'vt_area': 30.75,              # Vertical tail area [m^2] (e.g. 33.0)
    'vt_mac': 4.72,               # Vertical tail MAC [m] (e.g. 3.8)
    'vt_tc': 0.12,                # Thickness ratio t/c [-] (default: 0.10)
    'vt_sweep_mid_deg': 28.3,     # Mid-chord sweep [deg] (default: 35.0)
    'vt_swet': None,              # Vertical tail wetted area [m^2] (None = auto-calc)

    # --- Nacelles & Engine Geometry ---
    'num_engines': 2,             # Number of engines/nacelles [-] (default: 2)
    'nacelle_length': 6,       # Nacelle length [m] (e.g. 4.2)
    'nacelle_diameter': 2.3,     # Nacelle max diameter [m] (e.g. 2.0)
    'nacelle_swet': None,         # Nacelles wetted area [m^2] (None = auto-calc)

    # --- Interference & Miscellaneous Drag Factors ---
    'IF_wing': 1.00,              # Wing-fuselage interference factor Q (default: 1.00)
    'IF_fuse': 1.00,              # Fuselage interference factor Q (default: 1.00)
    'IF_ht': 1.05,                # Horizontal tail interference factor Q (default: 1.05)
    'IF_vt': 1.05,                # Vertical tail interference factor Q (default: 1.05)
    'IF_nacelle': 1.30,           # Nacelle interference factor Q (default: 1.20)
    'pct_excrescence': 0.03,      # Excrescence & leakage drag fraction (default: 0.03 = 3%) stated in slides n stuff
    
    # --- High-Lift Devices & Off-Design ---
    'flap_chord_ratio': 0.40,     # Flap chord ratio c_f/c [-] (default: 0.25)
    'flapped_wing_area': 72.84211,    # Flapped wing area [m^2] (None = 65% of wing_area) (0 for cruise)
    'flap_TO_deg': 15.0,          # Take-off flap deflection [deg] (default: 15.0)
    'flap_LD_deg': 40.0,          # Landing flap deflection [deg] (default: 35.0)
    'gear_drag_coeff': 0.0180,     # Landing gear drag coefficient delta_CD_gear (default: 0.0180)
    'gear_closed_drag_coeff': 0.0005, # Retracted gear / door leakage drag increment in cruise [-] (e.g. 5 counts)
}


# ==============================================================================
# AERODYNAMIC & ATMOSPHERIC CALCULATIONS
# ==============================================================================

def calculate_atmosphere(alt_m):
    """
    Standard Atmosphere Model up to 11,000m (Troposphere) & Lower Stratosphere.
    Returns density (kg/m^3), dynamic viscosity (kg/m*s), speed of sound (m/s).
    """
    T0 = 288.15
    p0 = 101325.0
    L = 0.0065  # Temperature lapse rate (K/m)
    R = 287.058
    g = 9.80665

    if alt_m <= 11000:
        T = T0 - L * alt_m
        p = p0 * (1 - L * alt_m / T0) ** (g / (R * L))
    else:
        T_11k = T0 - L * 11000
        p_11k = p0 * (1 - L * 11000 / T0) ** (g / (R * L))
        T = T_11k
        p = p_11k * math.exp(-g * (alt_m - 11000) / (R * T))

    rho = p / (R * T)
    a = math.sqrt(1.4 * R * T)
    # Sutherland's law for dynamic viscosity
    mu0 = 1.716e-5
    S_suth = 110.4
    mu = mu0 * ((T / 273.15) ** 1.5) * (273.15 + S_suth) / (T + S_suth)

    return rho, mu, a


def skin_friction_cf(Re, Mach, percent_laminar=0.05, surface_k=0.405e-5, char_length=1.0):

    if Re <= 0:
        return 0.0

    # Cutoff Reynolds number based on surface roughness k
    Re_cutoff = 38.21 * ((char_length / surface_k) ** 1.053)
    Re_effective = min(Re, Re_cutoff)

    # Laminar skin friction
    Cf_lam = 1.328 / math.sqrt(Re)

    # Turbulent skin friction with compressibility correction
    Cf_turb = 0.455 / ((math.log10(Re_effective) ** 2.58) * ((1.0 + 0.144 * Mach**2) ** 0.65))

    # Weighted average
    pct_turb = 1.0 - percent_laminar
    return percent_laminar * Cf_lam + pct_turb * Cf_turb


def form_factor_lifting(tc, x_tmax_c, sweep_mid_rad, Mach):
    """Form Factor (FF) for wings, tails, struts, and pylons."""
    term1 = 1.0 + (0.6 / x_tmax_c) * tc + 100.0 * (tc ** 4)
    term2 = 1.34 * (Mach ** 0.18) * (math.cos(sweep_mid_rad) ** 0.28)
    return term1 * term2


def form_factor_body(fineness_ratio):
    """Form Factor (FF) for fuselage and smooth canopy."""
    f = fineness_ratio
    return 1.0 + (60.0 / (f ** 3)) + (f / 400.0)


def form_factor_nacelle(fineness_ratio):
    """Form Factor (FF) for nacelles and external stores."""
    f = fineness_ratio
    return 1.0 + (0.35 / f)


def oswald_efficiency(AR, sweep_mid_deg,sweep_LE_deg=0.0,):
    """Oswald Efficiency Factor (e)."""
    e = 2.0 / (2.0 - AR + math.sqrt(4.0 + AR**2 * (1.0 + math.tan(math.radians(sweep_mid_deg))**2)))
    return e


def auto_fill_missing_inputs(p):
    """
    Populates any missing (None) values with sensible baseline defaults
    or calculated geometric estimates.
    """
    if p.get('altitude_m') is None: p['altitude_m'] = 10668.0
    if p.get('cruise_mach') is None: p['cruise_mach'] = 0.78
    if p.get('wing_area') is None: p['wing_area'] = 182.6
    if p.get('wingspan') is None: p['wingspan'] = 40.0
    if p.get('aspect_ratio') is None: 
        p['aspect_ratio'] = (p['wingspan']**2) / p['wing_area'] if p.get('wingspan') and p.get('wing_area') else 8.85
    if p.get('sweep_LE_deg') is None: p['sweep_LE_deg'] = 28.0
    if p.get('sweep_mid_deg') is None: p['sweep_mid_deg'] = 25.0
    if p.get('wing_mac') is None: p['wing_mac'] = math.sqrt(p['wing_area'] / p['aspect_ratio'])
    if p.get('wing_tc') is None: p['wing_tc'] = 0.11

    if p.get('fuse_length') is None: p['fuse_length'] = 37.5
    if p.get('fuse_diameter') is None: p['fuse_diameter'] = 4.0

    if p.get('ht_area') is None: p['ht_area'] = 0.25 * p['wing_area']
    if p.get('ht_mac') is None: p['ht_mac'] = 0.65 * p['wing_mac']

    if p.get('vt_area') is None: p['vt_area'] = 0.18 * p['wing_area']
    if p.get('vt_mac') is None: p['vt_mac'] = 0.75 * p['wing_mac']

    if p.get('nacelle_length') is None: p['nacelle_length'] = 4.2
    if p.get('nacelle_diameter') is None: p['nacelle_diameter'] = 2.0

    if p.get('flapped_wing_area') is None: p['flapped_wing_area'] = 0.65 * p['wing_area']

    return p


def perform_class_ii_drag_buildup(inputs):
    """Executes the full Class II Drag Build-Up procedure."""
    p = auto_fill_missing_inputs(inputs.copy())

    alt_m = p['altitude_m']
    Mach = p['cruise_mach']
    rho, mu, a = calculate_atmosphere(alt_m)
    V = Mach * a

    S_ref = p['wing_area']
    b = p['wingspan']
    AR = p['aspect_ratio']
    sweep_LE_deg = p['sweep_LE_deg']
    sweep_mid_deg = p['sweep_mid_deg']
    sweep_mid_rad = math.radians(sweep_mid_deg)

    # 1. Component Geometry & Wetted Areas
    MAC_w = p['wing_mac']
    tc_w = p['wing_tc']
    x_tmax_w = p.get('wing_x_tmax', 0.30)
    S_wet_w = p.get('wing_swet') or (1.07 * 2.0 * S_ref)

    L_fuse = p['fuse_length']
    D_fuse = p['fuse_diameter']
    L1 = p['fuse_nc_length']
    L3 = p['fuse_tc_length']
    L2 = L_fuse - p['fuse_tc_length'] - p['fuse_nc_length']
    f_fuse = L_fuse / D_fuse
    S_wet_fuse = p.get('fuse_swet') or (
        math.pi * D_fuse * L1 * 0.75
        + math.pi * D_fuse * L2
        + math.pi * (D_fuse / 2.0) * math.sqrt((D_fuse / 2.0) ** 2 + L3**2)
    )



    S_ht = p['ht_area']
    MAC_ht = p['ht_mac']
    tc_ht = p['ht_tc']
    sweep_ht_rad = math.radians(p['ht_sweep_mid_deg'])
    S_wet_ht = p.get('ht_swet') or (1.05 * 2.0 * S_ht)

    S_vt = p['vt_area']
    MAC_vt = p['vt_mac']
    tc_vt = p['vt_tc']
    sweep_vt_rad = math.radians(p['vt_sweep_mid_deg'])
    S_wet_vt = p.get('vt_swet') or (1.05 * 2.0 * S_vt)

    N_eng = p['num_engines']
    L_nac = p['nacelle_length']
    D_nac = p['nacelle_diameter']
    f_nac = L_nac / D_nac
    S_wet_nac = p.get('nacelle_swet') or (N_eng * math.pi * D_nac * L_nac)

    # 2. Build-Up Table Calculation
    components = [
        {'name': 'Wing', 'L': MAC_w, 'S_wet': S_wet_w, 'FF': form_factor_lifting(tc_w, x_tmax_w, sweep_mid_rad, Mach), 'IF': p.get('IF_wing', 1.00), 'pct_lam': 0.1},
        {'name': 'Fuselage', 'L': L_fuse, 'S_wet': S_wet_fuse, 'FF': form_factor_body(f_fuse), 'IF': p.get('IF_fuse', 1.00), 'pct_lam':0.05},
        {'name': 'Horizontal Tail', 'L': MAC_ht, 'S_wet': S_wet_ht, 'FF': form_factor_lifting(tc_ht, 0.30, sweep_ht_rad, Mach), 'IF': p.get('IF_ht', 1.05), 'pct_lam':0.1},
        {'name': 'Vertical Tail', 'L': MAC_vt, 'S_wet': S_wet_vt, 'FF': form_factor_lifting(tc_vt, 0.30, sweep_vt_rad, Mach), 'IF': p.get('IF_vt', 1.05), 'pct_lam':0.1},
        {'name': 'Nacelles & Pylons', 'L': L_nac, 'S_wet': S_wet_nac, 'FF': form_factor_nacelle(f_nac), 'IF': p.get('IF_nacelle', 1.20), 'pct_lam':0.05}
    ]

    CD0_clean = 0.0
    table_rows = []

    for comp in components:
        Re = (rho * V * comp['L']) / mu
        Cf = skin_friction_cf(Re, Mach, percent_laminar=comp['pct_lam'], surface_k=0.405e-5, char_length=comp['L'])
        CD0_comp = (Cf * comp['FF'] * comp['IF'] * comp['S_wet']) / S_ref
        CD0_clean += CD0_comp

        table_rows.append({
            'Component': comp['name'],
            'S_wet': comp['S_wet'],
            'Re': Re,
            'Cf': Cf,
            'FF': comp['FF'],
            'IF': comp['IF'],
            'CD0': CD0_comp
        })

    # 3. Excrescence, Upsweep, Base Drag
    CD0_excrescence = p.get('pct_excrescence', 0.03) * CD0_clean
    u_rad = math.radians(p.get('fuse_upsweep_deg', 3.0))
    A_max_fuse = (math.pi / 4.0) * (D_fuse ** 2)
    CD0_upsweep = 3.83 * (u_rad ** 2.5) * (A_max_fuse / S_ref)

    A_base = p.get('fuse_base_area', 0.0)
    CD0_base = (0.139 + 0.419 * ((Mach - 0.16) ** 2)) * (A_base / S_ref) if Mach > 0.16 else 0.0
    
    #  Closed Wheel well Gear Drag
    CD0_gear_closed = p.get('gear_closed_drag_coeff', 0.0)
    
    CD0_total = CD0_clean + CD0_excrescence + CD0_upsweep + CD0_base 

    e_clean = oswald_efficiency(AR, sweep_mid_deg, sweep_LE_deg)
    K_clean = 1.0 / (math.pi * AR * e_clean)
    LD_max_clean = 0.5 * math.sqrt(math.pi * AR * e_clean / CD0_total)

    # 4. Off-Design Configurations
    c_f_c = p.get('flap_chord_ratio', 0.25)
    S_flapped = p.get('flapped_wing_area', 0.65 * S_ref)
    F_flap = 0.0074  # Slotted flap coefficient

    delta_flap_TO = p.get('flap_TO_deg', 15.0)
    CD_flap_TO = F_flap * c_f_c * (S_flapped / S_ref) * (delta_flap_TO - 10.0) if delta_flap_TO > 10 else 0.0

    delta_flap_LD = p.get('flap_LD_deg', 35.0)
    CD_flap_LD = F_flap * c_f_c * (S_flapped / S_ref) * (delta_flap_LD - 10.0) if delta_flap_LD > 10 else 0.0

    CD_gear = p.get('gear_drag_coeff', 0.0180)

    CD0_TO = CD0_total + CD_flap_TO + CD_gear 
    CD0_LD = CD0_total + CD_flap_LD + CD_gear 

    e_TO = e_clean + 0.0026*p.get('flap_TO_deg')
    e_LD = e_clean + 0.0026*p.get('flap_LD_deg')

    return {
        'inputs_used': p,
        'table_rows': table_rows,
        'CD0_clean': CD0_clean,
        'CD0_excrescence': CD0_excrescence,
        'CD0_upsweep': CD0_upsweep,
        'CD0_base': CD0_base,
        'CD0_gear_closed': CD0_gear_closed,
        'CD0_total': CD0_total,
        'e_clean': e_clean,
        'K_clean': K_clean,
        'LD_max_clean': LD_max_clean,
        'CD0_TO': CD0_TO,
        'e_TO': e_TO,
        'CD0_LD': CD0_LD,
        'e_LD': e_LD,
        'CD_flap_TO': CD_flap_TO,
        'CD_flap_LD': CD_flap_LD,
        'CD_gear': CD_gear
        
    }


def prompt_user_interactively():
    """Prompts the user in terminal to enter custom parameter values."""
    print("\n--- INTERACTIVE AIRCRAFT PARAMETER INPUT ---")
    print("Press Enter to keep the [default value] for any parameter.\n")
    p = USER_INPUTS.copy()

    def get_input(prompt_text, key_name, default_val):
        val_str = input(f"{prompt_text} [{default_val}]: ").strip()
        if val_str == "":
            p[key_name] = default_val
        else:
            try:
                p[key_name] = float(val_str)
            except ValueError:
                print(f"Invalid input, using default: {default_val}")
                p[key_name] = default_val

    get_input("Altitude [m]", 'altitude_m', 10668.0)
    get_input("Cruise Mach", 'cruise_mach', 0.78)
    get_input("Wing Reference Area S_ref [m^2]", 'wing_area', 182.6)
    get_input("Wingspan b [m]", 'wingspan', 40.0)
    get_input("Aspect Ratio AR [-]", 'aspect_ratio', 8.85)
    get_input("Leading Edge Sweep [deg]", 'sweep_LE_deg', 28.0)
    get_input("Wing MAC [m]", 'wing_mac', 4.98)
    get_input("Wing Thickness Ratio t/c [-]", 'wing_tc', 0.11)
    get_input("Fuselage Length [m]", 'fuse_length', 37.5)
    get_input("Fuselage Diameter [m]", 'fuse_diameter', 4.0)

    return p


# ==============================================================================
# MAIN EXECUTION ENTRY POINT
# ==============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Class II Aircraft Drag Estimation Script")
    parser.add_argument('--interactive', action='store_true', help="Prompt interactively for input parameters")
    parser.add_argument('--json', type=str, help="Path to JSON file containing parameter inputs")
    args = parser.parse_args()

    current_inputs = USER_INPUTS.copy()

    if args.interactive:
        current_inputs = prompt_user_interactively()
    elif args.json:
        with open(args.json, 'r') as f:
            json_data = json.load(f)
            current_inputs.update(json_data)

    res = perform_class_ii_drag_buildup(current_inputs)
    p = res['inputs_used']

    print("\n" + "="*80)
    print("      CLASS II COMPONENT DRAG BUILD-UP raymer")
    print("="*80)
    print(f"Flight Condition: Mach {p['cruise_mach']} at {p['altitude_m']} m ({p['altitude_m']*3.28084:.0f} ft)")
    print(f"Wing Area S_ref: {p['wing_area']} m^2 | Wingspan: {p['wingspan']} m | Aspect Ratio AR: {p['aspect_ratio']:.2f}\n")

    print(f"{'Component':<20} | {'Swet (m^2)':<10} | {'Reynolds No.':<12} | {'Cf':<8} | {'FF':<7} | {'IF':<5} | {'CD0':<8} | {'Share (%)':<8}")
    print("-" * 90)
    for row in res['table_rows']:
        share = (row['CD0'] / res['CD0_clean']) * 100.0
        print(f"{row['Component']:<20} | {row['S_wet']:<10.1f} | {row['Re']:<12.2e} | {row['Cf']:<8.5f} | {row['FF']:<7.3f} | {row['IF']:<5.2f} | {row['CD0']:<8.5f} | {share:<8.1f}%")

    print("-" * 90)
    print(f"{'Sum Clean CD0':<20} : {res['CD0_clean']:.5f} ({res['CD0_clean']*1e4:.1f} drag counts)")
    print(f"{'Excrescence (+3%)':<20} : {res['CD0_excrescence']:.5f} ({res['CD0_excrescence']*1e4:.1f} drag counts)")
    print(f"{'Fuselage Upsweep':<20} : {res['CD0_upsweep']:.5f} ({res['CD0_upsweep']*1e4:.1f} drag counts)")
    print(f"{'TOTAL CRUISE CD0':<20} : {res['CD0_total']:.5f} ({res['CD0_total']*1e4:.1f} drag counts)")
    print(f"{'Oswald Efficiency (e)':<20} : {res['e_clean']:.3f}")
    print(f"{'Max (L/D)':<20} : {res['LD_max_clean']:.2f}\n")

    print("OFF-DESIGN CONFIGURATIONS:")
    print(f" CD_Flap_TO = {res['CD_flap_TO']:.5f}")
    print(f" CD_Flap_LD = {res['CD_flap_LD']:.5f}")
    print(f" CD_gear = {res['CD_gear']:.5f}")
    print(f"  * Take-Off Configuration  : CD0 = {res['CD0_TO']:.5f} ({res['CD0_TO']*1e4:.1f} counts), e = {res['e_TO']:.3f}")
    print(f"  * Landing Configuration   : CD0 = {res['CD0_LD']:.5f} ({res['CD0_LD']*1e4:.1f} counts), e = {res['e_LD']:.3f}")
    print("="*80 + "\n")
