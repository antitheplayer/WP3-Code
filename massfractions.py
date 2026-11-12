import numpy as np
import pandas as pd

aircraft_data = {
    "Airbus A330-300": {
        "M_MTO": 233000.0, "M_ZF": 175000.0, "M_OE": 124500.0,
        "b": 60.30, "S_w": 361.6, "Lambda_half": 27.5, "t_root": 1.38,
        "n_ult": 3.75, "K_w": 0.95,
        "l_t": 28.50, "b_f": 5.64, "h_f": 5.64, "S_gs": 927.0, "V_D": 195.0, "K_f": 1.08,
        "S_h": 72.9, "Lambda_h": 30.0, "K_h": 1.1,
        "S_v": 47.6, "Lambda_v": 38.0, "K_v": 1.0,
        "N_eng": 2, "M_eng_dry": 5094.0, "T_TO": 631600.0,
        "K_lg": 1.00,
        "A_mlg": 18.1, "B_mlg": 0.131, "C_mlg": 0.019, "D_mlg": 2.23e-5,
        "A_nlg": 9.1,  "B_nlg": 0.082, "C_nlg": 0.0,   "D_nlg": 2.97e-6,
        "K_fe": 0.115
    },
    "Boeing BBJ": {
        "M_MTO": 77564.0, "M_ZF": 57153.0, "M_OE": 42899.0,
        "b": 35.79, "S_w": 124.6, "Lambda_half": 23.0, "t_root": 0.72,
        "n_ult": 3.75, "K_w": 0.95,
        "l_t": 15.20, "b_f": 3.76, "h_f": 4.01, "S_gs": 322.2, "V_D": 200.0, "K_f": 1.08,
        "S_h": 32.8, "Lambda_h": 30.0, "K_h": 1.1,
        "S_v": 22.7, "Lambda_v": 35.0, "K_v": 1.0,
        "N_eng": 2, "M_eng_dry": 2366.0, "T_TO": 234000.0,
        "K_lg": 1.00,
        "A_mlg": 18.1, "B_mlg": 0.131, "C_mlg": 0.019, "D_mlg": 2.23e-5,
        "A_nlg": 9.1,  "B_nlg": 0.082, "C_nlg": 0.0,   "D_nlg": 2.97e-6,
        "K_fe": 0.145
    },
    "Douglas DC-8-72": {
        "M_MTO": 158757.0, "M_ZF": 111130.0, "M_OE": 68810.0,
        "b": 45.23, "S_w": 271.9, "Lambda_half": 28.0, "t_root": 1.05,
        "n_ult": 3.75, "K_w": 0.90,
        "l_t": 23.50, "b_f": 3.73, "h_f": 4.14, "S_gs": 565.0, "V_D": 205.0, "K_f": 1.08,
        "S_h": 54.2, "Lambda_h": 32.0, "K_h": 1.1,
        "S_v": 40.6, "Lambda_v": 35.0, "K_v": 1.0,
        "N_eng": 4, "M_eng_dry": 2100.0, "T_TO": 391400.0,
        "K_lg": 1.00,
        "A_mlg": 18.1, "B_mlg": 0.131, "C_mlg": 0.019, "D_mlg": 2.23e-5,
        "A_nlg": 9.1,  "B_nlg": 0.082, "C_nlg": 0.0,   "D_nlg": 2.97e-6,
        "K_fe": 0.105
    },
    "Boeing 707-320B": {
        "M_MTO": 151318.0, "M_ZF": 104326.0, "M_OE": 66406.0,
        "b": 44.42, "S_w": 283.4, "Lambda_half": 32.5, "t_root": 1.08,
        "n_ult": 3.75, "K_w": 0.90,
        "l_t": 20.10, "b_f": 3.76, "h_f": 4.32, "S_gs": 468.0, "V_D": 205.0, "K_f": 1.08,
        "S_h": 51.1, "Lambda_h": 35.0, "K_h": 1.1,
        "S_v": 37.7, "Lambda_v": 35.0, "K_v": 1.0,
        "N_eng": 4, "M_eng_dry": 2059.0, "T_TO": 320300.0,
        "K_lg": 1.00,
        "A_mlg": 18.1, "B_mlg": 0.131, "C_mlg": 0.019, "D_mlg": 2.23e-5,
        "A_nlg": 9.1,  "B_nlg": 0.082, "C_nlg": 0.0,   "D_nlg": 2.97e-6,
        "K_fe": 0.105
    },
    "Dassault Falcon 8X": {
        "M_MTO": 33112.0, "M_ZF": 18597.0, "M_OE": 16692.0,
        "b": 26.29, "S_w": 70.7, "Lambda_half": 28.5, "t_root": 0.52,
        "n_ult": 3.75, "K_w": 1.00,
        "l_t": 10.50, "b_f": 2.52, "h_f": 2.52, "S_gs": 156.0, "V_D": 190.0, "K_f": 1.15,
        "S_h": 16.5, "Lambda_h": 32.0, "K_h": 1.1,
        "S_v": 12.3, "Lambda_v": 38.0, "K_v": 1.0,
        "N_eng": 3, "M_eng_dry": 640.0, "T_TO": 89680.0,
        "K_lg": 1.00,
        "A_mlg": 15.0, "B_mlg": 0.033, "C_mlg": 0.021, "D_mlg": 0.0,
        "A_nlg": 5.4,  "B_nlg": 0.049, "C_nlg": 0.0,   "D_nlg": 0.0,
        "K_fe": 0.145
    }
}

def calc_wing_mass(d):
    """Torenbeek Eq. 8-12: Wing group mass [kg]"""
    lam = np.radians(d["Lambda_half"])
    b_s = d["b"] / np.cos(lam)
    return (
        0.00667 * d["K_w"] * d["M_ZF"]
        * (b_s ** 0.75) * (1.0 + np.sqrt(1.905 / b_s))
        * (d["n_ult"] ** 0.55)
        * ((d["b"] * d["S_w"]) / (d["t_root"] * d["M_ZF"] * np.cos(lam))) ** 0.30
    )

def calc_fuselage_mass(d):
    """Torenbeek Eq. 8-17: Pressurized fuselage shell mass [kg]"""
    return (
        0.23 * d["K_f"]
        * np.sqrt((d["V_D"] * d["l_t"]) / (d["b_f"] + d["h_f"]))
        * (d["S_gs"] ** 1.2)
    )

def calc_tail_mass(d):
    """Torenbeek Eq. 8-15 & 8-16: Horizontal + Vertical tail mass [kg]"""
    lam_h = np.radians(d["Lambda_h"])
    lam_v = np.radians(d["Lambda_v"])
    m_h = d["K_h"] * d["S_h"] * (62.0 * (d["S_h"]**0.2 * d["V_D"]) / (1000.0 * np.sqrt(np.cos(lam_h))) - 2.5)
    m_v = d["K_v"] * d["S_v"] * (62.0 * (d["S_v"]**0.2 * d["V_D"]) / (1000.0 * np.sqrt(np.cos(lam_v))) - 2.5)
    return m_h + m_v

def calc_propulsion_mass(d):
    """Torenbeek Eq. 8-24: Installed propulsion group mass [kg]"""
    return 1.35 * d["N_eng"] * d["M_eng_dry"]

def calc_nacelle_mass(d):
    """Torenbeek turbofan nacelle group mass (~0.065 lb/lbf) [kg]"""
    return (0.065 / 9.80665) * d["T_TO"]

def calc_landing_gear_mass(d):
    """Torenbeek Eq. 8-20: Main + Nose landing gear mass [kg]"""
    M = d["M_MTO"]
    m_mlg = d["K_lg"] * (d["A_mlg"] + d["B_mlg"] * (M**0.75) + d["C_mlg"] * M + d["D_mlg"] * (M**1.5))
    m_nlg = d["K_lg"] * (d["A_nlg"] + d["B_nlg"] * (M**0.75) + d["C_nlg"] * M + d["D_nlg"] * (M**1.5))
    return m_mlg + m_nlg

def calc_fixed_equipment_mass(d):
    """Torenbeek statistical fixed equipment mass [kg]"""
    return d["K_fe"] * d["M_MTO"]

masses = {}
mtow_fractions = {}
oew_fractions = {}

for name, d in aircraft_data.items():
    m_wing = calc_wing_mass(d)
    m_fuse = calc_fuselage_mass(d)
    m_tail = calc_tail_mass(d)
    m_prop = calc_propulsion_mass(d)
    m_nac  = calc_nacelle_mass(d)
    m_lg   = calc_landing_gear_mass(d)
    m_fe   = calc_fixed_equipment_mass(d)
    
    m_sum = m_wing + m_fuse + m_tail + m_prop + m_nac + m_lg + m_fe
    m_unacc = d["M_OE"] - m_sum
    
    group_dict = {
        "(operating) empty mass": d["M_OE"],
        "wing group": m_wing,
        "fuselage group": m_fuse,
        "tail group": m_tail,
        "propulsion group": m_prop,
        "nacelle group": m_nac,
        "landing gear group": m_lg,
        "fixed equipment group": m_fe,
        "unaccounted": m_unacc
    }
    masses[name] = group_dict
    mtow_fractions[name] = {k: v / d["M_MTO"] for k, v in group_dict.items()}
    oew_fractions[name] = {
        "(operating) empty mass": d["M_OE"] / d["M_MTO"],
        **{k: v / d["M_OE"] for k, v in group_dict.items() if k != "(operating) empty mass"}
    }

df_masses = pd.DataFrame(masses)
df_mtow = pd.DataFrame(mtow_fractions)
df_oew = pd.DataFrame(oew_fractions)

print("--- COMPONENT MASSES (kg) ---")
print(df_masses.round(1))
print("\n--- MASS FRACTIONS (RELATIVE TO MTOW) ---")
print(df_mtow.round(4))
print("\n--- MASS FRACTIONS (RELATIVE TO OEW) ---")
print(df_oew.round(4))