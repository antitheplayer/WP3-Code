"""
Aircraft Total Drag Polar Plotter (CD vs CL) - Version 3
AE2111-I Aircraft System Design

This script plots total drag polars (CD vs CL) for Cruise, Take-off, and Landing
configurations based on user-input Aspect Ratio (AR), Oswald Efficiency (e),
and Zero-Lift Parasite Drag (CD0). Off-design Oswald efficiency factors (e_TO, e_LD)
are automatically calculated using Obert's flap deflection correction if not explicitly provided.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ==============================================================================
# USER INPUT PARAMETERS
# Leave as None to use defaults, or replace with your calculated numbers.
# ==============================================================================

# Aircraft Geometry & Oswald Efficiency
AR = 12            # Aspect Ratio [-] (e.g. 8.76)
e_clean = 0.6       # Cruise/Clean Oswald efficiency [-] (e.g. 0.790)
e_TO = None          # Take-off Oswald efficiency [-] (None = calculated via flap correction)
e_LD = None          # Landing Oswald efficiency [-] (None = calculated via flap correction)

# Off-Design Flap Settings for Oswald Efficiency Correction (Obert Method)
delta_f_TO = 15.0    # Take-off flap deflection [deg] (default: 15.0°)
delta_f_LD = 40.0    # Landing flap deflection [deg] (default: 40.0°)
engine_mounting = 'wing'  # Engine mounting: 'wing' (0.0046*delta_f) or 'fuselage' (0.0026*delta_f)

# Zero-Lift Parasite Drag Coefficient (CD0) per configuration
CD0_clean = None     # Cruise/Clean CD0 [-] (e.g. 0.01626)
CD0_TO = None        # Take-off CD0 [-] (e.g. 0.04027)
CD0_LD = None        # Landing CD0 [-] (e.g. 0.06432)

# Range of Lift Coefficient (CL) for plotting
CL_min = -2.2
CL_max = 2.2         # Maximum CL range to display on plot [-]
num_points = 300

# File name for saving the output plot
output_filename = "drag_polar_plot.png"


# ==============================================================================
# FALLBACK DEFAULTS & OFF-DESIGN CORRECTIONS
# ==============================================================================
def apply_defaults():
    global AR, e_clean, e_TO, e_LD, CD0_clean, CD0_TO, CD0_LD
    
    if AR is None:
        AR = 8.76
        print("Note: 'AR' was set to None. Using default baseline AR = 8.76")
        
    if e_clean is None:
        e_clean = 0.790
        print("Note: 'e_clean' was set to None. Using default e_clean = 0.790")
        
    # Flap correction factor according to Obert:
    # delta_e = 0.0026 * delta_f (fuselage mounted) or 0.0046 * delta_f (wing mounted)
    k_flap = 0.0046 if engine_mounting.lower() == 'wing' else 0.0026

    if e_TO is None:
        e_TO = e_clean + k_flap * delta_f_TO
        print(f"Calculated Take-off e_TO = {e_TO:.4f} (e_clean {e_clean:.3f} + {k_flap} * {delta_f_TO}° flap)")
        
    if e_LD is None:
        e_LD = e_clean + k_flap * delta_f_LD
        print(f"Calculated Landing e_LD  = {e_LD:.4f} (e_clean {e_clean:.3f} + {k_flap} * {delta_f_LD}° flap)")

    if CD0_clean is None:
        CD0_clean = 0.01626
        print("Note: 'CD0_clean' was set to None. Using default CD0_clean = 0.01626")
        
    if CD0_TO is None:
        CD0_TO = 0.04027
        print("Note: 'CD0_TO' was set to None. Using default CD0_TO = 0.04027")
        
    if CD0_LD is None:
        CD0_LD = 0.06432
        print("Note: 'CD0_LD' was set to None. Using default CD0_LD = 0.06432")


# ==============================================================================
# PLOTTING ROUTINE
# ==============================================================================
def plot_drag_polars():
    apply_defaults()

    # Vector of CL values
    CL = np.linspace(CL_min, CL_max, num_points)

    # Compute induced drag factor K = 1 / (pi * AR * e)
    K_clean = 1.0 / (np.pi * AR * e_clean)
    K_TO = 1.0 / (np.pi * AR * e_TO)
    K_LD = 1.0 / (np.pi * AR * e_LD)

    # Compute total drag CD = CD0 + K * CL^2
    CD_clean = CD0_clean + K_clean * (CL ** 2)
    CD_TO = CD0_TO + K_TO * (CL ** 2)
    CD_LD = CD0_LD + K_LD * (CL ** 2)

    # Create Plot
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)

    # Plot configuration curves
    ax.plot(CD_clean, CL, label=f"Cruise (Clean): $C_{{D0}}={CD0_clean:.4f}, e={e_clean:.3f}$", 
            color='#1f77b4', linestyle='-', linewidth=2)
    ax.plot(CD_TO, CL, label=f"Take-off: $C_{{D0}}={CD0_TO:.4f}, e={e_TO:.3f}$", 
            color='#ff7f0e', linestyle='--', linewidth=2)
    ax.plot(CD_LD, CL, label=f"Landing: $C_{{D0}}={CD0_LD:.4f}, e={e_LD:.3f}$", 
            color='#2ca02c', linestyle='-.', linewidth=2)

    # Highlight maximum lift-to-drag ratio (L/D)_max on cruise curve
    CL_opt = np.sqrt(CD0_clean / K_clean)
    CD_opt = 2.0 * CD0_clean
    LD_max = 0.5 * np.sqrt(np.pi * AR * e_clean / CD0_clean)

    # Draw dot at (L/D)_max point (without text/arrow annotation)
    if CL_opt <= CL_max:
        ax.scatter([CD_opt], [CL_opt], color='#1f77b4', s=60, zorder=5)

    # Print (L/D)_max results directly to console
    print("\n" + "="*50)
    print("       MAXIMUM LIFT-TO-DRAG RATIO (L/D)_max RESULTS")
    print("="*50)
    print(f"(L/D)_max           : {LD_max:.2f}")
    print(f"Optimal CL (CL_opt) : {CL_opt:.4f}")
    print(f"Optimal CD (CD_opt) : {CD_opt:.4f} ({CD_opt*10000:.1f} counts)")
    print("="*50 + "\n")

    # Axes and Formatting
    ax.set_title(f"Aircraft Total Drag Polars ($AR = {AR:.2f}$)", fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel("Total Drag Coefficient $C_D$ [-]", fontsize=11)
    ax.set_ylabel("Lift Coefficient $C_L$ [-]", fontsize=11)
    ax.set_ylim(CL_min, CL_max)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')
    plt.close()

    print(f"Plot successfully generated and saved to '{output_filename}'")


if __name__ == '__main__':
    plot_drag_polars()
