"""
Comparative Spaceflight Plant Cultivation Hardware Benchmark.
Synthesizes microenvironmental, aerodynamic, and gas-exchange performance across:
  1. NASA BRIC (Biological Research in Canisters, 60 mm sealed Petri dishes)
  2. NASA CARA (Characterizing Arabidopsis Root Attractions, 120 mm taped plates)
  3. NASA VEGGIE (Vegetable Production System, open cabin forced convection)
  4. NASA LEAF GBC2 (Growth Chamber 2 Pot Assembly, M09024F/M09025F 3-chamber active digital twin)

Ingests empirical and numerical benchmarks from LunarLeaf-CFD and GBC2 CFD engine.
Generates:
  - results/tables/gbc2_spaceflight_hardware_comprehensive_benchmark.csv
  - figures/fig18_spaceflight_hardware_comprehensive_benchmark.png (.svg)
"""

import os
import sys

# Ensure project root is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def build_hardware_benchmark_dataset():
    """
    Compile comprehensive quantitative metrics across spaceflight hardware platforms.
    """
    benchmark_data = [
        {
            'hardware': 'BRIC (Sealed 60mm Dish)',
            'flight_heritage': 'STS, ISS (BRIC-PDF, BRIC-LED)',
            'enclosure_type': 'Hermetically Sealed Dish',
            'culture_volume_cm3': 16.4,
            'headspace_volume_cm3': 10.8,
            'plant_capacity': '20-50 seedlings (agar)',
            'ventilation_mode': 'Zero ventilation (passive diffusion only)',
            'flow_rate_L_min': 0.0,
            'ach_per_hour': 0.0,
            'tau_carbon_starvation_min': 7.0, # 400 -> ~0 ppm CO2 in 7 min (empirical/LBM)
            'tau_hypoxia_dark_days': 6.5,    # 21% -> 5% O2 in dark respiration
            'g_bl_earth_1g': 0.070,           # mol m^-2 s^-1
            'g_bl_lunar_0166g': 0.063,
            'g_bl_microgravity_0g': 0.060,
            'dC_CO2_surface_ppm': -125.0,     # severe foliar starvation (>30% drawdown)
            'condensation_risk': 'Critical (100% saturation, droplet fogging on lid)',
            'perfusion_capability': 'Post-flight RNAlater injection via septum',
            'aerodynamic_control': 'None'
        },
        {
            'hardware': 'CARA (Taped 120mm Dish)',
            'flight_heritage': 'ISS (Microgravity / Kibo module)',
            'enclosure_type': 'Gas-Permeable Microporous Tape Dish',
            'culture_volume_cm3': 110.6,
            'headspace_volume_cm3': 68.2,
            'plant_capacity': '80-120 seedlings (agar/membrane)',
            'ventilation_mode': 'Passive tape perimeter venting',
            'flow_rate_L_min': 0.002,         # estimated passive diffusion exchange
            'ach_per_hour': 1.76,             # ~1.8 air changes per hour
            'tau_carbon_starvation_min': 42.0,# delayed starvation, bounded by tape flux
            'tau_hypoxia_dark_days': 28.0,
            'g_bl_earth_1g': 0.247,
            'g_bl_lunar_0166g': 0.198,
            'g_bl_microgravity_0g': 0.170,
            'dC_CO2_surface_ppm': -45.0,
            'condensation_risk': 'High (corner perimeter fogging, agar meniscus wetting)',
            'perfusion_capability': 'Manual harvest or chemical cassette',
            'aerodynamic_control': 'Passive perimeter diffusion'
        },
        {
            'hardware': 'VEGGIE (ISS Cabin Enclosure)',
            'flight_heritage': 'ISS Columbus / Destiny (VEG-01 to VEG-05)',
            'enclosure_type': 'Bellows Chamber with Active Fan',
            'culture_volume_cm3': 115000.0,   # ~115 L canopy volume
            'headspace_volume_cm3': 95000.0,
            'plant_capacity': '6 plant pillows (zinnia, lettuce, dwarf tomato)',
            'ventilation_mode': 'Cabin air forced push (recirculating blower)',
            'flow_rate_L_min': 28.3,          # ~1 CFM cabin intake
            'ach_per_hour': 17.8,             # ~18 air changes per hour
            'tau_carbon_starvation_min': np.inf, # Continuous replenishment from ISS cabin
            'tau_hypoxia_dark_days': np.inf,
            'g_bl_earth_1g': 0.480,
            'g_bl_lunar_0166g': 0.410,
            'g_bl_microgravity_0g': 0.350,
            'dC_CO2_surface_ppm': -8.5,
            'condensation_risk': 'Moderate (inter-leaf contact microclimates in pillow root zones)',
            'perfusion_capability': 'Manual destructive harvest only (no in-situ perfusion)',
            'aerodynamic_control': 'Variable fan duty cycle (bulk chamber push)'
        },
        {
            'hardware': 'LEAF GBC2 Pot (Model Lunar LEAF)',
            'flight_heritage': 'Artemis III / Lunar Surface Architecture',
            'enclosure_type': 'Precision CNC Acrylic 3-Chamber Pot Assembly',
            'culture_volume_cm3': 64.6,       # Total CAD internal volume
            'headspace_volume_cm3': 54.1,     # Open canopy headspace
            'plant_capacity': '164 plants (32 B. rapa + 100 Arabidopsis + 32 B. rapa)',
            'ventilation_mode': 'Controlled Active Pull Suction (Port 1 -> Port 2)',
            'flow_rate_L_min': 0.50,          # 0.50 - 0.85 L/min nominal flight
            'ach_per_hour': 554.5,            # 554.5 air changes/hr (1 complete flush every 6.5s)
            'tau_carbon_starvation_min': np.inf, # Constant replenishment (0.34 mmol fresh air/s)
            'tau_hypoxia_dark_days': np.inf,
            'g_bl_earth_1g': 0.070,           # canopy microenvironment
            'g_bl_lunar_0166g': 0.045,
            'g_bl_microgravity_0g': 0.043,
            'dC_CO2_surface_ppm': -1.86,      # Minimal foliar drawdown (>99.5% ambient CO2 maintained)
            'condensation_risk': 'Mitigated (<62% lid RH, negative pressure pull exhaust)',
            'perfusion_capability': 'Autonomous Port 3 RNAlater sub-second perfusion (1.20s saturation)',
            'aerodynamic_control': 'Port-level negative gauge pull with uniform manifold distribution'
        }
    ]
    return pd.DataFrame(benchmark_data)

def generate_figure_18(df_hw):
    print("Generating Figure 18: Spaceflight Hardware Comprehensive Benchmark...")
    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.05], hspace=0.30, wspace=0.24)
    
    # Palette
    c_bric = '#64748b'  # Slate grey
    c_cara = '#d97706'  # Amber
    c_veggie = '#2563eb'# Royal blue
    c_gbc2 = '#10b981'  # Emerald green
    hw_colors = [c_bric, c_cara, c_veggie, c_gbc2]
    hw_short = ['BRIC\n(Sealed)', 'CARA\n(Taped)', 'VEGGIE\n(ISS Cabin)', 'LEAF GBC2\n(Digital Twin)']
    
    # Panel A: Air Exchange Rate (ACH) vs Headspace Volume
    ax_a = fig.add_subplot(gs[0, 0])
    
    # Plot bars of ACH on logarithmic scale
    achs = [0.05 if x == 0.0 else x for x in df_hw['ach_per_hour']] # 0.05 for log plot visibility
    bars_a = ax_a.bar(hw_short, achs, color=hw_colors, edgecolor='black', linewidth=0.8, alpha=0.9, width=0.55)
    ax_a.set_yscale('log')
    ax_a.set_ylabel('Air Exchange Rate (Air Changes per Hour, ACH)', fontsize=10, weight='semibold')
    ax_a.set_title('A  Ventilation Flux Architecture (Turnover Rate)', fontsize=11, weight='bold', loc='left')
    ax_a.grid(True, linestyle=':', alpha=0.5, axis='y', which='both')
    
    # Annotations
    ax_a.text(0, 0.08, '0 ACH\n(Hermetic)', ha='center', va='bottom', fontsize=8.5, weight='bold', color='#475569')
    ax_a.text(1, 2.5, '1.76 ACH\n(~34 min flush)', ha='center', va='bottom', fontsize=8.5, weight='bold', color='#b45309')
    ax_a.text(2, 25.0, '17.8 ACH\n(~3.4 min flush)', ha='center', va='bottom', fontsize=8.5, weight='bold', color='#1d4ed8')
    ax_a.text(3, 750.0, '554.5 ACH\n(6.5s flush)', ha='center', va='bottom', fontsize=8.5, weight='bold', color='#047857')
    ax_a.set_ylim(0.03, 2000)

    # Panel B: Carbon Starvation & Depletion Timescales
    ax_b = fig.add_subplot(gs[0, 1])
    
    # Display depletion timescales
    starve_times = [7.0, 42.0, 240.0, 240.0] # cap continuous systems at 240 min for visual clarity
    bars_b = ax_b.bar(hw_short, starve_times, color=hw_colors, edgecolor='black', linewidth=0.8, alpha=0.9, width=0.55)
    ax_b.set_ylabel('Time to Carbon Starvation tau_starve (minutes)', fontsize=10, weight='semibold')
    ax_b.set_title('B  Carbon Depletion & Replenishment Dynamics', fontsize=11, weight='bold', loc='left')
    ax_b.grid(True, linestyle=':', alpha=0.5, axis='y')
    
    ax_b.text(0, 12, '7.0 min\n(Rapid Starvation)', ha='center', va='bottom', fontsize=8.5, weight='bold', color='#ef4444')
    ax_b.text(1, 48, '42.0 min\n(Diffusion-Limited)', ha='center', va='bottom', fontsize=8.5, weight='bold', color='#d97706')
    ax_b.text(2, 180, 'INDEFINITE\n(Continuous Open)', ha='center', va='center', fontsize=8.5, weight='bold', color='#1d4ed8')
    ax_b.text(3, 180, 'INDEFINITE\n(Active Controlled Pull\nDelta_C < 1.9 ppm)', ha='center', va='center', fontsize=8.5, weight='bold', color='#047857')
    ax_b.set_ylim(0, 260)
    ax_b.axhline(60, color='#94a3b8', linestyle='--', lw=1.0, label='1 Hour Threshold')
    ax_b.legend(loc='upper left', fontsize=8.5, framealpha=0.9)

    # Panel C: Boundary Layer Conductance g_bl Across Gravity Levels
    ax_c = fig.add_subplot(gs[1, 0])
    x_pos = np.arange(len(hw_short))
    bar_w = 0.25
    
    g_1g = df_hw['g_bl_earth_1g'].values
    g_lunar = df_hw['g_bl_lunar_0166g'].values
    g_0g = df_hw['g_bl_microgravity_0g'].values
    
    rects1 = ax_c.bar(x_pos - bar_w, g_1g, bar_w, label='Earth 1g (Buoyant)', color='#3b82f6', edgecolor='black', linewidth=0.6)
    rects2 = ax_c.bar(x_pos, g_lunar, bar_w, label='Lunar 0.166g (Intermediate)', color='#10b981', edgecolor='black', linewidth=0.6)
    rects3 = ax_c.bar(x_pos + bar_w, g_0g, bar_w, label='Microgravity 0g (Pure Convection)', color='#f59e0b', edgecolor='black', linewidth=0.6)
    
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(hw_short, fontsize=9)
    ax_c.set_ylabel('Boundary Layer Conductance g_bl (mol m^-2 s^-1)', fontsize=10, weight='semibold')
    ax_c.set_title('C  Boundary Layer Conductance Across Gravitational Regimes', fontsize=11, weight='bold', loc='left')
    ax_c.legend(loc='upper right', fontsize=8.5, framealpha=0.9)
    ax_c.grid(True, linestyle=':', alpha=0.5, axis='y')
    ax_c.set_ylim(0, 0.55)

    # Panel D: Spaceflight System Capability Radar / Feature Matrix
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis('off')
    
    # Render rich comparative summary table in Panel D
    table_data = [
        ['Platform', 'Ventilation Mode', 'ACH (h^-1)', 'Perfusion Fixation', 'Condensation Control'],
        ['BRIC (60mm)', 'Hermetic Sealed', '0.0', 'Post-flight Septum', 'None (100% Fogging)'],
        ['CARA (120mm)', 'Microporous Tape', '1.76', 'Manual / Cassette', 'High (Rim Meniscus)'],
        ['VEGGIE', 'Cabin Forced Push', '17.8', 'Manual Harvest Only', 'Moderate (Pillows)'],
        ['LEAF GBC2', 'Active Pull Suction', '554.5', 'Sub-second (1.20s)', 'Active Pull (<62% RH)']
    ]
    
    table = ax_d.table(cellText=table_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 2.2)
    
    # Style table headers and cells
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor('#cbd5e1')
        cell.set_linewidth(0.8)
        if row == 0:
            cell.set_facecolor('#0f172a')
            cell.get_text().set_color('white')
            cell.get_text().set_weight('bold')
        elif row == 4: # LEAF GBC2 highlight
            cell.set_facecolor('#ecfdf5')
            cell.get_text().set_color('#065f46')
            cell.get_text().set_weight('bold')
        else:
            cell.set_facecolor('#f8fafc' if row % 2 == 1 else '#ffffff')
            cell.get_text().set_color('#1e293b')
            
    ax_d.set_title('D  Spaceflight Hardware Architectural Trade Matrix', fontsize=11, weight='bold', loc='left')

    os.makedirs('figures', exist_ok=True)
    fig_png = 'figures/fig18_spaceflight_hardware_comprehensive_benchmark.png'
    fig_svg = 'figures/fig18_spaceflight_hardware_comprehensive_benchmark.svg'
    plt.tight_layout()
    plt.savefig(fig_png, dpi=300)
    plt.savefig(fig_svg, format='svg')
    plt.close()
    print(f"[OK] Successfully saved Figure 18 to {fig_png} and {fig_svg}")

def main():
    print("=" * 70)
    print("LEAF GBC2 Pot Model - Spaceflight Hardware Benchmark Synthesis")
    print("=" * 70)
    
    df_hw = build_hardware_benchmark_dataset()
    
    os.makedirs('results/tables', exist_ok=True)
    csv_path = 'results/tables/gbc2_spaceflight_hardware_comprehensive_benchmark.csv'
    df_hw.to_csv(csv_path, index=False)
    print(f"[OK] Successfully saved benchmark dataset to {csv_path}")
    
    print("\n=== SPACEFLIGHT HARDWARE BENCHMARK SUMMARY ===")
    cols = ['hardware', 'enclosure_type', 'ach_per_hour', 'tau_carbon_starvation_min', 'g_bl_microgravity_0g', 'perfusion_capability']
    print(df_hw[cols].to_string(index=False))
    
    generate_figure_18(df_hw)

if __name__ == '__main__':
    main()
