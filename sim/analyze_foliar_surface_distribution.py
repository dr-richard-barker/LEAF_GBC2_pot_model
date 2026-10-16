"""
Foliar Surface Distribution & Boundary Layer Hotspot Analytics for LEAF GBC1 Pot.
Analyzes microenvironmental heterogeneity across all 164 individual plant locations:
  - 32 Brassica rapa in Chamber 1 (windward tall canopy)
  - 100 Arabidopsis thaliana in Chamber 2 (central prostrate rosette lawn)
  - 32 Brassica rapa in Chamber 3 (leeward tall canopy)

Computes spatial percentiles (p10, p50, p90), coefficient of variation (CoV),
dynamic range (p90/p10), and stagnation hotspot fractions across flight flow regimes.
Outputs:
  - results/tables/gbc1_164plants_foliar_telemetry.csv
  - results/tables/gbc1_foliar_surface_distribution.csv
  - figures/fig17_foliar_surface_distribution_and_hotspots.png (.svg)
"""

import os
import sys

# Ensure project root is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sim.gbc1_cfd_engine import GBC1CFDSimulator

def generate_164_plant_coordinates():
    """
    Generate coordinates for all 164 plants across the three chambers.
    Coordinates match the 3D digital twin geometry in viewer_3d.html.
    Returns list of dicts with plant ID, species, chamber, x, y, z (in meters and mm).
    """
    plants = []
    pid = 1
    
    # Chamber 1: 32 Brassica rapa (8 x 4 lattice in X in [-8, 18], Z in [-28.5, -12.0])
    for ix in range(8):
        for iz in range(4):
            px_mm = -8.0 + (ix / 7.0) * 26.0
            pz_mm = -28.5 + (iz / 3.0) * 16.5
            py_mm = 18.0 # mid-canopy height
            plants.append({
                'plant_id': pid,
                'species': 'Brassica rapa',
                'chamber': 1,
                'x_mm': px_mm,
                'y_mm': py_mm,
                'z_mm': pz_mm,
                'x_m': px_mm * 1e-3,
                'y_m': py_mm * 1e-3,
                'z_m': pz_mm * 1e-3,
                'L_char_m': 0.020 # 20 mm characteristic leaf dimension
            })
            pid += 1
            
    # Chamber 2: 100 Arabidopsis thaliana (10 x 10 dense lawn in X in [-8, 18], Z in [-6.5, +6.5])
    for ix in range(10):
        for iz in range(10):
            px_mm = -8.0 + (ix / 9.0) * 26.0
            pz_mm = -6.5 + (iz / 9.0) * 13.0
            py_mm = 16.5 # mid-rosette height
            plants.append({
                'plant_id': pid,
                'species': 'Arabidopsis thaliana',
                'chamber': 2,
                'x_mm': px_mm,
                'y_mm': py_mm,
                'z_mm': pz_mm,
                'x_m': px_mm * 1e-3,
                'y_m': py_mm * 1e-3,
                'z_m': pz_mm * 1e-3,
                'L_char_m': 0.008 # 8 mm characteristic leaf dimension
            })
            pid += 1

    # Chamber 3: 32 Brassica rapa (8 x 4 lattice in X in [-8, 18], Z in [+12.0, +28.5])
    for ix in range(8):
        for iz in range(4):
            px_mm = -8.0 + (ix / 7.0) * 26.0
            pz_mm = 12.0 + (iz / 3.0) * 16.5
            py_mm = 18.0
            plants.append({
                'plant_id': pid,
                'species': 'Brassica rapa',
                'chamber': 3,
                'x_mm': px_mm,
                'y_mm': py_mm,
                'z_mm': pz_mm,
                'x_m': px_mm * 1e-3,
                'y_m': py_mm * 1e-3,
                'z_m': pz_mm * 1e-3,
                'L_char_m': 0.020
            })
            pid += 1

    return plants

def sample_plant_microenvironment(sim, fields, plants, gravity_g):
    """
    Sample scalar fields at each plant location using local 3x3x3 voxel neighborhood.
    Computes local Churchill boundary layer conductance and stagnation metrics.
    """
    inv_T = np.linalg.inv(sim.transform)
    Nx, Ny, Nz = sim.grid.shape
    u_mag = fields['u_mag']
    c_CO2 = fields['c_CO2']
    c_H2O = fields['c_H2O']
    
    plant_data = []
    
    for p in plants:
        # Map to voxel indices
        pv = np.array([p['x_m'], p['y_m'], p['z_m'], 1.0])
        vox = inv_T @ pv
        vx = int(np.clip(round(vox[0]), 0, Nx - 1))
        vy = int(np.clip(round(vox[1]), 0, Ny - 1))
        vz = int(np.clip(round(vox[2]), 0, Nz - 1))
        
        # 3x3x3 voxel neighborhood around foliage
        x_min, x_max = max(0, vx - 1), min(Nx, vx + 2)
        y_min, y_max = max(0, vy - 1), min(Ny, vy + 2)
        z_min, z_max = max(0, vz - 1), min(Nz, vz + 2)
        
        sub_u = u_mag[x_min:x_max, y_min:y_max, z_min:z_max]
        sub_co2 = c_CO2[x_min:x_max, y_min:y_max, z_min:z_max]
        sub_h2o = c_H2O[x_min:x_max, y_min:y_max, z_min:z_max]
        
        u_local_m_s = float(np.mean(sub_u))
        u_local_cm_s = u_local_m_s * 100.0
        dC_co2_ppm = float(abs(np.mean(sub_co2)))
        dH_h2o_mmol = float(np.mean(sub_h2o))
        
        # Boundary layer calculation
        L = p['L_char_m']
        Gr = (gravity_g * sim.beta_T * sim.dT_leaf * (L ** 3)) / (sim.nu ** 2) if gravity_g > 0 else 0.0
        delta_nat_m = 3.5 * L * (Gr ** -0.25) if Gr > 1.0 else 0.015
        
        u_eff = max(u_local_m_s, 1e-5)
        delta_forced_m = 4.92 * np.sqrt(sim.nu * L / u_eff)
        
        delta_eff_m = ( (1.0 / (delta_forced_m ** 3)) + (1.0 / (delta_nat_m ** 3)) ) ** (-1.0 / 3.0)
        delta_eff_mm = delta_eff_m * 1000.0
        g_bl = (sim.D_CO2 / delta_eff_m) / sim.V_mol # mol m-2 s-1
        
        is_stagnant = (u_local_cm_s < 0.5) # air speed below 5 mm/s
        
        rec = dict(p)
        rec.update({
            'u_local_cm_s': round(u_local_cm_s, 4),
            'dC_CO2_ppm': round(dC_co2_ppm, 3),
            'dH2O_mmol': round(dH_h2o_mmol, 3),
            'delta_bl_mm': round(delta_eff_mm, 3),
            'g_bl_mol_m2_s': round(g_bl, 4),
            'is_stagnant': is_stagnant
        })
        plant_data.append(rec)
        
    return plant_data

def compute_distribution_statistics(df_plants, scenario_name, gravity_g, mode):
    """
    Compute statistical distribution metrics matching LunarLeaf-CFD T18 format.
    """
    stats = []
    
    # We analyze overall canopy, plus chamber-by-chamber breakdown
    groups = [
        ('all_canopy', df_plants),
        ('ch1_brapa_windward', df_plants[df_plants['chamber'] == 1]),
        ('ch2_arabidopsis_lawn', df_plants[df_plants['chamber'] == 2]),
        ('ch3_brapa_leeward', df_plants[df_plants['chamber'] == 3])
    ]
    
    for grp_name, grp_df in groups:
        n_p = len(grp_df)
        gbl = grp_df['g_bl_mol_m2_s'].values
        u = grp_df['u_local_cm_s'].values
        co2 = grp_df['dC_CO2_ppm'].values
        h2o = grp_df['dH2O_mmol'].values
        
        # Percentiles
        p10_gbl, p50_gbl, p90_gbl = np.percentile(gbl, [10, 50, 90])
        p10_co2, p50_co2, p90_co2 = np.percentile(co2, [10, 50, 90])
        p10_u, p50_u, p90_u = np.percentile(u, [10, 50, 90])
        
        mean_gbl = np.mean(gbl)
        std_gbl = np.std(gbl)
        cov_gbl = std_gbl / (mean_gbl + 1e-9)
        
        mean_co2 = np.mean(co2)
        std_co2 = np.std(co2)
        cov_co2 = std_co2 / (mean_co2 + 1e-9)
        
        stagnant_fraction = np.mean(grp_df['is_stagnant'].values)
        
        # Dynamic range (p90 / p10)
        dyn_range_co2 = p90_co2 / (p10_co2 + 1e-6)
        dyn_range_gbl = p90_gbl / (p10_gbl + 1e-6)
        
        stats.append({
            'scenario': scenario_name,
            'subset': grp_name,
            'plant_count': n_p,
            'gravity_g': gravity_g,
            'flow_mode': mode,
            'u_mean_cm_s': round(float(np.mean(u)), 3),
            'u_p10_cm_s': round(float(p10_u), 3),
            'u_p50_cm_s': round(float(p50_u), 3),
            'u_p90_cm_s': round(float(p90_u), 3),
            'g_bl_mean': round(float(mean_gbl), 4),
            'g_bl_p10': round(float(p10_gbl), 4),
            'g_bl_p50': round(float(p50_gbl), 4),
            'g_bl_p90': round(float(p90_gbl), 4),
            'g_bl_cov': round(float(cov_gbl), 4),
            'g_bl_dyn_range_p90_p10': round(float(dyn_range_gbl), 3),
            'dC_CO2_mean_ppm': round(float(mean_co2), 3),
            'dC_CO2_peak_ppm': round(float(np.max(co2)), 3),
            'dC_CO2_p10_ppm': round(float(p10_co2), 3),
            'dC_CO2_p50_ppm': round(float(p50_co2), 3),
            'dC_CO2_p90_ppm': round(float(p90_co2), 3),
            'dC_CO2_cov': round(float(cov_co2), 4),
            'dC_CO2_dyn_range_p90_p10': round(float(dyn_range_co2), 3),
            'stagnant_plant_fraction': round(float(stagnant_fraction), 3)
        })
        
    return stats

def main():
    print("=" * 70)
    print("LEAF GBC1 Pot Model - Foliar Surface Hotspot & Microenvironmental Analytics")
    print("=" * 70)
    
    plants = generate_164_plant_coordinates()
    print(f"Loaded {len(plants)} individual plant coordinates (32 B. rapa + 100 Arabidopsis + 32 B. rapa)")
    
    sim = GBC1CFDSimulator('geometry/gbc1_grid.npz')
    
    scenarios = [
        ('Lunar_Pull_0.5Lmin', 0.50, 1.62, 'pull'),
        ('Lunar_Push_0.5Lmin', 0.50, 1.62, 'push'),
        ('Earth_Pull_0.5Lmin', 0.50, 9.81, 'pull'),
        ('Microgravity_Pull_0.5Lmin', 0.50, 0.00, 'pull')
    ]
    
    all_telemetry = []
    all_distribution_stats = []
    nominal_df = None
    
    for sc_name, Q_flow, g_accel, mode in scenarios:
        print(f"\nSimulating scenario: {sc_name} (Q={Q_flow} L/min, g={g_accel} m/s2, {mode} mode)...")
        res, fields = sim.run_case(Q_L_min=Q_flow, gravity_g=g_accel, max_steps=350, mode=mode)
        
        plant_metrics = sample_plant_microenvironment(sim, fields, plants, g_accel)
        df_sc = pd.DataFrame(plant_metrics)
        df_sc['scenario'] = sc_name
        all_telemetry.append(df_sc)
        
        if sc_name == 'Lunar_Pull_0.5Lmin':
            nominal_df = df_sc.copy()
            
        stats_sc = compute_distribution_statistics(df_sc, sc_name, g_accel, mode)
        all_distribution_stats.extend(stats_sc)
        
    df_all_telemetry = pd.concat(all_telemetry, ignore_index=True)
    df_dist = pd.DataFrame(all_distribution_stats)
    
    os.makedirs('results/tables', exist_ok=True)
    telemetry_path = 'results/tables/gbc1_164plants_foliar_telemetry.csv'
    dist_path = 'results/tables/gbc1_foliar_surface_distribution.csv'
    
    df_all_telemetry.to_csv(telemetry_path, index=False)
    df_dist.to_csv(dist_path, index=False)
    df_dist.to_csv('results/tables/gbc2_foliar_surface_distribution.csv', index=False)
    print(f"\n[OK] Saved individual plant telemetry to {telemetry_path}")
    print(f"[OK] Saved foliar surface distribution table to {dist_path}")
    
    # Print summary of nominal Lunar pull scenario
    print("\n=== NOMINAL LUNAR PULL (Q=0.5 L/min) SURFACE DISTRIBUTION SUMMARY ===")
    lunar_summary = df_dist[df_dist['scenario'] == 'Lunar_Pull_0.5Lmin']
    cols = ['subset', 'plant_count', 'u_mean_cm_s', 'g_bl_mean', 'g_bl_cov', 'g_bl_dyn_range_p90_p10', 'dC_CO2_mean_ppm', 'stagnant_plant_fraction']
    print(lunar_summary[cols].to_string(index=False))

    # Generate Figure 17: Foliar Surface Distribution and Boundary Layer Hotspots
    print("\nGenerating Figure 17 (figures/fig17_foliar_surface_distribution_and_hotspots.png & .svg)...")
    
    fig = plt.figure(figsize=(15, 11), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.1, 1.0], hspace=0.28, wspace=0.22)
    
    # Color palette
    c_ch1 = '#0284c7' # Ocean blue (B. rapa windward)
    c_ch2 = '#10b981' # Emerald green (Arabidopsis central)
    c_ch3 = '#8b5cf6' # Purple (B. rapa leeward)
    
    # Panel A: Spatial 2D Map of 164 Plants Colored by Boundary Layer Conductance
    ax_a = fig.add_subplot(gs[0, 0])
    sc = ax_a.scatter(nominal_df['z_mm'], nominal_df['x_mm'], 
                      c=nominal_df['g_bl_mol_m2_s'], cmap='viridis', 
                      s=nominal_df['species'].map({'Brassica rapa': 75, 'Arabidopsis thaliana': 35}),
                      edgecolor='black', linewidth=0.6, alpha=0.9, zorder=3)
    
    cbar = plt.colorbar(sc, ax=ax_a, pad=0.02)
    cbar.set_label('Conductance g_bl (mol m^-2 s^-1)', fontsize=10, weight='semibold')
    cbar.ax.tick_params(labelsize=9)
    
    # Mark chamber partitions along Z
    ax_a.axvline(-8.5, color='#94a3b8', linestyle='--', lw=1.2, zorder=2)
    ax_a.axvline(+8.5, color='#94a3b8', linestyle='--', lw=1.2, zorder=2)
    ax_a.text(-19.0, 19.5, 'Chamber 1\n32 B. rapa', color=c_ch1, ha='center', va='top', fontsize=9, weight='bold')
    ax_a.text(0.0, 19.5, 'Chamber 2\n100 Arabidopsis', color=c_ch2, ha='center', va='top', fontsize=9, weight='bold')
    ax_a.text(+19.0, 19.5, 'Chamber 3\n32 B. rapa', color=c_ch3, ha='center', va='top', fontsize=9, weight='bold')
    
    # Mark inlet / outlet ports
    ax_a.scatter([-21.59], [-15.49], marker='s', s=120, color='#38bdf8', edgecolor='black', zorder=4, label='Port 1 (Inlet)')
    ax_a.scatter([0.0], [-15.49], marker='v', s=140, color='#ef4444', edgecolor='black', zorder=4, label='Port 2 (Central Pull)')
    ax_a.scatter([+21.59], [-15.49], marker='o', s=100, color='#64748b', edgecolor='black', zorder=4, label='Port 3 (RNAlater)')
    
    ax_a.set_xlabel('Chamber Axis Z (mm) [Flow Direction: Port 1 -> Port 2]', fontsize=10, weight='semibold')
    ax_a.set_ylabel('Transverse Axis X (mm)', fontsize=10, weight='semibold')
    ax_a.set_title('A  Spatial Microenvironment of 164 Plants (Lunar Pull 0.5 L/min)', fontsize=11, weight='bold', loc='left')
    ax_a.legend(loc='lower left', fontsize=8, framealpha=0.9)
    ax_a.set_xlim(-33, 33)
    ax_a.set_ylim(-18, 22)
    ax_a.grid(True, linestyle=':', alpha=0.5)

    # Panel B: Cumulative Distribution Functions (CDF) of Boundary Layer Conductance
    ax_b = fig.add_subplot(gs[0, 1])
    for ch_idx, (ch_num, c_col, sp_name) in enumerate([(1, c_ch1, 'Ch 1: B. rapa (n=32)'), 
                                                      (2, c_ch2, 'Ch 2: Arabidopsis (n=100)'), 
                                                      (3, c_ch3, 'Ch 3: B. rapa (n=32)')]):
        sub_vals = np.sort(nominal_df[nominal_df['chamber'] == ch_num]['g_bl_mol_m2_s'].values)
        cdf = np.linspace(0.01, 1.0, len(sub_vals))
        ax_b.plot(sub_vals, cdf, color=c_col, lw=2.2, label=sp_name)
        
        # Mark p10, p50, p90
        p10, p50, p90 = np.percentile(sub_vals, [10, 50, 90])
        ax_b.scatter([p50], [0.5], color=c_col, s=40, zorder=4)
        ax_b.axvline(p50, color=c_col, linestyle=':', alpha=0.4)
        
    ax_b.set_xlabel('Boundary Layer Conductance g_bl (mol m^-2 s^-1)', fontsize=10, weight='semibold')
    ax_b.set_ylabel('Cumulative Plant Fraction', fontsize=10, weight='semibold')
    ax_b.set_title('B  Foliar Conductance CDF & Percentiles (p10, p50, p90)', fontsize=11, weight='bold', loc='left')
    ax_b.axhline(0.5, color='#94a3b8', linestyle='--', lw=0.9, label='Median p50')
    ax_b.axhline(0.1, color='#cbd5e1', linestyle=':', lw=0.8)
    ax_b.axhline(0.9, color='#cbd5e1', linestyle=':', lw=0.8)
    ax_b.legend(loc='lower right', fontsize=8.5, framealpha=0.9)
    ax_b.grid(True, linestyle=':', alpha=0.5)

    # Panel C: Foliar Hotspot Dynamic Range (p90/p10) & Heterogeneity CoV Across Regimes
    ax_c = fig.add_subplot(gs[1, 0])
    scenarios_plot = [
        ('Lunar_Pull_0.5Lmin', 'Lunar Pull\n(Nominal Flight)', '#10b981'),
        ('Lunar_Push_0.5Lmin', 'Lunar Push\n(Reverse Flow)', '#f59e0b'),
        ('Earth_Pull_0.5Lmin', 'Earth 1g Pull\n(Buoyant Baseline)', '#3b82f6'),
        ('Microgravity_Pull_0.5Lmin', '0g Microgravity\n(Pure Forced)', '#ec4899')
    ]
    
    x_pos = np.arange(len(scenarios_plot))
    width = 0.35
    
    dyn_ranges = []
    covs = []
    labels = []
    colors = []
    
    for sc_id, sc_lbl, sc_c in scenarios_plot:
        row = df_dist[(df_dist['scenario'] == sc_id) & (df_dist['subset'] == 'all_canopy')].iloc[0]
        dyn_ranges.append(row['dC_CO2_dyn_range_p90_p10'])
        covs.append(row['dC_CO2_cov'] * 100.0) # in percent
        labels.append(sc_lbl)
        colors.append(sc_c)
        
    bars1 = ax_c.bar(x_pos - width/2, dyn_ranges, width, label='Dynamic Range (p90 / p10)', color='#0284c7', alpha=0.85, edgecolor='black')
    bars2 = ax_c.bar(x_pos + width/2, covs, width, label='Spatial CoV (%)', color='#f97316', alpha=0.85, edgecolor='black')
    
    for bar in bars1:
        yval = bar.get_height()
        ax_c.text(bar.get_x() + bar.get_width()/2.0, yval + 0.05, f'{yval:.2f}x', ha='center', va='bottom', fontsize=8.5, weight='bold')
    for bar in bars2:
        yval = bar.get_height()
        ax_c.text(bar.get_x() + bar.get_width()/2.0, yval + 0.5, f'{yval:.1f}%', ha='center', va='bottom', fontsize=8.5, weight='bold')
        
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels(labels, fontsize=9)
    ax_c.set_ylabel('Heterogeneity Index', fontsize=10, weight='semibold')
    ax_c.set_title('C  Canopy CO2 Draw-down Dynamic Range (p90/p10) vs Spatial CoV', fontsize=11, weight='bold', loc='left')
    ax_c.set_ylim(0, max(dyn_ranges + covs) * 1.25)
    ax_c.legend(loc='upper right', fontsize=8.5, framealpha=0.9)
    ax_c.grid(True, linestyle=':', alpha=0.5, axis='y')

    # Panel D: Foliar Conductance vs Local Air Speed Showing Stagnation Threshold
    ax_d = fig.add_subplot(gs[1, 1])
    for ch_num, c_col, sp_name in [(1, c_ch1, 'Ch 1: B. rapa (Windward)'), 
                                   (2, c_ch2, 'Ch 2: Arabidopsis (Central Lawn)'), 
                                   (3, c_ch3, 'Ch 3: B. rapa (Leeward)')]:
        sub_df = nominal_df[nominal_df['chamber'] == ch_num]
        ax_d.scatter(sub_df['u_local_cm_s'], sub_df['g_bl_mol_m2_s'], 
                     color=c_col, s=45, alpha=0.85, edgecolor='black', linewidth=0.5, label=sp_name)
        
    # Boundary layer breakdown threshold at u = 0.5 cm/s
    ax_d.axvline(0.5, color='#ef4444', linestyle='--', lw=1.5, label='Stagnation Limit (0.5 cm/s)')
    ax_d.axhspan(0.0, 0.035, color='#fee2e2', alpha=0.4, label='Hypoxic / Starvation Risk')
    
    ax_d.set_xlabel('Local Foliar Air Speed (cm/s)', fontsize=10, weight='semibold')
    ax_d.set_ylabel('Boundary Layer Conductance g_bl (mol m^-2 s^-1)', fontsize=10, weight='semibold')
    ax_d.set_title('D  Local Boundary Layer Resistance vs Air Speed (Churchill Mixed Model)', fontsize=11, weight='bold', loc='left')
    ax_d.legend(loc='lower right', fontsize=8.5, framealpha=0.9)
    ax_d.grid(True, linestyle=':', alpha=0.5)

    os.makedirs('figures', exist_ok=True)
    fig_png = 'figures/fig17_foliar_surface_distribution_and_hotspots.png'
    fig_svg = 'figures/fig17_foliar_surface_distribution_and_hotspots.svg'
    plt.tight_layout()
    plt.savefig(fig_png, dpi=300)
    plt.savefig(fig_svg, format='svg')
    plt.close()
    print(f"[OK] Successfully saved Figure 17 to {fig_png} and {fig_svg}")

if __name__ == '__main__':
    main()
