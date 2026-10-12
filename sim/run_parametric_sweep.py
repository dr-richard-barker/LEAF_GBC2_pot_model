"""
Parametric Sweep Runner for LEAF GBC2 Pot CFD Model with Enhanced Publication Figures.
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sim.gbc2_cfd_engine import GBC2CFDSimulator

def run_sweep():
    os.makedirs('results/tables', exist_ok=True)
    os.makedirs('results/figures', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    sim = GBC2CFDSimulator()
    
    gravities = [
        ('1g_Earth', 9.81),
        ('Lunar_0.166g', 1.62),
        ('Microgravity_0g', 0.0)
    ]
    
    pump_speeds = [0.05, 0.10, 0.25, 0.50, 1.00, 2.00, 4.00] # L/min
    
    records = []
    field_snapshots = {}
    
    print("=================================================================")
    print(f"Starting LEAF GBC2 Parametric CFD Sweep: {len(gravities)} gravities x {len(pump_speeds)} speeds")
    print("=================================================================")
    
    for g_label, g_val in gravities:
        for Q in pump_speeds:
            res, fields = sim.run_case(Q_L_min=Q, gravity_g=g_val, max_steps=400)
            
            row = {
                'scenario': f"{g_label}_Q{Q:.2f}",
                'gravity_condition': g_label,
                'gravity_m_s2': g_val,
                'pump_flow_L_min': Q,
                'inlet_velocity_m_s': res['U_inlet_m_s'],
                'mean_headspace_vel_cm_s': res['U_mean_headspace_cm_s'],
                'max_velocity_m_s': res['U_max_m_s'],
                'porous_vel_cm_s': res['U_porous_cm_s'],
                'stagnant_volume_pct': res['stagnant_volume_pct'],
                'delta_P_Pa': res['delta_P_Pa'],
                'flow_uniformity_gamma': res['flow_uniformity_gamma'],
                'dC_CO2_mean_ppm': res['dC_CO2_mean_ppm'],
                'dC_H2O_mean_mmol_m3': res['dC_H2O_mean_mmol'],
                'delta_bl_mm': res['delta_bl_mean_mm'],
                'g_bl_mol_m2_s': res['g_bl_mean_mol_m2_s'],
                'comp1_u_cm_s': res['compartments']['comp_1']['u_mean_cm_s'],
                'comp2_u_cm_s': res['compartments']['comp_2']['u_mean_cm_s'],
                'comp3_u_cm_s': res['compartments']['comp_3']['u_mean_cm_s'],
                'comp1_g_bl': res['compartments']['comp_1']['g_bl_mol_m2_s'],
                'comp2_g_bl': res['compartments']['comp_2']['g_bl_mol_m2_s'],
                'comp3_g_bl': res['compartments']['comp_3']['g_bl_mol_m2_s'],
            }
            records.append(row)
            
            if Q in [0.10, 0.50, 2.00]:
                field_snapshots[f"{g_label}_Q{Q:.2f}"] = {
                    'u': fields['u'], 'v': fields['v'], 'w': fields['w'],
                    'u_mag': fields['u_mag'],
                    'c_CO2': fields['c_CO2'],
                    'c_H2O': fields['c_H2O']
                }

    df = pd.DataFrame(records)
    table_path = 'results/tables/gbc2_cfd_sweep_results.tsv'
    df.to_csv(table_path, sep='\t', index=False)
    print(f"\nSaved full sweep results table to: {table_path}")
    
    # Save OSDR / LunarLeaf-CFD covariate format table
    covar_df = df[['scenario', 'gravity_condition', 'gravity_m_s2', 'pump_flow_L_min',
                   'inlet_velocity_m_s', 'g_bl_mol_m2_s', 'delta_bl_mm', 'dC_CO2_mean_ppm',
                   'stagnant_volume_pct', 'flow_uniformity_gamma']].copy()
    covar_path = 'results/tables/gbc2_covariates_for_omics.tsv'
    covar_df.to_csv(covar_path, sep='\t', index=False)
    print(f"Saved omics covariates table to: {covar_path}")

    # -------------------------------------------------------------
    # PUBLICATION FIGURES
    # -------------------------------------------------------------
    print("\nGenerating publication figures...")
    
    # Figure 1: Stagnation & Boundary Layer vs Pump Speed
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    colors = {'1g_Earth': '#1f77b4', 'Lunar_0.166g': '#ff7f0e', 'Microgravity_0g': '#2ca02c'}
    markers = {'1g_Earth': 'o', 'Lunar_0.166g': 's', 'Microgravity_0g': '^'}
    labels = {'1g_Earth': '1g Earth (9.81 m/s²)', 'Lunar_0.166g': 'Lunar (1.62 m/s²)', 'Microgravity_0g': 'Microgravity (0 m/s²)'}
    
    # Panel A: Stagnant Volume % vs Pump Flow
    ax = axes[0, 0]
    for g_label, grp in df.groupby('gravity_condition'):
        ax.plot(grp['pump_flow_L_min'], grp['stagnant_volume_pct'],
                color=colors[g_label], marker=markers[g_label], lw=2.5, ms=7, label=labels[g_label])
    ax.axhline(10.0, color='red', ls='--', alpha=0.7, label='Target Threshold (<10% Stagnant)')
    ax.set_xscale('log')
    ax.set_xlabel('Pump Flow Rate Q (L/min)', fontsize=12)
    ax.set_ylabel('Stagnant Volume Fraction (%)', fontsize=12)
    ax.set_title('A. Stagnation Clearance vs. Pump Speed', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=True, fontsize=10)
    
    # Panel B: Boundary Layer Conductance g_bl vs Pump Flow
    ax = axes[0, 1]
    for g_label, grp in df.groupby('gravity_condition'):
        ax.plot(grp['pump_flow_L_min'], grp['g_bl_mol_m2_s'],
                color=colors[g_label], marker=markers[g_label], lw=2.5, ms=7, label=labels[g_label])
    ax.set_xscale('log')
    ax.set_xlabel('Pump Flow Rate Q (L/min)', fontsize=12)
    ax.set_ylabel('Boundary Layer Conductance g_bl (mol m⁻² s⁻¹)', fontsize=12)
    ax.set_title('B. Boundary Layer Conductance (g_bl)', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=True, fontsize=10)
    
    # Panel C: CO2 Drawdown at Leaf Surface vs Pump Flow
    ax = axes[1, 0]
    for g_label, grp in df.groupby('gravity_condition'):
        ax.plot(grp['pump_flow_L_min'], grp['dC_CO2_mean_ppm'],
                color=colors[g_label], marker=markers[g_label], lw=2.5, ms=7, label=labels[g_label])
    ax.set_xscale('log')
    ax.set_xlabel('Pump Flow Rate Q (L/min)', fontsize=12)
    ax.set_ylabel('Mean Leaf CO₂ Drawdown ΔC_CO₂ (ppm)', fontsize=12)
    ax.set_title('C. Foliar CO₂ Depletion Gradient', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=True, fontsize=10)
    
    # Panel D: Boundary Layer Thickness delta vs Pump Flow
    ax = axes[1, 1]
    for g_label, grp in df.groupby('gravity_condition'):
        ax.plot(grp['pump_flow_L_min'], grp['delta_bl_mm'],
                color=colors[g_label], marker=markers[g_label], lw=2.5, ms=7, label=labels[g_label])
    ax.set_xscale('log')
    ax.set_xlabel('Pump Flow Rate Q (L/min)', fontsize=12)
    ax.set_ylabel('Boundary Layer Thickness δ (mm)', fontsize=12)
    ax.set_title('D. Diffusive Boundary Layer Thickness (δ)', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=True, fontsize=10)
    
    plt.tight_layout()
    plt.savefig('figures/fig1_pump_speed_vs_aerodynamics.png', dpi=200)
    plt.savefig('figures/fig1_pump_speed_vs_aerodynamics.svg')
    print("Saved figures/fig1_pump_speed_vs_aerodynamics.png and .svg")
    
    # Figure 2: Longitudinal Vertical Airflow Velocity Vectors & Recirculation
    # Slice at X = -5 mm (index ~ 20)
    fig, axes = plt.subplots(3, 2, figsize=(15, 14))
    i_slice = 20
    
    cases_to_plot = [
        ('Q = 0.10 L/min (Low Flow)', '1g_Earth_Q0.10', 'Lunar_0.166g_Q0.10'),
        ('Q = 0.50 L/min (Nominal Flow)', '1g_Earth_Q0.50', 'Lunar_0.166g_Q0.50'),
        ('Q = 2.00 L/min (High Flow)', '1g_Earth_Q2.00', 'Lunar_0.166g_Q2.00')
    ]
    
    for row_idx, (title_speed, k_earth, k_lunar) in enumerate(cases_to_plot):
        snap_e = field_snapshots[k_earth]
        snap_l = field_snapshots[k_lunar]
        
        # Earth slice
        u_mag_e = snap_e['u_mag'][i_slice, :, :] * 100 # cm/s
        v_e = snap_e['v'][i_slice, :, :] * 100
        w_e = snap_e['w'][i_slice, :, :] * 100
        
        # Lunar slice
        u_mag_l = snap_l['u_mag'][i_slice, :, :] * 100 # cm/s
        v_l = snap_l['v'][i_slice, :, :] * 100
        w_l = snap_l['w'][i_slice, :, :] * 100
        
        # Earth plot
        ax_e = axes[row_idx, 0]
        vmax_val = max(np.percentile(u_mag_e, 95), 1.0)
        im_e = ax_e.imshow(u_mag_e, origin='lower', cmap='viridis', extent=[-34, 34, 0, 38], vmin=0, vmax=vmax_val)
        # Quiver arrows
        skip = 3
        Z_grid, Y_grid = np.meshgrid(np.linspace(-34, 34, u_mag_e.shape[1]), np.linspace(0, 38, u_mag_e.shape[0]))
        ax_e.quiver(Z_grid[::skip, ::skip], Y_grid[::skip, ::skip],
                   w_e[::skip, ::skip], v_e[::skip, ::skip], color='white', alpha=0.7, scale=vmax_val*15)
        ax_e.set_title(f'1g Earth: {title_speed}\nLongitudinal Airflow Vector Field (X = -5mm)', fontsize=11, fontweight='bold')
        ax_e.set_xlabel('Z Length (mm)')
        ax_e.set_ylabel('Y Height (mm)')
        plt.colorbar(im_e, ax=ax_e, fraction=0.046, pad=0.04, label='Air Speed (cm/s)')
        
        # Lunar plot
        ax_l = axes[row_idx, 1]
        vmax_val_l = max(np.percentile(u_mag_l, 95), 1.0)
        im_l = ax_l.imshow(u_mag_l, origin='lower', cmap='viridis', extent=[-34, 34, 0, 38], vmin=0, vmax=vmax_val_l)
        ax_l.quiver(Z_grid[::skip, ::skip], Y_grid[::skip, ::skip],
                   w_l[::skip, ::skip], v_l[::skip, ::skip], color='white', alpha=0.7, scale=vmax_val_l*15)
        ax_l.set_title(f'Lunar (0.166g): {title_speed}\nLongitudinal Airflow Vector Field (X = -5mm)', fontsize=11, fontweight='bold')
        ax_l.set_xlabel('Z Length (mm)')
        ax_l.set_ylabel('Y Height (mm)')
        plt.colorbar(im_l, ax=ax_l, fraction=0.046, pad=0.04, label='Air Speed (cm/s)')

        for s_ax in [ax_e, ax_l]:
            s_ax.add_patch(patches.Rectangle((-30.5, 7.62), 21.0, 7.0, linewidth=1.2, edgecolor='#e67e22', facecolor='none', ls='--'))
            s_ax.add_patch(patches.Rectangle((-7.5, 7.62), 15.0, 7.0, linewidth=1.2, edgecolor='#e67e22', facecolor='none', ls='--'))
            s_ax.add_patch(patches.Rectangle((9.5, 7.62), 21.0, 7.0, linewidth=1.2, edgecolor='#e67e22', facecolor='none', ls='--', label='Substrates (Y=7.6-14.6mm)'))
            s_ax.plot([-9.5, -9.5], [7.62, 24.0], color='white', lw=1.2, alpha=0.5)
            s_ax.plot([9.5, 9.5], [7.62, 24.0], color='white', lw=1.2, alpha=0.5)
        
    plt.tight_layout()
    plt.savefig('figures/fig2_velocity_fields_earth_vs_lunar.png', dpi=200)
    plt.savefig('figures/fig2_velocity_fields_earth_vs_lunar.svg')
    print("Saved figures/fig2_velocity_fields_earth_vs_lunar.png and .svg")
    
    # Figure 3: Lunar Compensatory Pump Speed Curve
    fig, ax = plt.subplots(figsize=(9, 6))
    
    earth_data = df[df['gravity_condition'] == '1g_Earth'].sort_values('pump_flow_L_min')
    lunar_data = df[df['gravity_condition'] == 'Lunar_0.166g'].sort_values('pump_flow_L_min')
    
    earth_Q_vals = earth_data['pump_flow_L_min'].values
    earth_gbl_vals = earth_data['g_bl_mol_m2_s'].values
    lunar_Q_vals = lunar_data['pump_flow_L_min'].values
    lunar_gbl_vals = lunar_data['g_bl_mol_m2_s'].values
    
    # To match Earth g_bl at low pump speed, what Q on Moon is required?
    # Log-log interpolation
    lunar_Q_needed = np.exp(np.interp(np.log(earth_gbl_vals), np.log(lunar_gbl_vals), np.log(lunar_Q_vals)))
    
    ax.plot(earth_Q_vals, earth_Q_vals, 'k--', lw=2, label='1:1 Line (Zero Gravity Penalty)')
    ax.plot(earth_Q_vals, lunar_Q_needed, 'r-o', lw=2.5, ms=8, label='Required Lunar Pump Flow Rate')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Baseline Earth 1g Pump Flow Rate Q_Earth (L/min)', fontsize=12)
    ax.set_ylabel('Required Lunar Pump Flow Rate Q_Lunar (L/min)', fontsize=12)
    ax.set_title('Lunar Convection Penalty: Compensatory Pump Sizing Curve', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=True, fontsize=11)
    
    penalty_ratio = lunar_Q_needed[0] / earth_Q_vals[0]
    ax.text(0.08, 0.78, f'Lunar Convection Deficit:\nAt low ventilation (Q = 0.05 - 0.1 L/min),\nLunar hardware requires {penalty_ratio:.1f}x higher pump flow\nto equal 1g Earth boundary-layer conductance.',
            transform=ax.transAxes, bbox=dict(boxstyle='round,pad=0.6', facecolor='#ffffcc', edgecolor='#cccc99'), fontsize=10)
            
    plt.tight_layout()
    plt.savefig('figures/fig3_lunar_compensatory_pump_curve.png', dpi=200)
    plt.savefig('figures/fig3_lunar_compensatory_pump_curve.svg')
    print("Saved figures/fig3_lunar_compensatory_pump_curve.png and .svg")
    
    print("\nAll simulations and figures completed successfully!")
    return df

if __name__ == '__main__':
    run_sweep()
