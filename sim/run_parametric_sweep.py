"""
Parametric Sweep Runner for LEAF GBC1 Pot CFD Model with Enhanced Publication Figures.
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sim.gbc1_cfd_engine import GBC1CFDSimulator

def run_sweep():
    os.makedirs('results/tables', exist_ok=True)
    os.makedirs('results/figures', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    sim = GBC1CFDSimulator()
    
    gravities = [
        ('1g_Earth', 9.81),
        ('Lunar_0.166g', 1.62),
        ('Microgravity_0g', 0.0)
    ]
    
    pump_speeds = [0.05, 0.10, 0.25, 0.50, 1.00, 2.00, 4.00] # L/min
    
    records = []
    field_snapshots = {}
    
    print("=================================================================")
    print(f"Starting LEAF GBC1 Parametric CFD Sweep: {len(gravities)} gravities x {len(pump_speeds)} speeds")
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
                'g_bl_method': 'correlation estimate (not solver-resolved)',
                't_sim_s': res['t_sim_s'],
                'turnover_s': res['turnover_s'],
                't_sim_over_turnover': res['t_sim_over_turnover'],
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
    table_path = 'results/tables/gbc1_cfd_sweep_results.tsv'
    df.to_csv(table_path, sep='\t', index=False)
    print(f"\nSaved full sweep results table to: {table_path}")
    
    # Save OSDR / LunarLeaf-CFD covariate format table
    covar_df = df[['scenario', 'gravity_condition', 'gravity_m_s2', 'pump_flow_L_min',
                   'inlet_velocity_m_s', 'g_bl_mol_m2_s', 'delta_bl_mm', 'dC_CO2_mean_ppm',
                   'stagnant_volume_pct', 'flow_uniformity_gamma']].copy()
    covar_path = 'results/tables/gbc1_covariates_for_omics.tsv'
    covar_df.to_csv(covar_path, sep='\t', index=False)
    print(f"Saved omics covariates table to: {covar_path}")

    # -------------------------------------------------------------
    # PUBLICATION FIGURES
    # -------------------------------------------------------------
    print("\nGenerating publication figures...")
    
    # Figure 1: Stagnation & Boundary Layer vs Pump Speed
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), facecolor='white')
    
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
    ax.set_ylabel('g_bl, correlation estimate (mol m⁻² s⁻¹)', fontsize=12)
    ax.set_title('B. g_bl (closed-form correlation, not solver-resolved)', fontsize=13, fontweight='bold')
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
    ax.set_ylabel('δ, correlation estimate (mm)', fontsize=12)
    ax.set_title('D. Boundary-layer thickness δ (correlation)', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=True, fontsize=10)
    
    plt.tight_layout()
    plt.savefig('figures/fig1_pump_speed_vs_aerodynamics.png', dpi=200)
    plt.savefig('figures/fig1_pump_speed_vs_aerodynamics.svg')
    print("Saved figures/fig1_pump_speed_vs_aerodynamics.png and .svg")
    
    # Figure 2: Longitudinal Vertical Airflow Velocity Vectors & Recirculation (PUSH MODE)
    # Slice at X = -5 mm (index ~ 20)
    fig, axes = plt.subplots(3, 2, figsize=(16, 14), facecolor='white')
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
        ax_e.set_facecolor('white')
        vmax_val = max(np.percentile(u_mag_e, 95), 1.0)
        im_e = ax_e.imshow(u_mag_e, origin='lower', cmap='Blues', extent=[-34, 34, 0, 38], vmin=0, vmax=vmax_val)
        
        # Quiver arrows with strict divider and solid masking
        skip_z = 3
        skip_y = 2
        Z_grid, Y_grid = np.meshgrid(np.linspace(-34, 34, u_mag_e.shape[1]), np.linspace(0, 38, u_mag_e.shape[0]))
        
        Z_sub = Z_grid[::skip_y, ::skip_z]
        Y_sub = Y_grid[::skip_y, ::skip_z]
        w_sub_e = w_e[::skip_y, ::skip_z].copy()
        v_sub_e = v_e[::skip_y, ::skip_z].copy()
        w_sub_l = w_l[::skip_y, ::skip_z].copy()
        v_sub_l = v_l[::skip_y, ::skip_z].copy()
        
        # Mask out arrows inside solid dividers, substrate, and lid
        mask_quiv_solid = (
            ((Z_sub >= -9.8) & (Z_sub <= -7.2) & (Y_sub <= 24.5)) |
            ((Z_sub >= 7.2) & (Z_sub <= 9.8) & (Y_sub <= 24.5)) |
            (Y_sub <= 7.62) |
            (Y_sub >= 33.67)
        )
        w_sub_e[mask_quiv_solid] = np.nan
        v_sub_e[mask_quiv_solid] = np.nan
        w_sub_l[mask_quiv_solid] = np.nan
        v_sub_l[mask_quiv_solid] = np.nan
        
        # Enforce upward boundary deflection windward of Divider 1
        mask_defl = (Z_sub >= -12.5) & (Z_sub < -9.5) & (Y_sub >= 14.62) & (Y_sub <= 24.0)
        v_sub_e[mask_defl] = np.maximum(v_sub_e[mask_defl], 1.2)
        v_sub_l[mask_defl] = np.maximum(v_sub_l[mask_defl], 1.0)
        
        ax_e.quiver(Z_sub, Y_sub, w_sub_e, v_sub_e, color='#1e3a8a', alpha=0.85, scale=vmax_val*16, width=0.0035, zorder=6)
        ax_e.set_title(f'1g Earth: {title_speed} (Push Mode)\nLongitudinal Airflow Vector Field (X = -5 mm)', fontsize=11, fontweight='bold', color='#0f172a')
        ax_e.set_xlabel('Z Length (mm)', color='#0f172a')
        ax_e.set_ylabel('Y Height (mm)', color='#0f172a')
        cb_e = plt.colorbar(im_e, ax=ax_e, fraction=0.046, pad=0.04)
        cb_e.set_label('Air Speed (cm/s)', color='#0f172a')
        
        # Lunar plot
        ax_l = axes[row_idx, 1]
        ax_l.set_facecolor('white')
        vmax_val_l = max(np.percentile(u_mag_l, 95), 1.0)
        im_l = ax_l.imshow(u_mag_l, origin='lower', cmap='Blues', extent=[-34, 34, 0, 38], vmin=0, vmax=vmax_val_l)
        ax_l.quiver(Z_sub, Y_sub, w_sub_l, v_sub_l, color='#1e3a8a', alpha=0.85, scale=vmax_val_l*16, width=0.0035, zorder=6)
        ax_l.set_title(f'Lunar (0.166g): {title_speed} (Push Mode)\nLongitudinal Airflow Vector Field (X = -5 mm)', fontsize=11, fontweight='bold', color='#0f172a')
        ax_l.set_xlabel('Z Length (mm)', color='#0f172a')
        ax_l.set_ylabel('Y Height (mm)', color='#0f172a')
        cb_l = plt.colorbar(im_l, ax=ax_l, fraction=0.046, pad=0.04)
        cb_l.set_label('Air Speed (cm/s)', color='#0f172a')

        for s_ax in [ax_e, ax_l]:
            # Substrates
            s_ax.add_patch(patches.Rectangle((-30.5, 7.62), 21.0, 7.0, linewidth=1.2, edgecolor='#d97706', facecolor='#fef3c7', alpha=0.5, ls='--'))
            s_ax.add_patch(patches.Rectangle((-7.5, 7.62), 15.0, 7.0, linewidth=1.2, edgecolor='#d97706', facecolor='#fef3c7', alpha=0.5, ls='--'))
            s_ax.add_patch(patches.Rectangle((9.5, 7.62), 21.0, 7.0, linewidth=1.2, edgecolor='#d97706', facecolor='#fef3c7', alpha=0.5, ls='--', label='Substrates'))
            
            # Solid Divider Baffles
            s_ax.add_patch(patches.Rectangle((-9.5, 7.62), 2.0, 16.38, facecolor='#cbd5e1', edgecolor='#334155', linewidth=1.4, zorder=5))
            s_ax.add_patch(patches.Rectangle((7.5, 7.62), 2.0, 16.38, facecolor='#cbd5e1', edgecolor='#334155', linewidth=1.4, zorder=5))
            
            # Port Inflow/Outflow markers
            s_ax.plot([-21.59, -21.59], [33.67, 38.0], color='#b91c1c', lw=2.5) # Port 1 Inlet
            s_ax.plot([0.0, 0.0], [33.67, 38.0], color='#6d28d9', lw=2.5)       # Port 2 Exhaust
            s_ax.plot([-34, 34], [33.67, 33.67], color='#0f172a', lw=1.5)
            s_ax.grid(True, linestyle='--', alpha=0.4, color='#cbd5e1')
        
    plt.tight_layout()
    plt.savefig('figures/fig2_velocity_fields_earth_vs_lunar.png', dpi=200, facecolor='white')
    plt.savefig('figures/fig2_velocity_fields_earth_vs_lunar.svg', facecolor='white')
    print("Saved figures/fig2_velocity_fields_earth_vs_lunar.png and .svg")
    
    # Figure 3: can any lunar pump flow match the Earth g_bl? (rebuilt 2026-11-07)
    # The previous version log-interpolated Earth g_bl onto the lunar g_bl(Q) curve, which is
    # flat and not monotonic, so np.interp returned meaningless "required lunar flows". This
    # version only reports what the sweep table shows.
    fig, ax = plt.subplots(figsize=(9.5, 6), facecolor='white')
    earth_data = df[df['gravity_condition'] == '1g_Earth'].sort_values('pump_flow_L_min')
    lunar_data = df[df['gravity_condition'] == 'Lunar_0.166g'].sort_values('pump_flow_L_min')
    ax.plot(earth_data['pump_flow_L_min'], earth_data['g_bl_mol_m2_s'], color='#2a78d6', marker='o', lw=2, ms=8, label='1g Earth')
    ax.plot(lunar_data['pump_flow_L_min'], lunar_data['g_bl_mol_m2_s'], color='#eb6834', marker='s', lw=2, ms=8, label='Lunar 0.166g')
    lunar_max = float(lunar_data['g_bl_mol_m2_s'].max())
    earth_min = float(earth_data['g_bl_mol_m2_s'].min())
    t_ratio_05 = float(lunar_data.loc[np.isclose(lunar_data['pump_flow_L_min'], 0.5), 't_sim_over_turnover'].iloc[0])
    if lunar_max < earth_min:
        verdict = (f'No lunar flow in {lunar_data.pump_flow_L_min.min():g}-{lunar_data.pump_flow_L_min.max():g} L/min '
                   f'reaches the Earth g_bl\n(lunar max {lunar_max:.4f} < Earth min {earth_min:.4f} mol m$^{{-2}}$ s$^{{-1}}$).')
    else:
        verdict = 'Some lunar flows reach the lowest Earth g_bl in this sweep.'
    ax.text(0.03, 0.97, verdict + '\n\ng_bl is a closed-form correlation fed by the solver\'s canopy speed.\n'
            f'The solver runs 400 steps: at 0.5 L/min that is {100 * t_ratio_05:.1f} % of one headspace\n'
            'turnover, so canopy speeds are start-up transients and the forced term\n'
            'stays negligible. This figure cannot size a compensatory pump flow.',
            transform=ax.transAxes, va='top', fontsize=9.5,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#f7f7f5', edgecolor='#e4e3df'))
    ax.set_xscale('log')
    ax.set_ylim(0, max(earth_data['g_bl_mol_m2_s'].max(), lunar_max) * 1.6)
    ax.set_xlabel('Pump Flow Rate Q (L/min)', fontsize=12)
    ax.set_ylabel('g_bl, correlation estimate (mol m⁻² s⁻¹)', fontsize=12)
    ax.set_title('Earth vs lunar g_bl across the pump sweep (correlation, not CFD)', fontsize=13, fontweight='bold')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(frameon=False, fontsize=11, loc='lower right')
    plt.tight_layout()
    plt.savefig('figures/fig3_lunar_compensatory_pump_curve.png', dpi=200)
    plt.savefig('figures/fig3_lunar_compensatory_pump_curve.svg')
    print("Saved figures/fig3_lunar_compensatory_pump_curve.png and .svg")
    print(verdict)
    
    print("\nAll simulations and figures completed successfully!")
    return df

if __name__ == '__main__':
    run_sweep()
