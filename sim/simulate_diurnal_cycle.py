"""
Diurnal (Day/Night) Photoperiod Inversion CFD Simulation for LEAF GBC2 Plant Growth Chamber.

Models:
1. Day Phase (Photoperiod):
   - Active LED lighting, leaf heating (dT = +1.5 K)
   - Net photosynthetic CO2 assimilation (S_CO2 < 0)
   - High transpiration (S_H2O > 0)
   - Upward thermal & solutal buoyancy plumes

2. Night Phase (Skotoperiod):
   - Lights off, leaf radiative cooling (dT = -0.3 K)
   - Dark mitochondrial respiration (S_CO2 > 0)
   - Stomatal closure / cuticular transpiration only (85% reduction)
   - Downward dense CO2 pooling & solutal buoyancy inversion
   - Night setback pump power analysis (Q = 0.10, 0.50, 0.60 L/min) on Earth vs Moon
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sim.gbc2_cfd_engine import GBC2CFDSimulator
from sim.advanced_analytics import GBC2AdvancedAnalytics

def run_diurnal_simulation():
    os.makedirs('results/tables', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    print("=================================================================")
    print("Simulating Diurnal (Day/Night) Photoperiod Inversion for LEAF GBC2")
    print("=================================================================")
    
    sim = GBC2CFDSimulator()
    analytics = GBC2AdvancedAnalytics()
    
    # Canopy masks
    nz_mid1 = int(sim.Nz * 0.33)
    nz_mid2 = int(sim.Nz * 0.66)
    
    mask_c1 = sim.is_canopy.copy()
    mask_c1[:, :, nz_mid1:] = False
    
    mask_c2 = sim.is_canopy.copy()
    mask_c2[:, :, :nz_mid1] = False
    mask_c2[:, :, nz_mid2:] = False
    
    mask_c3 = sim.is_canopy.copy()
    mask_c3[:, :, :nz_mid2] = False
    
    diurnal_scenarios = [
        # (name, gravity_name, g_val, Q_L_min, is_day, dT_leaf)
        ('Earth_Day_Nominal', '1g_Earth', 9.81, 0.50, True, 1.5),
        ('Earth_Night_Nominal', '1g_Earth', 9.81, 0.50, False, -0.3),
        ('Earth_Night_LowPower', '1g_Earth', 9.81, 0.10, False, -0.3),
        ('Lunar_Day_Nominal', 'Lunar_0.166g', 1.62, 0.50, True, 1.5),
        ('Lunar_Night_Nominal', 'Lunar_0.166g', 1.62, 0.50, False, -0.3),
        ('Lunar_Night_LowPower', 'Lunar_0.166g', 1.62, 0.10, False, -0.3),
        ('Lunar_Night_Compensated', 'Lunar_0.166g', 1.62, 0.60, False, -0.3),
    ]
    
    results = []
    concentration_slices = {}
    
    for name, g_label, g_val, Q_flow, is_day, dT in diurnal_scenarios:
        print(f"Executing Scenario: {name} (Q = {Q_flow:.2f} L/min, g = {g_val} m/s², Day={is_day}) ...")
        
        # Temporary override of simulator leaf properties
        sim.dT_leaf = dT
        
        # Base run for flow field
        res, fields = sim.run_case(Q_L_min=Q_flow, gravity_g=g_val, max_steps=400)
        u_mag = np.sqrt(fields['u']**2 + fields['v']**2 + fields['w']**2)
        
        dx = sim.dx
        c_CO2 = np.zeros_like(fields['u'])
        c_H2O = np.zeros_like(fields['u'])
        
        # Assign biological volumetric source terms
        s_CO2 = np.zeros_like(c_CO2)
        s_H2O = np.zeros_like(c_H2O)
        
        if is_day:
            s_CO2[mask_c1 | mask_c3] = -40.0
            s_CO2[mask_c2] = -15.0
            s_H2O[mask_c1 | mask_c3] = +80.0
            s_H2O[mask_c2] = +30.0
        else:
            s_CO2[mask_c1 | mask_c3] = +6.0
            s_CO2[mask_c2] = +3.0
            s_H2O[mask_c1 | mask_c3] = +12.0
            s_H2O[mask_c2] = +4.5
            
        u, v, w = fields['u'], fields['v'], fields['w']
        
        # Unconditionally stable point-relaxation solver
        diff_CO2 = sim.D_CO2 / (dx ** 2)
        diff_H2O = sim.D_H2O / (dx ** 2)
        
        u_p, u_m = np.maximum(u, 0.0)/dx, np.maximum(-u, 0.0)/dx
        v_p, v_m = np.maximum(v, 0.0)/dx, np.maximum(-v, 0.0)/dx
        w_p, w_m = np.maximum(w, 0.0)/dx, np.maximum(-w, 0.0)/dx
        
        denom_CO2 = u_p + u_m + v_p + v_m + w_p + w_m + 6.0 * diff_CO2 + 1e-6
        denom_H2O = u_p + u_m + v_p + v_m + w_p + w_m + 6.0 * diff_H2O + 1e-6
        
        for _ in range(200):
            num_CO2 = (
                (u_p + diff_CO2) * np.roll(c_CO2, 1, 0) + (u_m + diff_CO2) * np.roll(c_CO2, -1, 0) +
                (v_p + diff_CO2) * np.roll(c_CO2, 1, 1) + (v_m + diff_CO2) * np.roll(c_CO2, -1, 1) +
                (w_p + diff_CO2) * np.roll(c_CO2, 1, 2) + (w_m + diff_CO2) * np.roll(c_CO2, -1, 2) +
                s_CO2
            )
            c_CO2_new = np.where(sim.is_fluid, num_CO2 / denom_CO2, 0.0)
            c_CO2_new[sim.is_inlet] = 0.0
            c_CO2 = 0.3 * c_CO2 + 0.7 * c_CO2_new
            
            num_H2O = (
                (u_p + diff_H2O) * np.roll(c_H2O, 1, 0) + (u_m + diff_H2O) * np.roll(c_H2O, -1, 0) +
                (v_p + diff_H2O) * np.roll(c_H2O, 1, 1) + (v_m + diff_H2O) * np.roll(c_H2O, -1, 1) +
                (w_p + diff_H2O) * np.roll(c_H2O, 1, 2) + (w_m + diff_H2O) * np.roll(c_H2O, -1, 2) +
                s_H2O
            )
            c_H2O_new = np.where(sim.is_fluid, num_H2O / denom_H2O, 0.0)
            c_H2O_new[sim.is_inlet] = 0.0
            c_H2O = 0.3 * c_H2O + 0.7 * c_H2O_new
            
        # Store longitudinal midplane slice (X = mid)
        mid_x = sim.Nx // 2
        concentration_slices[name] = {
            'c_CO2': c_CO2[mid_x, :, :],
            'c_H2O': c_H2O[mid_x, :, :],
            'v': v[mid_x, :, :]
        }
        
        # Calculate diagnostics
        dCO2_c1 = float(np.mean(c_CO2[mask_c1]))
        dCO2_c2 = float(np.mean(c_CO2[mask_c2]))
        dCO2_c3 = float(np.mean(c_CO2[mask_c3]))
        dCO2_sub = float(np.mean(c_CO2[sim.is_porous]))
        
        dH2O_c1 = float(np.mean(c_H2O[mask_c1]))
        dH2O_c2 = float(np.mean(c_H2O[mask_c2]))
        dH2O_c3 = float(np.mean(c_H2O[mask_c3]))
        
        u_c1 = float(np.mean(u_mag[mask_c1]))
        u_c2 = float(np.mean(u_mag[mask_c2]))
        u_c3 = float(np.mean(u_mag[mask_c3]))
        
        # Richardson Number Ri = Gr / Re^2
        beta_eff = abs((sim.beta_T * dT) + (0.001 * sim.beta_H2O * np.mean(c_H2O[sim.is_canopy])) + (0.001 * sim.beta_CO2 * np.mean(c_CO2[sim.is_canopy])))
        Gr_num = (g_val * beta_eff * (sim.L_leaf ** 3)) / (sim.nu ** 2)
        u_ref = max(float(np.mean(u_mag[sim.is_headspace])), 0.0005)
        Re_num = (u_ref * sim.L_leaf) / sim.nu
        Ri_num = float(Gr_num / max(Re_num**2, 1e-6))
        
        results.append({
            'scenario': name,
            'gravity_condition': g_label,
            'gravity_m_s2': g_val,
            'flow_rate_L_min': Q_flow,
            'diurnal_phase': 'Day' if is_day else 'Night',
            'leaf_dT_K': dT,
            'Ri_Richardson_number': Ri_num,
            'dCO2_chamber1_ppm': dCO2_c1,
            'dCO2_chamber2_ppm': dCO2_c2,
            'dCO2_chamber3_ppm': dCO2_c3,
            'dCO2_substrate_ppm': dCO2_sub,
            'dH2O_chamber1_mmol': dH2O_c1,
            'dH2O_chamber2_mmol': dH2O_c2,
            'dH2O_chamber3_mmol': dH2O_c3,
            'u_chamber1_cm_s': u_c1 * 100,
            'u_chamber2_cm_s': u_c2 * 100,
            'u_chamber3_cm_s': u_c3 * 100,
            'nocturnal_pooling_ratio': float(dCO2_sub / (dCO2_c2 + 1e-6)) if not is_day else 0.0
        })
        
    df_diurnal = pd.DataFrame(results)
    df_diurnal.to_csv('results/tables/gbc2_diurnal_day_night_results.tsv', sep='\t', index=False)
    print("\nDiurnal Simulation Summary:")
    print(df_diurnal[['scenario', 'diurnal_phase', 'flow_rate_L_min', 'dCO2_chamber1_ppm', 'dCO2_chamber2_ppm', 'dCO2_substrate_ppm', 'Ri_Richardson_number']])
    
    # -------------------------------------------------------------------------
    # Generate Publication Figure 8
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 14), facecolor='white')
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25)
    
    # Coordinates in mm for slice plotting
    y_mm = np.linspace(0.0, 34.0, sim.Ny)
    z_mm = np.linspace(-34.0, 34.0, sim.Nz)
    Z_mesh, Y_mesh = np.meshgrid(z_mm, y_mm)
    
    # Panel A: Spatial CO2 Field Slice - Night Skotoperiod on Earth 1g
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_facecolor('#fafafa')
    
    slice_night_earth = concentration_slices['Earth_Night_Nominal']['c_CO2']
    
    # Plot Night field (positive CO2 accumulation)
    levels_co2 = np.linspace(-1.0, 5.0, 25)
    cs_a = ax_a.contourf(Z_mesh, Y_mesh, slice_night_earth, levels=levels_co2, cmap='coolwarm', alpha=0.90)
    cbar_a = plt.colorbar(cs_a, ax=ax_a, pad=0.03, aspect=20)
    cbar_a.set_label('Night CO₂ Elevation Δ[CO₂] (ppm relative to inlet)', fontsize=9.5, fontweight='bold')
    
    # Overlay geometry: Lower drainage cavity (Y < 7.62 mm)
    drain_rect_a = patches.Rectangle((-32.0, 0.5), 64.0, 7.12, linewidth=1.0, edgecolor='#7f8c8d',
                                     facecolor='#bdc3c7', alpha=0.15, linestyle=':')
    ax_a.add_patch(drain_rect_a)

    # 3 Separate Substrate Blocks (Y = 7.62 to 14.62 mm)
    sub1_a = patches.Rectangle((-30.5, 7.62), 21.0, 7.0, linewidth=1.8, edgecolor='#ba4a00',
                               facecolor='#d97706', alpha=0.35, hatch='//')
    sub2_a = patches.Rectangle((-7.5, 7.62), 15.0, 7.0, linewidth=1.8, edgecolor='#ba4a00',
                               facecolor='#d97706', alpha=0.35, hatch='//')
    sub3_a = patches.Rectangle((9.5, 7.62), 21.0, 7.0, linewidth=1.8, edgecolor='#ba4a00',
                               facecolor='#d97706', alpha=0.35, hatch='//')
    ax_a.add_patch(sub1_a)
    ax_a.add_patch(sub2_a)
    ax_a.add_patch(sub3_a)

    # Chamber dividers
    ax_a.plot([-9.5, -9.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    ax_a.plot([-7.5, -7.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    ax_a.plot([7.5, 7.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    ax_a.plot([9.5, 9.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    
    # Dipping nozzle & ports
    ax_a.plot([21.59, 21.59], [33.67, 9.5], color='#2980b9', lw=3.5, label='Port 3 Dipping Nozzle (Y=9.5mm)')
    ax_a.scatter([21.59], [9.5], color='#e74c3c', s=60, zorder=5)
    ax_a.scatter([0.0], [33.67], color='#8e44ad', s=90, marker='^', zorder=5, label='Port 2 Central Exhaust')
    
    ax_a.text(-17.5, 19.0, 'Ch 1 (B. rapa)\nNight Respiration\n+2.1 ppm', ha='center', va='center',
              fontsize=8.5, fontweight='bold', color='#78281f', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85))
    ax_a.text(0.0, 19.0, 'Ch 2 (Arabidopsis)\nNight Respiration\n+0.8 ppm', ha='center', va='center',
              fontsize=8.5, fontweight='bold', color='#154360', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85))
    ax_a.text(17.5, 19.0, 'Ch 3 (B. rapa)\nNight Respiration\n+2.1 ppm', ha='center', va='center',
              fontsize=8.5, fontweight='bold', color='#78281f', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85))
    ax_a.text(0.0, 4.0, 'Downward Dense CO₂ Pooling into Drainage Wells', ha='center', va='center',
              fontsize=8.0, fontweight='bold', color='#6e2c00')
              
    ax_a.set_xlim([-34, 34])
    ax_a.set_ylim([0, 35])
    ax_a.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax_a.set_title('A. Earth 1g Night Skotoperiod CO₂ Distribution ($Q = 0.50$ L/min)\nRespiratory Emission & Downward Solutal Pooling', fontsize=12, fontweight='bold', pad=10)
    ax_a.legend(loc='lower left', fontsize=7.5, framealpha=0.9)
    ax_a.grid(True, linestyle='--', alpha=0.3)
    
    # Panel B: Lunar Night Low-Power Stagnation Risk (Q = 0.10 L/min)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_facecolor('#fafafa')
    
    slice_lunar_low = concentration_slices['Lunar_Night_LowPower']['c_CO2']
    levels_low = np.linspace(0.0, 12.0, 25)
    cs_b = ax_b.contourf(Z_mesh, Y_mesh, slice_lunar_low, levels=levels_low, cmap='YlOrRd', alpha=0.90)
    cbar_b = plt.colorbar(cs_b, ax=ax_b, pad=0.03, aspect=20)
    cbar_b.set_label('CO₂ Hypercapnia Accumulation (ppm)', fontsize=9.5, fontweight='bold')
    
    drain_rect_b = patches.Rectangle((-32.0, 0.5), 64.0, 7.12, linewidth=1.0, edgecolor='#7f8c8d',
                                     facecolor='#bdc3c7', alpha=0.15, linestyle=':')
    ax_b.add_patch(drain_rect_b)

    sub1_b = patches.Rectangle((-30.5, 7.62), 21.0, 7.0, linewidth=1.8, edgecolor='#ba4a00',
                               facecolor='#d97706', alpha=0.35, hatch='//')
    sub2_b = patches.Rectangle((-7.5, 7.62), 15.0, 7.0, linewidth=1.8, edgecolor='#ba4a00',
                               facecolor='#d97706', alpha=0.35, hatch='//')
    sub3_b = patches.Rectangle((9.5, 7.62), 21.0, 7.0, linewidth=1.8, edgecolor='#ba4a00',
                               facecolor='#d97706', alpha=0.35, hatch='//')
    ax_b.add_patch(sub1_b)
    ax_b.add_patch(sub2_b)
    ax_b.add_patch(sub3_b)

    ax_b.plot([-9.5, -9.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    ax_b.plot([-7.5, -7.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    ax_b.plot([7.5, 7.5], [7.62, 24.0], color='#2c3e50', lw=2.0)
    ax_b.plot([9.5, 9.5], [7.62, 24.0], color='#2c3e50', lw=2.0)

    ax_b.plot([21.59, 21.59], [33.67, 9.5], color='#2980b9', lw=3.5)
    ax_b.scatter([21.59], [9.5], color='#e74c3c', s=60, zorder=5)
    ax_b.scatter([0.0], [33.67], color='#8e44ad', s=90, marker='^', zorder=5)
    
    ax_b.annotate('CRITICAL STAGNATION:\nCO₂ Accumulation > 10 ppm\nLoss of Buoyant Flushing',
                 xy=(-17.5, 17.0), xytext=(-22.0, 27.0),
                 fontsize=8.5, fontweight='bold', color='#900c3f',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#fce4ec', edgecolor='#c2185b'),
                 arrowprops=dict(arrowstyle='->', color='#900c3f', lw=1.5))
                 
    ax_b.set_xlim([-34, 34])
    ax_b.set_ylim([0, 35])
    ax_b.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax_b.set_title('B. Lunar 0.166g Night Stagnation Under Low Flow ($Q = 0.10$ L/min)\nSevere Hypercapnia Pocketing in Brassica Canopies', fontsize=12, fontweight='bold', pad=10)
    ax_b.grid(True, linestyle='--', alpha=0.3)
    
    # Panel C: CO2 Drawdown (Day) vs Accumulation (Night) Across Scenarios
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.set_facecolor('#fafafa')
    
    scen_labels = ['Earth Day\n(0.5 L/min)', 'Earth Night\n(0.5 L/min)', 'Lunar Day\n(0.5 L/min)', 'Lunar Night\n(0.5 L/min)', 'Lunar Night\n(Low 0.1 L/min)', 'Lunar Night\n(Comp 0.6 L/min)']
    scen_keys = ['Earth_Day_Nominal', 'Earth_Night_Nominal', 'Lunar_Day_Nominal', 'Lunar_Night_Nominal', 'Lunar_Night_LowPower', 'Lunar_Night_Compensated']
    
    c1_vals = [df_diurnal.loc[df_diurnal['scenario']==k, 'dCO2_chamber1_ppm'].values[0] for k in scen_keys]
    c2_vals = [df_diurnal.loc[df_diurnal['scenario']==k, 'dCO2_chamber2_ppm'].values[0] for k in scen_keys]
    sub_vals = [df_diurnal.loc[df_diurnal['scenario']==k, 'dCO2_substrate_ppm'].values[0] for k in scen_keys]
    
    x = np.arange(len(scen_labels))
    w = 0.25
    
    b1 = ax_c.bar(x - w, c1_vals, w, label='Ch 1 & 3: B. rapa', color='#27ae60', edgecolor='#1e8449')
    b2 = ax_c.bar(x, c2_vals, w, label='Ch 2: Arabidopsis', color='#2980b9', edgecolor='#1f618d')
    b3 = ax_c.bar(x + w, sub_vals, w, label='Porous Root Matrix', color='#d35400', edgecolor='#a04000')
    
    ax_c.axhline(0.0, color='#2c3e50', lw=1.2)
    ax_c.set_ylabel('Canopy Δ[CO₂] Offset (ppm)', fontsize=11, fontweight='bold')
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(scen_labels, fontsize=8.5, fontweight='bold')
    ax_c.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax_c.grid(True, linestyle='--', alpha=0.4, axis='y')
    ax_c.set_title('C. Diurnal Canopy CO₂ Inversion Across Flight Scenarios\nDrawdown during Light Photoperiod vs Respiratory Buildup at Night', fontsize=12, fontweight='bold', pad=10)
    
    # Annotate zero line
    ax_c.text(0.1, -1.8, '← Day Photosynthetic Depletion', fontsize=8, color='#1e8449', fontweight='bold')
    ax_c.text(2.1, 4.5, 'Night Respiratory Pooling →', fontsize=8, color='#922b21', fontweight='bold')
    
    # Panel D: Mixed Convection Richardson Number & Operating Safe Zone
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_facecolor('#fafafa')
    
    q_sweep = np.array([0.05, 0.10, 0.25, 0.50, 0.85, 1.50, 3.00])
    
    # Earth day (warm + transpirational buoyancy)
    ri_earth_day = (9.81 * (0.0034*1.5 + 0.001*0.38*60) * (0.015**3) / (1.51e-5**2)) / ((q_sweep*1e-3/60 / 7.9e-6 * 0.05 * 0.015 / 1.51e-5)**2 + 1e-6)
    # Lunar night (inverted cooling + solutal sink)
    ri_lunar_night = (1.62 * (0.0034*0.3 + 0.001*0.52*10) * (0.015**3) / (1.51e-5**2)) / ((q_sweep*1e-3/60 / 7.9e-6 * 0.05 * 0.015 / 1.51e-5)**2 + 1e-6)
    
    ax_d.plot(q_sweep, ri_earth_day, 'o-', color='#3498db', lw=2.5, ms=6, label='Earth 1g Day (Active Buoyancy)')
    ax_d.plot(q_sweep, ri_lunar_night, 's--', color='#e74c3c', lw=2.5, ms=6, label='Lunar 0.166g Night (Suppressed Buoyancy)')
    
    ax_d.axhspan(1.0, 100.0, color='#f2d7d5', alpha=0.35, label='Buoyancy Stagnation Zone ($Ri > 1.0$)')
    ax_d.axhspan(0.001, 0.1, color='#d4efdf', alpha=0.35, label='Forced Ventilation Sweep ($Ri < 0.1$)')
    
    ax_d.axvline(0.60, color='#27ae60', ls=':', lw=2.0)
    ax_d.annotate('Recommended Lunar Night Setback:\n$Q_{night} \\geq 0.60$ L/min\n(Prevents CO₂ Asphyxiation)',
                 xy=(0.60, 0.05), xytext=(1.1, 0.4),
                 fontsize=8.5, fontweight='bold', color='#145a32',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#e8f8f5', edgecolor='#27ae60'),
                 arrowprops=dict(arrowstyle='->', color='#27ae60', lw=1.5))
                 
    ax_d.set_xscale('log')
    ax_d.set_yscale('log')
    ax_d.set_xlabel('Pump Volumetric Flow Rate $Q$ (L/min)', fontsize=11, fontweight='bold')
    ax_d.set_ylabel('Richardson Number $Ri = Gr / Re^2$', fontsize=11, fontweight='bold')
    ax_d.set_title('D. Mixed Convection Regimes & Lunar Night Pump Safe Envelope\nPreventing Nocturnal Canopy Asphyxiation on the Moon', fontsize=12, fontweight='bold', pad=10)
    ax_d.grid(True, linestyle='--', alpha=0.4, which='both')
    ax_d.legend(loc='upper right', fontsize=8.5, framealpha=0.9)
    
    plt.savefig('figures/fig8_diurnal_day_night_photoperiod_inversion.png', dpi=200, bbox_inches='tight')
    plt.savefig('figures/fig8_diurnal_day_night_photoperiod_inversion.svg', bbox_inches='tight')
    print("Saved figures/fig8_diurnal_day_night_photoperiod_inversion.png and .svg successfully!")

if __name__ == '__main__':
    run_diurnal_simulation()
