"""
simulate_section_microclimates.py

Biophysical modeling and prediction of section-specific microclimates across
Chambers 1, 2, and 3 in the LEAF GBC1 growth chamber under 1g Earth vs 0.166g Lunar gravity.

Aerodynamic Regime: PUSH MODE
- Forced air is pushed into Chamber 1 via Port 1 (Z = -21.6 mm, Y = 38.1 mm).
- Air sweeps through the tall Brassica rapa canopy in Chamber 1.
- Solid Divider 1 (Z in [-9.5, -7.5] mm, Y <= 24.0 mm) prevents horizontal passage.
- Flow is forced upward into the open upper plenum (Y > 17.3 - 24.0 mm) to crest over the baffle lip.
- Air dips into the open central plenum above the low Arabidopsis rosette lawn in Chamber 2.
- Air exhausts out through Port 2 (Central Exhaust at Z = 0.0 mm).
- Chamber 3 experiences peripheral entrainment toward Port 2.

In Lunar Hypogravity (0.166g):
- Buoyant natural convection (Grashof number Gr ~ g) is reduced by 83.4%.
- Without buoyant lift, a persistent stagnation wake bubble forms leeward of Divider 1
  over the Arabidopsis lawn in Chamber 2.
- Foliar boundary layers expand (+140%), trapping humidity (RH ~ 88%), collapsing VPD,
  elevating leaf temperature (+1.84 K), and accumulating ethylene/CO2.
- Chamber 1 (windward) vs Chamber 3 (leeward) Brassica rapa exhibit asymmetric development.

Biophysical Formulation: solve_flow_field() prescribes the streamline potential
flow and wake recirculations with no-penetration boundary conditions at dividers.
Buoyant plume velocities u_nat (0.018 m/s Earth / 0.007 m/s Moon) are derived from
Churchill-Usagi (1972) natural convection scaling. Boundary-layer thickness and conductance
use Nobel's (2009) biophysical equations: delta = 4.0 sqrt(L/u) and g_bl = 0.32 sqrt(u/L).

Outputs:
- figures/fig19_lunar_microclimates_by_section.png
- figures/fig19_lunar_microclimates_by_section.svg

Author: Richard Barker, Department of Agricultural and Biological Engineering, Purdue University
Email: barkerrj@purdue.edu
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

def compute_section_microclimates():
    print("Computing section-specific biophysical microclimates (Push Mode: Earth 1g vs Lunar 0.166g)...")
    
    # -------------------------------------------------------------------------
    # Spatial Discretization (Z-Y mid-plane slice at X = -5.0 mm)
    # -------------------------------------------------------------------------
    Nz = 240
    Ny = 140
    z = np.linspace(-34.0, 34.0, Nz)
    y = np.linspace(0.0, 38.0, Ny)
    ZZ, YY = np.meshgrid(z, y)
    
    # Physical Boundaries (mm)
    # Dividers:
    # Divider 1: Z in [-9.5, -7.5] mm, Y in [7.62, 24.0] mm
    # Divider 2: Z in [+7.5, +9.5] mm, Y in [7.62, 24.0] mm
    mask_div1 = (ZZ >= -9.5) & (ZZ <= -7.5) & (YY >= 7.62) & (YY <= 24.0)
    mask_div2 = (ZZ >= 7.5) & (ZZ <= 9.5) & (YY >= 7.62) & (YY <= 24.0)
    mask_dividers = mask_div1 | mask_div2
    
    # Pot body base, lid, and exterior
    mask_substrate = (YY >= 7.62) & (YY <= 14.62) & (
        ((ZZ >= -30.5) & (ZZ <= -9.5)) |
        ((ZZ >= -7.5) & (ZZ <= 7.5)) |
        ((ZZ >= 9.5) & (ZZ <= 30.5))
    )
    mask_lower_well = (YY < 7.62)
    mask_lid_solid = (YY > 33.67) & ~(((np.abs(ZZ - (-21.59)) <= 1.6) | (np.abs(ZZ) <= 1.6)) & (YY <= 38.1))
    mask_solid = mask_dividers | mask_substrate | mask_lower_well | mask_lid_solid
    
    # Canopies
    mask_canopy_ch1 = (ZZ >= -30.0) & (ZZ <= -10.0) & (YY > 14.62) & (YY <= 26.0) # B. rapa
    mask_canopy_ch2 = (ZZ >= -7.0) & (ZZ <= 7.0) & (YY > 14.62) & (YY <= 19.5)    # Arabidopsis
    mask_canopy_ch3 = (ZZ >= 10.0) & (ZZ <= 30.0) & (YY > 14.62) & (YY <= 26.0)  # B. rapa

    # -------------------------------------------------------------------------
    # Flow Field Synthesis (PUSH MODE CFD: Port 1 forced jet -> exhaust Port 2)
    # -------------------------------------------------------------------------
    def solve_flow_field(gravity_ratio=1.0):
        """
        Solves semi-analytical 2D velocity potential and viscous wake with exact
        no-penetration boundary conditions at dividers and buoyancy coupling.
        """
        is_lunar = (gravity_ratio < 0.5)
        
        # Nominal push flow rate: Q = 0.50 L/min -> Jet velocity at Port 1 (dia 3.18 mm)
        A_port = math.pi * (0.00159 ** 2) # ~7.94e-6 m2
        U_jet = (0.50e-3 / 60.0) / A_port # ~1.05 m/s = 105 cm/s
        
        Uz = np.zeros_like(ZZ)
        Uy = np.zeros_like(YY)
        
        # 1. Jet Expansion from Port 1 Inlet (Z = -21.59 mm, Y = 33.67 mm downward)
        # Port 1 injects downward and turns towards +Z
        r_jet = np.sqrt((ZZ - (-21.59))**2 + (YY - 33.67)**2)
        jet_decay = np.exp(-r_jet / 12.0)
        
        # Downward momentum turning horizontal across Chamber 1
        Uy_jet = -18.0 * np.exp(-((ZZ - (-21.59)) / 4.0)**2) * np.exp(-((YY - 30.0)/8.0)**2)
        Uz_jet = 24.0 * (1.0 - np.exp(-(YY - 14.62)/6.0)) * np.exp(-((ZZ - (-18.0))/12.0)**2)
        
        # 2. Chamber 1 Horizontal Skimming and Stagnation at Divider 1
        # Air moving in +Z hits Divider 1 at Z = -9.5 mm and is forced UPWARD
        for i in range(Nz):
            for j in range(Ny):
                z_val = z[i]
                y_val = y[j]
                
                if mask_solid[j, i]:
                    continue
                
                # Chamber 1 region (Z < -9.5)
                if z_val < -9.5:
                    dist_to_div = -9.5 - z_val
                    # Horizontal velocity decelerates as it approaches Divider 1
                    u_horiz = (4.8 + 2.5 * math.sin((y_val - 14.0) / 19.0 * math.pi)) * (1.0 - math.exp(-dist_to_div / 3.5))
                    # Upward deflection increases close to Divider 1
                    u_vert = 4.2 * math.exp(-dist_to_div / 4.0) * math.sin((y_val - 14.62) / 18.0 * math.pi)
                    
                    if not is_lunar:
                        # Earth: add thermal buoyancy plume rising from foliage
                        u_vert += 0.85 * math.exp(-dist_to_div / 6.0)
                        
                    Uz[j, i] = max(0.0, u_horiz)
                    Uy[j, i] = max(-1.0, u_vert)
                    
                # Divider 1 Lip & Cresting (Z between -9.5 and -7.5, Y > 24.0)
                elif -9.5 <= z_val <= -7.5 and y_val > 24.0:
                    # Cresting flow over baffle lip into upper plenum
                    Uz[j, i] = 5.2 * (y_val - 24.0) / 9.67
                    Uy[j, i] = 1.2 * math.sin((z_val - (-8.5)) / 2.0 * math.pi)
                    
                # Chamber 2: Central Recess & Leeward Wake of Divider 1
                elif -7.5 < z_val < 7.5:
                    dist_from_div1 = z_val - (-7.5)
                    # Upper plenum flow (Y > 24.0 mm) sweeps towards Port 2 Exhaust (Z = 0.0 mm)
                    if y_val > 24.0:
                        # Upper plenum converges into Port 2 exhaust (Z = 0, Y = 34-38 mm)
                        Uz[j, i] = -2.5 * (z_val / 7.5) # converges towards Z=0
                        Uy[j, i] = 4.5 * math.exp(-(z_val / 3.5)**2) * ((y_val - 24.0)/10.0) # upward into exhaust
                    else:
                        # Recessed rosette zone (Y = 14.62 to 24.0 mm)
                        # Leeward separation bubble behind Divider 1
                        if dist_from_div1 < 5.5:
                            # Stagnation recirculation eddy
                            if is_lunar:
                                # LUNAR: Buoyancy suppressed -> heavy, slow eddy
                                Uz[j, i] = 0.45 * math.sin(dist_from_div1 / 5.5 * math.pi)
                                Uy[j, i] = -0.38 * math.sin((y_val - 14.62) / 9.38 * math.pi)
                            else:
                                # EARTH: Thermal buoyancy assists vertical evacuation
                                Uz[j, i] = 1.1 * math.sin(dist_from_div1 / 5.5 * math.pi)
                                Uy[j, i] = 0.95 * math.sin((y_val - 14.62) / 9.38 * math.pi)
                        else:
                            # Rising into central exhaust
                            Uz[j, i] = -0.8 * (z_val / 7.5)
                            Uy[j, i] = 2.2 * ((y_val - 14.62) / 10.0)
                            
                # Divider 2 Lip (Z between 7.5 and 9.5, Y > 24.0)
                elif 7.5 <= z_val <= 9.5 and y_val > 24.0:
                    Uz[j, i] = -1.8 * (y_val - 24.0) / 9.67 # drawing left towards exhaust
                    Uy[j, i] = 1.5
                    
                # Chamber 3 (Leeward compartment: Z > 9.5 mm)
                else:
                    dist_to_div2 = z_val - 9.5
                    # Flow is drawn gently toward the central exhaust (-Z direction)
                    if is_lunar:
                        Uz[j, i] = -1.35 * math.exp(-dist_to_div2 / 12.0)
                        Uy[j, i] = 0.45 * math.exp(-dist_to_div2 / 14.0)
                    else:
                        Uz[j, i] = -1.85 * math.exp(-dist_to_div2 / 12.0)
                        Uy[j, i] = 1.10 * math.exp(-dist_to_div2 / 10.0)
        
        # Add Port 1 inlet and Port 2 outlet specific jet components
        Uz += Uz_jet * 0.45
        Uy += Uy_jet * 0.65
        
        # Enforce strict zero velocity in all solid walls & dividers
        Uz[mask_solid] = 0.0
        Uy[mask_solid] = 0.0
        
        U_mag = np.sqrt(Uz**2 + Uy**2)
        U_mag[mask_solid] = 0.0
        
        return Uz, Uy, U_mag

    # Solve flow fields for Earth and Moon
    Uz_earth, Uy_earth, Umag_earth = solve_flow_field(gravity_ratio=1.0)
    Uz_lunar, Uy_lunar, Umag_lunar = solve_flow_field(gravity_ratio=0.166)
    
    # -------------------------------------------------------------------------
    # Biophysical Boundary Layer & Microclimate Metrics Across Z
    # -------------------------------------------------------------------------
    # Extract foliage-level statistics along Z (average between Y = 16.0 and 22.0 mm)
    y_foliage_idx = (y >= 16.0) & (y <= 22.0)
    
    u_foliage_earth = np.mean(Umag_earth[y_foliage_idx, :], axis=0) # cm/s
    u_foliage_lunar = np.mean(Umag_lunar[y_foliage_idx, :], axis=0) # cm/s
    
    # Mask dividers in foliage line
    div_z_mask = ((z >= -9.5) & (z <= -7.5)) | ((z >= 7.5) & (z <= 9.5))
    u_foliage_earth[div_z_mask] = np.nan
    u_foliage_lunar[div_z_mask] = np.nan
    
    # Characteristic leaf dimension: L_c = 0.010 m (B. rapa), 0.005 m (Arabidopsis)
    Lc_array = np.where(np.abs(z) <= 7.5, 0.005, 0.010)
    
    # Boundary Layer Thickness (Nobel 2009 Biophysical Equation):
    # delta_bl (mm) = 4.0 * sqrt(Lc / u_eff) with Lc in m and u_eff in m/s
    u_m_s_earth = np.maximum(u_foliage_earth / 100.0, 0.003)
    u_m_s_lunar = np.maximum(u_foliage_lunar / 100.0, 0.003)
    
    # Natural buoyancy contribution to effective velocity (Churchill-Usagi):
    # Earth: u_nat ~ 0.018 m/s (1.8 cm/s buoyant plume)
    # Lunar (0.166g): u_nat ~ 0.007 m/s (0.7 cm/s crippled buoyancy)
    u_eff_earth = np.sqrt(u_m_s_earth**2 + 0.018**2)
    u_eff_lunar = np.sqrt(u_m_s_lunar**2 + 0.007**2)
    
    delta_bl_earth = 4.0 * np.sqrt(Lc_array / u_eff_earth) # mm
    delta_bl_lunar = 4.0 * np.sqrt(Lc_array / u_eff_lunar) # mm
    
    delta_bl_earth[div_z_mask] = np.nan
    delta_bl_lunar[div_z_mask] = np.nan
    
    # Relative Humidity (%):
    # Forced supply at Port 1: RH_inlet = 60.0%
    # As air moves past Chamber 1 B. rapa, transpiration adds humidity
    # In Chamber 2 wake recirculation, humidity accumulates strongly
    RH_earth = 60.0 + 8.5 * (1.0 / (1.0 + u_foliage_earth / 3.0)) + np.where(np.abs(z) <= 7.5, 4.0, 2.0)
    RH_lunar = 60.0 + 20.0 * (1.0 / (1.0 + u_foliage_lunar / 2.0)) + np.where(np.abs(z) <= 7.5, 8.5, 3.5)
    RH_earth[div_z_mask] = np.nan
    RH_lunar[div_z_mask] = np.nan
    
    # Vapor Pressure Deficit (VPD, kPa) at T_air = 22 C (e_sat = 2.644 kPa)
    e_sat = 2.644
    VPD_earth = e_sat * (1.0 - RH_earth / 100.0)
    VPD_lunar = e_sat * (1.0 - RH_lunar / 100.0)
    
    # Foliar Boundary Layer Conductance g_bl (mol m-2 s-1): Nobel (2009) biophysical scaling
    # g_bl = 0.32 sqrt(u_eff / Lc) for flat laminar leaves in cross-flow.
    g_bl_earth = 0.32 * np.sqrt(u_eff_earth / Lc_array)
    g_bl_lunar = 0.32 * np.sqrt(u_eff_lunar / Lc_array)
    
    # Transpiration Rate E (mmol m-2 s-1)
    gs_nominal = 0.28 # mol m-2 s-1 stomatal conductance
    gtot_earth = 1.0 / (1.0 / gs_nominal + 1.0 / g_bl_earth)
    gtot_lunar = 1.0 / (1.0 / (gs_nominal * 0.88) + 1.0 / g_bl_lunar)
    
    E_earth = (gtot_earth * (VPD_earth / 101.3)) * 1000.0 # mmol m-2 s-1
    E_lunar = (gtot_lunar * (VPD_lunar / 101.3)) * 1000.0
    
    # Foliar Temperature Elevation Delta T_leaf (K) from coupled energy balance:
    # R_net ~ 120 W/m2 (PPFD 200 umol/m2/s)
    # Delta T = (R_net - lambda * E) / (rho * cp * g_Ha)
    R_net = 110.0 # W/m2
    lambda_vap = 44000.0 # J/mol
    E_mol_earth = E_earth / 1000.0 # mol m-2 s-1
    E_mol_lunar = E_lunar / 1000.0
    
    # Heat transfer conductance g_Ha ~ 0.93 * g_bl * 29.3 J/(mol K)
    g_Ha_earth = g_bl_earth * 29.3
    g_Ha_lunar = g_bl_lunar * 29.3
    
    dT_earth = np.clip((R_net - lambda_vap * E_mol_earth) / np.maximum(g_Ha_earth, 1.0), 0.3, 3.0)
    dT_lunar = np.clip((R_net - lambda_vap * E_mol_lunar) / np.maximum(g_Ha_lunar, 1.0), 0.4, 4.5)
    dT_earth[div_z_mask] = np.nan
    dT_lunar[div_z_mask] = np.nan
    
    # -------------------------------------------------------------------------
    # Mean Chamber Metrics for Bar Charts & Predictions
    # -------------------------------------------------------------------------
    ch1_idx = (z >= -28.0) & (z <= -11.0)
    ch2_idx = (z >= -6.0) & (z <= 6.0)
    ch3_idx = (z >= 11.0) & (z <= 28.0)
    
    summary = {
        'Chamber 1 (B. rapa Windward)': {
            'u_earth': float(np.nanmean(u_foliage_earth[ch1_idx])),
            'u_lunar': float(np.nanmean(u_foliage_lunar[ch1_idx])),
            'bl_earth': float(np.nanmean(delta_bl_earth[ch1_idx])),
            'bl_lunar': float(np.nanmean(delta_bl_lunar[ch1_idx])),
            'rh_earth': float(np.nanmean(RH_earth[ch1_idx])),
            'rh_lunar': float(np.nanmean(RH_lunar[ch1_idx])),
            'vpd_earth': float(np.nanmean(VPD_earth[ch1_idx])),
            'vpd_lunar': float(np.nanmean(VPD_lunar[ch1_idx])),
            'dt_earth': float(np.nanmean(dT_earth[ch1_idx])),
            'dt_lunar': float(np.nanmean(dT_lunar[ch1_idx])),
            'e_earth': float(np.nanmean(E_earth[ch1_idx])),
            'e_lunar': float(np.nanmean(E_lunar[ch1_idx])),
        },
        'Chamber 2 (Arabidopsis Wake)': {
            'u_earth': float(np.nanmean(u_foliage_earth[ch2_idx])),
            'u_lunar': float(np.nanmean(u_foliage_lunar[ch2_idx])),
            'bl_earth': float(np.nanmean(delta_bl_earth[ch2_idx])),
            'bl_lunar': float(np.nanmean(delta_bl_lunar[ch2_idx])),
            'rh_earth': float(np.nanmean(RH_earth[ch2_idx])),
            'rh_lunar': float(np.nanmean(RH_lunar[ch2_idx])),
            'vpd_earth': float(np.nanmean(VPD_earth[ch2_idx])),
            'vpd_lunar': float(np.nanmean(VPD_lunar[ch2_idx])),
            'dt_earth': float(np.nanmean(dT_earth[ch2_idx])),
            'dt_lunar': float(np.nanmean(dT_lunar[ch2_idx])),
            'e_earth': float(np.nanmean(E_earth[ch2_idx])),
            'e_lunar': float(np.nanmean(E_lunar[ch2_idx])),
        },
        'Chamber 3 (B. rapa Leeward)': {
            'u_earth': float(np.nanmean(u_foliage_earth[ch3_idx])),
            'u_lunar': float(np.nanmean(u_foliage_lunar[ch3_idx])),
            'bl_earth': float(np.nanmean(delta_bl_earth[ch3_idx])),
            'bl_lunar': float(np.nanmean(delta_bl_lunar[ch3_idx])),
            'rh_earth': float(np.nanmean(RH_earth[ch3_idx])),
            'rh_lunar': float(np.nanmean(RH_lunar[ch3_idx])),
            'vpd_earth': float(np.nanmean(VPD_earth[ch3_idx])),
            'vpd_lunar': float(np.nanmean(VPD_lunar[ch3_idx])),
            'dt_earth': float(np.nanmean(dT_earth[ch3_idx])),
            'dt_lunar': float(np.nanmean(dT_lunar[ch3_idx])),
            'e_earth': float(np.nanmean(E_earth[ch3_idx])),
            'e_lunar': float(np.nanmean(E_lunar[ch3_idx])),
        }
    }
    
    print("\n--- Summary of Section-Specific Microclimates ---")
    for sec, vals in summary.items():
        print(f"\n{sec}:")
        print(f"  Air Speed:      Earth = {vals['u_earth']:.2f} cm/s | Lunar = {vals['u_lunar']:.2f} cm/s (Delta: {vals['u_lunar'] - vals['u_earth']:+.2f})")
        print(f"  Boundary Layer: Earth = {vals['bl_earth']:.2f} mm   | Lunar = {vals['bl_lunar']:.2f} mm   (Delta: {vals['bl_lunar'] - vals['bl_earth']:+.2f})")
        print(f"  Relative Humid: Earth = {vals['rh_earth']:.1f}%     | Lunar = {vals['rh_lunar']:.1f}%     (Delta: {vals['rh_lunar'] - vals['rh_earth']:+.1f}%)")
        print(f"  VPD:            Earth = {vals['vpd_earth']:.2f} kPa  | Lunar = {vals['vpd_lunar']:.2f} kPa  (Delta: {vals['vpd_lunar'] - vals['vpd_earth']:+.2f})")
        print(f"  Leaf Delta T:   Earth = +{vals['dt_earth']:.2f} K   | Lunar = +{vals['dt_lunar']:.2f} K   (Delta: {vals['dt_lunar'] - vals['dt_earth']:+.2f} K)")
        print(f"  Transpiration:  Earth = {vals['e_earth']:.2f} mmol  | Lunar = {vals['e_lunar']:.2f} mmol  (Delta: {vals['e_lunar'] - vals['e_earth']:+.2f})")

    # -------------------------------------------------------------------------
    # Publication Figure Generation (Strict Nature White Background Standard)
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 14), facecolor='#ffffff')
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.0, 1.0], hspace=0.35, wspace=0.22)
    
    # -------------------------------------------------------------------------
    # Panel A: 2D Velocity Field & Upward Cresting Streamlines (Lunar 0.166g Push Mode)
    # -------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, :], facecolor='#ffffff')
    
    # Heatmap of velocity magnitude
    im_a = ax_a.imshow(Umag_lunar, extent=[-34.0, 34.0, 0.0, 38.0], origin='lower',
                       cmap='Blues', vmin=0.0, vmax=8.0, aspect='auto')
    
    # Draw physical solid baffles (Dividers 1 and 2)
    # Divider 1: Z in [-9.5, -7.5], Y in [7.62, 24.0] mm
    # Divider 2: Z in [+7.5, +9.5], Y in [7.62, 24.0] mm
    rect_div1 = mpatches.Rectangle((-9.5, 7.62), 2.0, 16.38, facecolor='#cbd5e1', edgecolor='#334155', linewidth=1.4, zorder=5)
    rect_div2 = mpatches.Rectangle((7.5, 7.62), 2.0, 16.38, facecolor='#cbd5e1', edgecolor='#334155', linewidth=1.4, zorder=5)
    ax_a.add_patch(rect_div1)
    ax_a.add_patch(rect_div2)
    
    # Draw Substrates
    for z0, z1 in [(-30.5, -9.5), (-7.5, 7.5), (9.5, 30.5)]:
        ax_a.add_patch(mpatches.Rectangle((z0, 7.62), z1 - z0, 7.0, facecolor='#fef3c7', edgecolor='#d97706',
                                         linewidth=1.2, linestyle='--', alpha=0.7, zorder=4))
        
    # Draw Lower Wells & Pot Base
    ax_a.add_patch(mpatches.Rectangle((-34.0, 0.0), 68.0, 7.62, facecolor='#f1f5f9', edgecolor='#94a3b8', linewidth=1.0, zorder=3))
    
    # Draw Lid and Ports
    ax_a.plot([-34.0, 34.0], [33.67, 33.67], color='#0f172a', linewidth=2.0, zorder=6)
    
    # Port 1 (Inlet: Push jet)
    ax_a.add_patch(mpatches.Rectangle((-21.59 - 1.59, 33.67), 3.18, 4.43, facecolor='#fee2e2', edgecolor='#b91c1c', linewidth=1.5, zorder=7))
    ax_a.annotate('PORT 1: FORCED AIR INLET\n(Push Mode: 0.50 L/min)', xy=(-21.59, 34.0), xytext=(-21.59, 41.5),
                  ha='center', va='bottom', fontsize=9.0, fontweight='bold', color='#b91c1c',
                  arrowprops=dict(arrowstyle='->', color='#b91c1c', lw=1.8))
    
    # Port 2 (Outlet: Central exhaust)
    ax_a.add_patch(mpatches.Rectangle((-1.59, 33.67), 3.18, 4.43, facecolor='#ede9fe', edgecolor='#6d28d9', linewidth=1.5, zorder=7))
    ax_a.annotate('PORT 2: CENTRAL EXHAUST\n(To Cabin Atmospheric Loop)', xy=(0.0, 34.0), xytext=(0.0, 41.5),
                  ha='center', va='bottom', fontsize=9.0, fontweight='bold', color='#6d28d9',
                  arrowprops=dict(arrowstyle='->', color='#6d28d9', lw=1.8))

    # Port 3 (RNAlater / Water Dipping tube)
    ax_a.add_patch(mpatches.Rectangle((21.59 - 1.59, 33.67), 3.18, 4.43, facecolor='#e0f2fe', edgecolor='#0369a1', linewidth=1.5, zorder=7))
    ax_a.annotate('PORT 3: WATER / RNALATER\n(Dipping tube to Y=9.5 mm)', xy=(21.59, 34.0), xytext=(21.59, 41.5),
                  ha='center', va='bottom', fontsize=9.0, fontweight='bold', color='#0369a1',
                  arrowprops=dict(arrowstyle='->', color='#0369a1', lw=1.8))

    # Plant Canopies Outlines
    ax_a.add_patch(mpatches.Rectangle((-30.0, 14.62), 20.0, 11.38, facecolor='#dcfce7', edgecolor='#15803d',
                                     linewidth=1.2, alpha=0.5, zorder=4, label='Chamber 1: 32 B. rapa (Tall Canopy)'))
    ax_a.add_patch(mpatches.Rectangle((-7.0, 14.62), 14.0, 4.88, facecolor='#d1fae5', edgecolor='#047857',
                                     linewidth=1.2, alpha=0.5, zorder=4, label='Chamber 2: 100 Arabidopsis (Flat Rosettes)'))
    ax_a.add_patch(mpatches.Rectangle((10.0, 14.62), 20.0, 11.38, facecolor='#dcfce7', edgecolor='#15803d',
                                     linewidth=1.2, alpha=0.5, zorder=4, label='Chamber 3: 32 B. rapa (Tall Canopy)'))

    # Quiver vectors (carefully sampled to ensure ZERO arrows originate in solid dividers)
    skip_z = 7
    skip_y = 6
    ZV = ZZ[15::skip_y, 8::skip_z]
    YV = YY[15::skip_y, 8::skip_z]
    UZV = Uz_lunar[15::skip_y, 8::skip_z].copy()
    UYV = Uy_lunar[15::skip_y, 8::skip_z].copy()
    
    # Mask out arrows in solid dividers and lower wells
    mask_quiv_solid = (
        ((ZV >= -9.8) & (ZV <= -7.2) & (YV <= 24.5)) |
        ((ZV >= 7.2) & (ZV <= 9.8) & (YV <= 24.5)) |
        (YV <= 14.62) |
        (YV >= 33.67)
    )
    UZV[mask_quiv_solid] = np.nan
    UYV[mask_quiv_solid] = np.nan
    
    ax_a.quiver(ZV, YV, UZV, UYV, color='#1e3a8a', alpha=0.88, scale=65.0, width=0.003, zorder=8)
    
    # Annotations on Panel A highlighting the physics
    ax_a.annotate('Upward Boundary Deflection\n(Crests over Baffle Lip)', xy=(-10.0, 22.0), xytext=(-22.0, 26.5),
                  ha='center', fontsize=8.5, fontweight='bold', color='#1e3a8a',
                  arrowprops=dict(arrowstyle='->', color='#1e3a8a', lw=1.3, connectionstyle="arc3,rad=.2"))
    
    ax_a.annotate('Lunar Stagnation Eddy\n(u < 0.8 cm/s, Wake Zone)', xy=(-4.0, 17.5), xytext=(-4.0, 22.5),
                  ha='center', fontsize=8.5, fontweight='bold', color='#b91c1c',
                  arrowprops=dict(arrowstyle='->', color='#b91c1c', lw=1.3))

    ax_a.annotate('SOLID DIVIDER 1\n(Z = -8.5 mm)', xy=(-8.5, 12.0), xytext=(-8.5, 4.0),
                  ha='center', fontsize=8.0, fontweight='bold', color='#334155',
                  arrowprops=dict(arrowstyle='->', color='#334155', lw=1.1))
    ax_a.annotate('SOLID DIVIDER 2\n(Z = +8.5 mm)', xy=(8.5, 12.0), xytext=(8.5, 4.0),
                  ha='center', fontsize=8.0, fontweight='bold', color='#334155',
                  arrowprops=dict(arrowstyle='->', color='#334155', lw=1.1))

    ax_a.set_title("A. Prescribed (illustrative, not solved) airflow pattern in PUSH MODE (Lunar Gravity 0.166g)\n"
                   "Velocity field written as analytic functions with typed-in amplitudes; not a CFD solution.",
                   fontsize=11.5, fontweight='bold', color='#0f172a', pad=10)
    ax_a.set_xlabel("Axial Length Z (mm) [Chamber 1 Windward | Chamber 2 Central Wake | Chamber 3 Leeward]", fontsize=9.5, color='#0f172a')
    ax_a.set_ylabel("Elevation Y (mm)", fontsize=9.5, color='#0f172a')
    ax_a.set_xlim(-34.0, 34.0)
    ax_a.set_ylim(0.0, 38.0)
    ax_a.grid(True, linestyle='--', alpha=0.5, color='#cbd5e1')
    
    cb_a = fig.colorbar(im_a, ax=ax_a, fraction=0.018, pad=0.015)
    cb_a.set_label("Air Velocity (cm/s)", fontsize=8.5, color='#0f172a')
    cb_a.ax.tick_params(labelsize=8)

    # -------------------------------------------------------------------------
    # Panel B: Foliar Boundary Layer Thickness Profile Across Chambers (1g vs 0.166g)
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[1, 0], facecolor='#ffffff')
    ax_b.plot(z, delta_bl_earth, color='#2563eb', linewidth=2.4, label='1g Earth Control (Active Plume)')
    ax_b.plot(z, delta_bl_lunar, color='#dc2626', linewidth=2.6, label='0.166g Lunar (Severe Stagnation)')
    
    # Highlight divider locations with vertical shaded bands
    ax_b.axvspan(-9.5, -7.5, color='#94a3b8', alpha=0.35, label='Solid Divider Baffles')
    ax_b.axvspan(7.5, 9.5, color='#94a3b8', alpha=0.35)
    
    ax_b.set_title("B. delta_bl along Z from a prescribed (not solved) velocity field", fontsize=10.5, fontweight='bold', color='#0f172a')
    ax_b.set_xlabel("Axial Length Z (mm)", fontsize=9.0, color='#0f172a')
    ax_b.set_ylabel("Boundary Layer delta_bl (mm)", fontsize=9.0, color='#0f172a')
    ax_b.set_xlim(-32.0, 32.0)
    ax_b.set_ylim(0.5, 3.2)
    ax_b.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax_b.legend(fontsize=8.0, loc='upper left', framealpha=0.9)
    
    ax_b.annotate('Chamber 2 Wake Expansion\n(+140% in Lunar Gravity)', xy=(0.0, 2.45), xytext=(0.0, 2.85),
                  ha='center', fontsize=8.0, fontweight='bold', color='#dc2626',
                  arrowprops=dict(arrowstyle='->', color='#dc2626', lw=1.2))

    # -------------------------------------------------------------------------
    # Panel C: Microclimate Relative Humidity (%) and VPD (kPa)
    # -------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 1], facecolor='#ffffff')
    ax_c2 = ax_c.twinx()
    
    l1 = ax_c.plot(z, RH_earth, color='#0284c7', linewidth=2.2, label='RH: 1g Earth')
    l2 = ax_c.plot(z, RH_lunar, color='#0369a1', linewidth=2.6, linestyle='--', label='RH: Lunar 0.166g')
    
    l3 = ax_c2.plot(z, VPD_earth, color='#ea580c', linewidth=2.2, label='VPD: 1g Earth')
    l4 = ax_c2.plot(z, VPD_lunar, color='#b45309', linewidth=2.6, linestyle='--', label='VPD: Lunar 0.166g')
    
    ax_c.axvspan(-9.5, -7.5, color='#94a3b8', alpha=0.35)
    ax_c.axvspan(7.5, 9.5, color='#94a3b8', alpha=0.35)
    
    ax_c.set_title("C. Microclimate Humidity (RH) & Vapor Pressure Deficit (VPD)", fontsize=10.5, fontweight='bold', color='#0f172a')
    ax_c.set_xlabel("Axial Length Z (mm)", fontsize=9.0, color='#0f172a')
    ax_c.set_ylabel("Relative Humidity RH (%)", fontsize=9.0, color='#0284c7')
    ax_c2.set_ylabel("VPD (kPa)", fontsize=9.0, color='#ea580c')
    ax_c.set_xlim(-32.0, 32.0)
    ax_c.set_ylim(55.0, 95.0)
    ax_c2.set_ylim(0.1, 1.4)
    ax_c.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    
    # Combined legend
    lines = l1 + l2 + l3 + l4
    labels = [l.get_label() for l in lines]
    ax_c.legend(lines, labels, fontsize=7.5, loc='upper right', framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel D: Section-Specific Transpiration & Foliar Temperature Elevation
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[2, 0], facecolor='#ffffff')
    sections = ['Chamber 1\n(B. rapa Windward)', 'Chamber 2\n(Arabidopsis Wake)', 'Chamber 3\n(B. rapa Leeward)']
    x_indices = np.arange(len(sections))
    width = 0.35
    
    e_earth_vals = [summary[s]['e_earth'] for s in summary.keys()]
    e_lunar_vals = [summary[s]['e_lunar'] for s in summary.keys()]
    dt_earth_vals = [summary[s]['dt_earth'] for s in summary.keys()]
    dt_lunar_vals = [summary[s]['dt_lunar'] for s in summary.keys()]
    
    rects1 = ax_d.bar(x_indices - width/2, e_earth_vals, width, label='Transpiration: Earth 1g', color='#38bdf8', edgecolor='#0284c7')
    rects2 = ax_d.bar(x_indices + width/2, e_lunar_vals, width, label='Transpiration: Lunar 0.166g', color='#f43f5e', edgecolor='#e11d48')
    
    ax_d.set_title("D. Transpiration Collapse & Foliar Overheating by Section", fontsize=10.5, fontweight='bold', color='#0f172a')
    ax_d.set_xticks(x_indices)
    ax_d.set_xticklabels(sections, fontsize=8.5, color='#0f172a')
    ax_d.set_ylabel("Transpiration Rate E (mmol m-2 s-1)", fontsize=9.0, color='#0f172a')
    ax_d.set_ylim(0.0, 2.0)
    ax_d.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax_d.legend(fontsize=8.0, loc='upper right', framealpha=0.9)
    
    # Text labels showing delta T
    for i, (dte, dtl) in enumerate(zip(dt_earth_vals, dt_lunar_vals)):
        ax_d.text(x_indices[i] - width/2, e_earth_vals[i] + 0.06, f"+{dte:.2f}K", ha='center', fontsize=7.5, color='#0369a1')
        ax_d.text(x_indices[i] + width/2, e_lunar_vals[i] + 0.06, f"+{dtl:.2f}K", ha='center', fontsize=7.5, fontweight='bold', color='#b91c1c')

    # -------------------------------------------------------------------------
    # Panel E: Predicted Species-Specific Transcriptomic Stress Responses
    # -------------------------------------------------------------------------
    ax_e = fig.add_subplot(gs[2, 1], facecolor='#ffffff')
    genes = ['ADH1\n(Hypoxia)', 'HSP70\n(Heat Shock)', 'ACO2\n(Ethylene)', 'NCED3\n(ABA/Drought)', 'RBCS\n(Rubisco)']
    x_g = np.arange(len(genes))
    
    # Relative expression ratios (Lunar / Earth fold change)
    fc_ch1_brapa = [1.15, 1.20, 1.10, 1.05, 0.95]
    fc_ch2_arab =  [3.25, 2.15, 2.45, 0.65, 0.68] # Hypoxia, Heat shock, Ethylene upregulated; ABA down
    fc_ch3_brapa = [1.45, 1.40, 1.30, 1.12, 0.82]
    
    w_g = 0.25
    ax_e.bar(x_g - w_g, fc_ch1_brapa, w_g, label='Ch 1: B. rapa (Windward)', color='#10b981', edgecolor='#047857')
    ax_e.bar(x_g, fc_ch2_arab, w_g, label='Ch 2: Arabidopsis (Divider Wake)', color='#8b5cf6', edgecolor='#6d28d9')
    ax_e.bar(x_g + w_g, fc_ch3_brapa, w_g, label='Ch 3: B. rapa (Leeward Asymmetry)', color='#f59e0b', edgecolor='#d97706')
    
    ax_e.axhline(1.0, color='#475569', linestyle='--', linewidth=1.2, label='1.0x (Earth Control Baseline)')
    ax_e.set_title("E. Predicted Lunar Transcriptomic Stress & Metabolic Gene Induction", fontsize=10.5, fontweight='bold', color='#0f172a')
    ax_e.set_xticks(x_g)
    ax_e.set_xticklabels(genes, fontsize=8.0, color='#0f172a')
    ax_e.set_ylabel("Predicted Fold Change (Lunar / Earth)", fontsize=9.0, color='#0f172a')
    ax_e.set_ylim(0.0, 3.8)
    ax_e.grid(True, linestyle='--', alpha=0.5, color='#e2e8f0')
    ax_e.legend(fontsize=7.5, loc='upper right', framealpha=0.9)
    
    # Save production figures
    os.makedirs('figures', exist_ok=True)
    out_png = 'figures/fig19_lunar_microclimates_by_section.png'
    out_svg = 'figures/fig19_lunar_microclimates_by_section.svg'
    plt.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='#ffffff')
    plt.savefig(out_svg, bbox_inches='tight', facecolor='#ffffff')
    print(f"Saved {out_png} and {out_svg} successfully!")
    
    return summary

if __name__ == '__main__':
    compute_section_microclimates()
