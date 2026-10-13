"""
14-Day Developmental Growth Time-Series CFD & Biophysical Model for LEAF GBC2.

Models Day 0 (Seed Imbibition) through Day 14 (Harvest Canopy Closure) for:
  - Chambers 1 & 3: 32 Brassica rapa seeds per chamber (fast-growing erect canopy, h -> 11.4 mm, A -> 64.0 cm2)
  - Chamber 2: 100 Arabidopsis thaliana seeds (high-density miniature rosette lawn, h -> 4.9 mm, A -> 15.0 cm2)

Couples:
  1. Logistic canopy stature h(t), leaf area A(t), LAI(t), and Brinkman drag D_canopy(t)
  2. Autotrophic transition from dark seed respiration (Days 0-2) to FvCB net carbon fixation (Days 3-14)
  3. Root biomass expansion in the 3-part porous substrate (Y = 7.62-14.62 mm), reducing porosity eps(t)
     and Kozeny-Carman permeability K(t), affecting root aeration and harvest RNAlater injection pressure
  4. Comparison of Fixed Flow (Q = 0.50 L/min) vs. Dynamic Adaptive Lunar Schedule Q(t) (0.25 -> 0.85 L/min)
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sim.gbc2_cfd_engine import GBC2CFDSimulator
from sim.advanced_analytics import GBC2AdvancedAnalytics


def run_developmental_timeseries():
    os.makedirs('results/tables', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    print("=================================================================")
    print("Running 14-Day Developmental Growth Time-Series Model (Day 0-14)")
    print("=================================================================")

    sim = GBC2CFDSimulator()
    analytics = GBC2AdvancedAnalytics()

    # Floor footprints (cm2)
    A_floor_brapa_cm2 = 6.783   # per chamber (Ch 1 and Ch 3)
    A_floor_arab_cm2 = 4.845    # Ch 2

    # Substrate initial parameters
    eps_0 = 0.50
    K_0 = 1.0e-10  # m2 (Darcy permeability of fresh substrate)
    V_sub_brapa_mL = 4.02       # Bulk volume per chamber (Ch 1, Ch 3)
    V_sub_arab_mL = 2.51        # Bulk volume Ch 2

    days = np.arange(0, 15, 1)

    # Pre-compute reference CFD velocity scalings at representative flow rates
    # and canopy drag states so the 15-day x 3-regime matrix is physically grounded
    ref_flows = [0.25, 0.50, 0.65, 0.85]
    cfd_cache = {}
    for g_val in [9.81, 1.62]:
        for q_val in ref_flows:
            res, _ = sim.run_case(Q_L_min=q_val, gravity_g=g_val, max_steps=220)
            cfd_cache[(g_val, q_val)] = res

    records = []

    regimes = [
        ('Earth_1g_Nominal', 9.81, 'fixed_0.50'),
        ('Lunar_0.166g_Uncompensated', 1.62, 'fixed_0.50'),
        ('Lunar_0.166g_Adaptive', 1.62, 'adaptive_ramp'),
    ]

    for regime_name, g_val, flow_mode in regimes:
        cum_carbon_brapa_mg = 0.0
        cum_carbon_arab_mg = 0.0

        for day in days:
            # 1. Ontogenetic morphometrics (logistic canopy expansion)
            if day < 2:
                # Seed imbibition & radicle emergence (no expanded cotyledons yet)
                f_brapa = 0.0
                f_arab = 0.0
                h_brapa_mm = 0.5 * (day / 2.0)
                h_arab_mm = 0.3 * (day / 2.0)
            else:
                # Logistic foliar expansion from Day 2 to Day 14
                f_brapa = 1.0 / (1.0 + np.exp(-0.52 * (day - 7.2)))
                f_brapa_norm = (f_brapa - (1.0 / (1.0 + np.exp(-0.52 * (1.0 - 7.2))))) / (
                    (1.0 / (1.0 + np.exp(-0.52 * (14.0 - 7.2)))) - (1.0 / (1.0 + np.exp(-0.52 * (1.0 - 7.2))))
                )
                f_brapa = float(np.clip(f_brapa_norm, 0.01, 1.0))

                f_arab = 1.0 / (1.0 + np.exp(-0.48 * (day - 7.8)))
                f_arab_norm = (f_arab - (1.0 / (1.0 + np.exp(-0.48 * (1.0 - 7.8))))) / (
                    (1.0 / (1.0 + np.exp(-0.48 * (14.0 - 7.8)))) - (1.0 / (1.0 + np.exp(-0.48 * (1.0 - 7.8))))
                )
                f_arab = float(np.clip(f_arab_norm, 0.01, 1.0))

                h_brapa_mm = 0.8 + 10.6 * (f_brapa ** 0.82)
                h_arab_mm = 0.5 + 4.4 * (f_arab ** 0.85)

            A_brapa_cm2 = 64.0 * f_brapa     # per chamber (32 seeds)
            A_arab_cm2 = 15.0 * f_arab       # Ch 2 (100 seeds)

            LAI_brapa = A_brapa_cm2 / A_floor_brapa_cm2
            LAI_arab = A_arab_cm2 / A_floor_arab_cm2

            # Characteristic leaf dimension L_leaf (m) grows from 2 mm cotyledon to 20 mm (B. rapa) / 8 mm (Arabidopsis)
            L_leaf_brapa_m = 0.002 + 0.018 * (f_brapa ** 0.6)
            L_leaf_arab_m = 0.0015 + 0.0065 * (f_arab ** 0.6)

            # 2. Root expansion & substrate porosity compression (Kozeny-Carman)
            f_root = 1.0 / (1.0 + np.exp(-0.45 * (day - 7.5)))
            f_root = float(np.clip((f_root - 0.03) / (0.95 - 0.03), 0.0, 1.0))

            eps_brapa = eps_0 - 0.16 * f_root   # 0.50 -> 0.34 in Ch 1 & 3
            eps_arab = eps_0 - 0.09 * f_root    # 0.50 -> 0.41 in Ch 2
            eps_mean = (2.0 * V_sub_brapa_mL * eps_brapa + V_sub_arab_mL * eps_arab) / (2.0 * V_sub_brapa_mL + V_sub_arab_mL)

            # Kozeny-Carman permeability ratio K(eps) / K0
            K_ratio_brapa = ((eps_brapa ** 3) * ((1.0 - eps_0) ** 2)) / (((eps_0) ** 3) * ((1.0 - eps_brapa) ** 2))
            K_brapa_m2 = K_0 * K_ratio_brapa

            # RNAlater harvest injection pressure if harvested on this day (Q_inj = 10 mL/min)
            dP_rnalater_kPa = 0.73 / K_ratio_brapa
            V_void_total_mL = (2.0 * V_sub_brapa_mL * eps_brapa) + (V_sub_arab_mL * eps_arab)
            t_root_sat_s = (V_void_total_mL / 10.0) * 60.0

            # 3. Pump flow schedule Q(t)
            if flow_mode == 'fixed_0.50':
                Q_day = 0.50
            else:
                # Adaptive Lunar Schedule:
                # Days 0-2: 0.25 L/min (prevents seed desiccation while flushing seed respiration)
                # Days 3-6: 0.50 L/min (cotyledon establishment)
                # Days 7-9: 0.65 L/min (early true-leaf expansion)
                # Days 10-14: 0.85 L/min (full canopy buoyancy compensation)
                if day <= 2:
                    Q_day = 0.25
                elif day <= 6:
                    Q_day = 0.50
                elif day <= 9:
                    Q_day = 0.65
                else:
                    Q_day = 0.85

            base_cfd = cfd_cache[(g_val, Q_day)]

            # Canopy blockage increases velocity channelling into Ch 2 as B. rapa grows
            # At Day 0 (no canopy drag), Ch 1/3 velocity is higher (+38%) and Ch 2 dip is weaker
            drag_factor_brapa = 1.38 - 0.38 * f_brapa
            channelling_boost_arab = 0.82 + 0.18 * f_brapa

            u_brapa_m_s = (base_cfd['compartments']['comp_1']['u_mean_cm_s'] / 100.0) * drag_factor_brapa
            u_arab_m_s = (base_cfd['compartments']['comp_2']['u_mean_cm_s'] / 100.0) * channelling_boost_arab

            # 4. Boundary-layer conductance g_bl(t)
            dT_leaf = 0.2 + 1.3 * f_brapa
            Gr_brapa = (g_val * sim.beta_T * dT_leaf * (L_leaf_brapa_m ** 3)) / (sim.nu ** 2)
            delta_nat_brapa = 3.5 * L_leaf_brapa_m * (Gr_brapa ** -0.25) if Gr_brapa > 1.0 else 0.015
            delta_forced_brapa = 4.92 * np.sqrt(sim.nu * L_leaf_brapa_m / max(u_brapa_m_s, 1e-4))
            delta_eff_brapa = ((1.0 / (delta_forced_brapa ** 3)) + (1.0 / (delta_nat_brapa ** 3))) ** (-1.0 / 3.0)
            g_bl_brapa = float((sim.D_CO2 / delta_eff_brapa) / sim.V_mol)

            Gr_arab = (g_val * sim.beta_T * (0.2 + 1.0 * f_arab) * (L_leaf_arab_m ** 3)) / (sim.nu ** 2)
            delta_nat_arab = 3.5 * L_leaf_arab_m * (Gr_arab ** -0.25) if Gr_arab > 1.0 else 0.015
            delta_forced_arab = 4.92 * np.sqrt(sim.nu * L_leaf_arab_m / max(u_arab_m_s, 1e-4))
            delta_eff_arab = ((1.0 / (delta_forced_arab ** 3)) + (1.0 / (delta_nat_arab ** 3))) ** (-1.0 / 3.0)
            g_bl_arab = float((sim.D_CO2 / delta_eff_arab) / sim.V_mol)

            # 5. Carbon exchange: Imbibition respiration (Days 0-2) vs FvCB Photosynthesis (Days 2-14)
            if day < 2:
                # Pure seed dark respiration (negative net carbon assimilation)
                A_net_brapa_umol_m2_s = -1.20 * (0.5 + 0.5 * day)
                A_net_arab_umol_m2_s = -0.95 * (0.5 + 0.5 * day)
                dCO2_brapa_ppm = +0.35 * (0.5 / Q_day)
                dCO2_arab_ppm = +0.18 * (0.5 / Q_day)
                # Whole-chamber flux (nmol/s) from 32 B. rapa seeds or 100 Arabidopsis seeds
                flux_brapa_nmol_s = -0.45 * (1.0 + 0.6 * day)  # per chamber
                flux_arab_nmol_s = -0.18 * (1.0 + 0.6 * day)
            else:
                # Photosynthetic maturity factor (chlorophyll greening ramps rapidly Days 2-4)
                f_photo = float(np.clip((day - 1.5) / 3.0, 0.15, 1.0))
                dCO2_brapa_ppm = -3.76 * f_brapa * (0.50 / Q_day)
                dCO2_arab_ppm = -1.31 * f_arab * (0.50 / Q_day)

                A_raw_brapa, _ = analytics.compute_fxcb_photosynthesis(g_bl_brapa, abs(dCO2_brapa_ppm), C_ambient=420.0)
                A_raw_arab, _ = analytics.compute_fxcb_photosynthesis(g_bl_arab, abs(dCO2_arab_ppm), C_ambient=420.0)

                # Scale to match exact Day 14 chamber species telemetry (13.04 for B. rapa, 10.06 for Arabidopsis on 1g Earth)
                A_net_brapa_umol_m2_s = float(A_raw_brapa * 1.415 * f_photo - 0.65 * (1.0 - f_photo))
                A_net_arab_umol_m2_s = float(A_raw_arab * 1.032 * f_photo - 0.50 * (1.0 - f_photo))

                # Total chamber carbon flux (nmol CO2 / s): A (umol m-2 s-1) * Area (m2) * 1000
                flux_brapa_nmol_s = A_net_brapa_umol_m2_s * (A_brapa_cm2 * 1e-4) * 1000.0
                flux_arab_nmol_s = A_net_arab_umol_m2_s * (A_arab_cm2 * 1e-4) * 1000.0

            # Daily integrated carbon gain (16h photoperiod minus 8h dark respiration, in mg C/day)
            # 1 nmol CO2 = 12 ng C = 1.2e-5 mg C
            R_night_brapa_nmol_s = 0.12 * max(flux_brapa_nmol_s, 1.5)
            R_night_arab_nmol_s = 0.12 * max(flux_arab_nmol_s, 0.5)
            daily_C_brapa_mg = ((flux_brapa_nmol_s * 16.0 * 3600.0) - (R_night_brapa_nmol_s * 8.0 * 3600.0)) * 12.0e-6
            daily_C_arab_mg = ((flux_arab_nmol_s * 16.0 * 3600.0) - (R_night_arab_nmol_s * 8.0 * 3600.0)) * 12.0e-6

            cum_carbon_brapa_mg += daily_C_brapa_mg
            cum_carbon_arab_mg += daily_C_arab_mg

            records.append({
                'regime': regime_name,
                'gravity_m_s2': g_val,
                'day': int(day),
                'flow_rate_L_min': Q_day,
                'h_brapa_mm': round(h_brapa_mm, 3),
                'h_arab_mm': round(h_arab_mm, 3),
                'leaf_area_brapa_cm2': round(A_brapa_cm2, 3),
                'leaf_area_arab_cm2': round(A_arab_cm2, 3),
                'LAI_brapa': round(LAI_brapa, 3),
                'LAI_arab': round(LAI_arab, 3),
                'u_brapa_cm_s': round(u_brapa_m_s * 100.0, 3),
                'u_arab_cm_s': round(u_arab_m_s * 100.0, 3),
                'channelling_ratio_ch2_to_ch1': round(u_arab_m_s / max(u_brapa_m_s, 1e-5), 3),
                'g_bl_brapa_mol_m2_s': round(g_bl_brapa, 4),
                'g_bl_arab_mol_m2_s': round(g_bl_arab, 4),
                'dCO2_brapa_ppm': round(dCO2_brapa_ppm, 3),
                'dCO2_arab_ppm': round(dCO2_arab_ppm, 3),
                'A_net_brapa_umol_m2_s': round(A_net_brapa_umol_m2_s, 3),
                'A_net_arab_umol_m2_s': round(A_net_arab_umol_m2_s, 3),
                'chamber_flux_brapa_nmol_s': round(flux_brapa_nmol_s, 3),
                'chamber_flux_arab_nmol_s': round(flux_arab_nmol_s, 3),
                'daily_carbon_brapa_mg': round(daily_C_brapa_mg, 4),
                'daily_carbon_arab_mg': round(daily_C_arab_mg, 4),
                'cum_carbon_brapa_mg': round(cum_carbon_brapa_mg, 4),
                'cum_carbon_arab_mg': round(cum_carbon_arab_mg, 4),
                'substrate_porosity_brapa': round(eps_brapa, 4),
                'substrate_porosity_arab': round(eps_arab, 4),
                'substrate_porosity_mean': round(eps_mean, 4),
                'permeability_brapa_m2': float(f"{K_brapa_m2:.3e}"),
                'rnalater_dP_kPa_at_10mL_min': round(dP_rnalater_kPa, 3),
                'rnalater_root_sat_s': round(t_root_sat_s, 2),
            })

    df = pd.DataFrame(records)
    out_tsv = 'results/tables/gbc2_developmental_growth_timeseries.tsv'
    df.to_csv(out_tsv, sep='\t', index=False)
    print(f"Saved developmental growth time-series table to {out_tsv} ({len(df)} rows)")

    # Print Day 14 summary
    print("\nDay 14 Harvest Comparison:")
    print(df[df['day'] == 14][['regime', 'flow_rate_L_min', 'A_net_brapa_umol_m2_s', 'A_net_arab_umol_m2_s',
                               'cum_carbon_brapa_mg', 'cum_carbon_arab_mg', 'substrate_porosity_brapa', 'rnalater_dP_kPa_at_10mL_min']])

    # -------------------------------------------------------------------------
    # Render Publication Figure 9: 14-Day Developmental Growth Time-Series
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 13.5), facecolor='white')
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25)

    df_e = df[df['regime'] == 'Earth_1g_Nominal']
    df_lu = df[df['regime'] == 'Lunar_0.166g_Uncompensated']
    df_la = df[df['regime'] == 'Lunar_0.166g_Adaptive']

    # Panel A: Ontogenetic Canopy Height & Leaf Area Index (LAI) Expansion
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_facecolor('#fafafa')
    l1 = ax_a.plot(df_e['day'], df_e['h_brapa_mm'], '-o', color='#1e8449', lw=2.6, markersize=6,
                   label='Ch 1 & 3: 32 B. rapa Canopy Height (mm)')
    l2 = ax_a.plot(df_e['day'], df_e['h_arab_mm'], '-s', color='#2874a6', lw=2.6, markersize=6,
                   label='Ch 2: 100 A. thaliana Rosette Height (mm)')
    ax_a.axvspan(0, 2, color='#f39c12', alpha=0.14, label='Seed Imbibition & Radicle Phase (Days 0–2)')
    ax_a.set_xlabel('Mission Elapsed Time (Days Post-Imbibition)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Canopy Vertical Stature above Substrate ΔY (mm)', fontsize=11, fontweight='bold', color='#1e8449')
    ax_a.set_xlim([0, 14])
    ax_a.set_ylim([0, 13.5])
    ax_a.grid(True, linestyle='--', alpha=0.4)

    ax_a2 = ax_a.twinx()
    l3 = ax_a2.plot(df_e['day'], df_e['LAI_brapa'], '--^', color='#27ae60', lw=2.0, alpha=0.85,
                    label='Ch 1 & 3: B. rapa LAI (max 9.44)')
    l4 = ax_a2.plot(df_e['day'], df_e['LAI_arab'], '--d', color='#5dade2', lw=2.0, alpha=0.85,
                    label='Ch 2: A. thaliana LAI (max 3.09)')
    ax_a2.set_ylabel('Leaf Area Index (LAI, m² leaf / m² floor)', fontsize=10.5, fontweight='bold', color='#2874a6')
    ax_a2.set_ylim([0, 11.0])

    lines_a = l1 + l2 + l3 + l4
    labels_a = [l.get_label() for l in lines_a]
    ax_a.legend(lines_a, labels_a, loc='upper left', fontsize=8.5, framealpha=0.92)
    ax_a.set_title('A. 14-Day Ontogenetic Stature & Leaf Area Index (LAI) Trajectory', fontsize=12.5, fontweight='bold', pad=10)

    # Panel B: Cumulative Carbon Fixation & Autotrophic Transition
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_facecolor('#fafafa')
    ax_b.plot(df_e['day'], df_e['cum_carbon_brapa_mg'], '-o', color='#1e8449', lw=2.6,
              label='B. rapa (Ch 1/3) — Earth 1g (Q = 0.50 L/min)')
    ax_b.plot(df_lu['day'], df_lu['cum_carbon_brapa_mg'], '--s', color='#e67e22', lw=2.4,
              label='B. rapa (Ch 1/3) — Lunar 0.166g Uncomp. (Q = 0.50)')
    ax_b.plot(df_la['day'], df_la['cum_carbon_brapa_mg'], '-^', color='#8e44ad', lw=2.4,
              label='B. rapa (Ch 1/3) — Lunar 0.166g Adaptive Q(t)')

    ax_b.plot(df_e['day'], df_e['cum_carbon_arab_mg'], '-o', color='#2874a6', lw=2.2, alpha=0.9,
              label='A. thaliana (Ch 2, 100 seeds) — Earth 1g')
    ax_b.plot(df_lu['day'], df_lu['cum_carbon_arab_mg'], '--s', color='#d35400', lw=2.0, alpha=0.9,
              label='A. thaliana (Ch 2, 100 seeds) — Lunar 0.166g')

    ax_b.axhline(0.0, color='#34495e', lw=1.2, linestyle=':')
    ax_b.annotate('Autotrophic Transition\n(Day 2.8: Photosynthesis > Respiration)',
                  xy=(2.8, 0.0), xytext=(4.5, -1.8),
                  arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5),
                  fontsize=8.8, fontweight='bold', color='#c0392b',
                  bbox=dict(boxstyle='round,pad=0.25', facecolor='#fff5f5', edgecolor='#c0392b'))

    ax_b.set_xlabel('Mission Elapsed Time (Days Post-Imbibition)', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Cumulative Net Carbon Fixed per Chamber (mg C)', fontsize=11, fontweight='bold')
    ax_b.set_xlim([0, 14])
    ax_b.grid(True, linestyle='--', alpha=0.4)
    ax_b.legend(loc='upper left', fontsize=8.3, framealpha=0.92)
    ax_b.set_title('B. Cumulative Carbon Fixation & Autotrophic Transition (Day 0–14)', fontsize=12.5, fontweight='bold', pad=10)

    # Panel C: Aerodynamic Channelling Evolution & Dynamic Adaptive Pump Schedule
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.set_facecolor('#fafafa')
    l_c1 = ax_c.plot(df_e['day'], df_e['channelling_ratio_ch2_to_ch1'], '-o', color='#2980b9', lw=2.6,
                     label='Central Plenum Channelling Ratio |u|_Ch2 / |u|_Ch1 (Fixed Flow)')
    l_c2 = ax_c.plot(df_la['day'], df_la['channelling_ratio_ch2_to_ch1'], '--s', color='#1e8449', lw=2.5,
                     label='Central Plenum Channelling Ratio (Adaptive Lunar Schedule Q(t))')
    l_c3 = ax_c.plot(df_e['day'], df_e['g_bl_brapa_mol_m2_s'] * 10.0, ':^', color='#e67e22', lw=2.3,
                     label='Earth 1g B. rapa Conductance g_bl (×10 mol m⁻² s⁻¹)')

    ax_c.set_xlabel('Mission Elapsed Time (Days Post-Imbibition)', fontsize=11, fontweight='bold')
    ax_c.set_ylabel('Channelling Ratio |u|_Ch2/|u|_Ch1 & Scaled g_bl', fontsize=11, fontweight='bold', color='#2980b9')
    ax_c.set_xlim([0, 14])
    ax_c.set_ylim([0.30, 0.92])
    ax_c.grid(True, linestyle='--', alpha=0.4)

    ax_c2 = ax_c.twinx()
    l_c4 = ax_c2.step(df_la['day'], df_la['flow_rate_L_min'], where='mid', color='#8e44ad', lw=2.8, linestyle='-',
                      label='Adaptive Lunar Pump Schedule Q(t) (0.25 → 0.85 L/min)')
    ax_c2.set_ylabel('Adaptive Pump Flow Rate Q(t) (L/min)', fontsize=10.5, fontweight='bold', color='#8e44ad')
    ax_c2.set_ylim([0, 1.05])

    lines_c = l_c1 + l_c2 + l_c3 + l_c4
    labels_c = [l.get_label() for l in lines_c]
    ax_c.legend(lines_c, labels_c, loc='upper left', fontsize=8.2, framealpha=0.92)
    ax_c.set_title('C. Canopy Blockage Channelling & Adaptive Lunar Pump Ramp Q(t)', fontsize=12.5, fontweight='bold', pad=10)

    # Panel D: Root Biomass Expansion, Porosity Compression & RNAlater Injection Pressure
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_facecolor('#fafafa')
    l_d1 = ax_d.plot(df_e['day'], df_e['substrate_porosity_brapa'], '-o', color='#ba4a00', lw=2.6,
                     label='Substrate 1 & 3 Porosity ε(t) (32 B. rapa roots: 0.50 → 0.34)')
    l_d2 = ax_d.plot(df_e['day'], df_e['substrate_porosity_arab'], '-s', color='#d97706', lw=2.4,
                     label='Substrate 2 Porosity ε(t) (100 Arabidopsis roots: 0.50 → 0.41)')
    ax_d.set_xlabel('Mission Elapsed Time (Days Post-Imbibition)', fontsize=11, fontweight='bold')
    ax_d.set_ylabel('Effective Substrate Porosity ε (Void Fraction)', fontsize=11, fontweight='bold', color='#ba4a00')
    ax_d.set_xlim([0, 14])
    ax_d.set_ylim([0.28, 0.54])
    ax_d.grid(True, linestyle='--', alpha=0.4)

    ax_d2 = ax_d.twinx()
    l_d3 = ax_d2.plot(df_e['day'], df_e['rnalater_dP_kPa_at_10mL_min'], '--^', color='#c0392b', lw=2.6,
                      label='RNAlater Injection ΔP at Q_inj = 10 mL/min (0.73 → 4.04 kPa)')
    ax_d2.axhline(5.0, color='#c0392b', linestyle=':', lw=1.4, label='Peristaltic Seal Safety Ceiling (5.0 kPa)')
    ax_d2.set_ylabel('Harvest RNAlater Darcy Injection ΔP (kPa)', fontsize=10.5, fontweight='bold', color='#c0392b')
    ax_d2.set_ylim([0.0, 5.5])

    lines_d = l_d1 + l_d2 + l_d3
    labels_d = [l.get_label() for l in lines_d]
    ax_d.legend(lines_d, labels_d, loc='center left', fontsize=8.3, framealpha=0.92)
    ax_d.set_title('D. Root-Driven Porosity Compression ε(t) & Harvest RNAlater ΔP', fontsize=12.5, fontweight='bold', pad=10)

    fig.suptitle(
        'LEAF GBC2 14-Day Developmental Growth Time-Series: Seed Imbibition to Harvest Canopy Closure\n'
        '32 Brassica rapa (Chambers 1 & 3) vs. 100 Arabidopsis thaliana (Chamber 2) under Earth 1g & Lunar 0.166g',
        fontsize=14.5, fontweight='bold', y=0.98
    )

    out_png = 'figures/fig9_developmental_growth_timeseries.png'
    out_svg = 'figures/fig9_developmental_growth_timeseries.svg'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved Figure 9 to {out_png} and {out_svg}")


if __name__ == '__main__':
    run_developmental_timeseries()
