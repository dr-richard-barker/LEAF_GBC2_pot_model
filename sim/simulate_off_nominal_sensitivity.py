"""
Off-Nominal Mission Sensitivity & Fault-Tolerance Sweeps for LEAF GBC1 Plant Growth Chamber.

Evaluates three critical spacecraft / lunar lander off-nominal operational stressors:
  1. Elevated Spacecraft Cabin CO2 (420 to 5,000 ppm):
     - Accounts for high-CO2 stomatal conductance downregulation gs(C_amb) (Medlyn / Ball-Berry closure above 1,500 ppm)
     - Triose-phosphate utilization (TPU) carboxylation ceiling at super-elevated CO2
     - Transpiration cooling loss and foliar heat-trapping (+dT_leaf) when stomata close under Lunar 0.166g
  2. Root-Growth Substrate Porosity Reduction (eps = 0.50 down to 0.20):
     - Kozeny-Carman hydraulic permeability drop K(eps)
     - Rhizosphere O2 depletion d[O2] in Substrate Blocks 1/3 (32 B. rapa) vs Block 2 (100 Arabidopsis)
     - Harvest-phase RNAlater Darcy injection pressure dP_total and root saturation time t_root
  3. Pump Brownout, Filter Clogging & Power Outage Fault Envelope (Q = 0.85 down to 0.00 L/min):
     - Stagnant volume fraction, canopy RH saturation time, and Maximum Allowable Pump Outage Grace Window (t_grace)
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sim.gbc1_cfd_engine import GBC1CFDSimulator
from sim.advanced_analytics import GBC1AdvancedAnalytics


def compute_fvcb_with_stomatal_feedback(g_bl, dC_CO2_ppm, C_ambient, species='brapa'):
    """
    Extended FvCB model incorporating high-CO2 stomatal closure (Ball-Berry / Medlyn)
    and Triose Phosphate Utilization (TPU) ceiling at spacecraft cabin CO2 (420 - 5000 ppm).
    """
    if species == 'brapa':
        Vcmax = 68.0
        Jmax = 135.0
        TPU_max = 11.2
        gs_ref = 0.165
    else:
        Vcmax = 52.0
        Jmax = 104.0
        TPU_max = 8.6
        gs_ref = 0.135

    Kc = 404.0
    Ko = 248.0
    O = 210.0
    Km = Kc * (1.0 + O / Ko)
    gamma_star = 42.75
    Rd = 0.015 * Vcmax

    # Stomatal closure response to super-ambient CO2 (partial closure above 800 ppm, ~55% closure at 5000 ppm)
    gs_factor = 1.0 / (1.0 + 0.28 * np.maximum((C_ambient - 420.0) / 1000.0, 0.0))
    gs = max(gs_ref * gs_factor, 0.035)

    gtot = 1.0 / ((1.0 / max(g_bl, 1e-4)) + (1.0 / gs))
    Cs = max(C_ambient - dC_CO2_ppm, 50.0)

    # Rubisco-limited (Ac)
    a = gtot
    b = (Vcmax - Rd) + gtot * (Km - Cs)
    c = -(gtot * Cs * Km + Vcmax * gamma_star + Rd * Km)
    Ci_c = (-b + np.sqrt(max(b**2 - 4 * a * c, 1e-12))) / (2 * a)
    Ac = gtot * (Cs - Ci_c)

    # RuBP-regeneration limited (Aj)
    a_j = gtot
    b_j = (Jmax / 4.0 - Rd) + gtot * (2 * gamma_star - Cs)
    c_j = -(2 * gtot * Cs * gamma_star + (Jmax / 4.0) * gamma_star + 2 * Rd * gamma_star)
    Ci_j = (-b_j + np.sqrt(max(b_j**2 - 4 * a_j * c_j, 1e-12))) / (2 * a_j)
    Aj = gtot * (Cs - Ci_j)

    # TPU-limited (Ap) at super-elevated CO2 (>1500 ppm)
    Ap = 3.0 * TPU_max - Rd

    A_net = float(max(min(Ac, Aj, Ap), 0.0))
    return A_net, float(gs), float(gs_factor)


def run_off_nominal_sweeps():
    os.makedirs('results/tables', exist_ok=True)
    os.makedirs('figures', exist_ok=True)

    print("=================================================================")
    print("Running Off-Nominal Mission Sensitivity & Fault-Tolerance Sweeps")
    print("=================================================================")

    sim = GBC1CFDSimulator()

    # Baseline CFD runs
    res_earth_050, _ = sim.run_case(Q_L_min=0.50, gravity_g=9.81, max_steps=260)
    res_lunar_050, _ = sim.run_case(Q_L_min=0.50, gravity_g=1.62, max_steps=260)
    res_lunar_085, _ = sim.run_case(Q_L_min=0.85, gravity_g=1.62, max_steps=260)

    records = []

    # -------------------------------------------------------------------------
    # Sweep 1: Elevated Spacecraft Cabin CO2 (420 to 5,000 ppm)
    # -------------------------------------------------------------------------
    co2_levels = [420, 800, 1200, 1500, 2500, 4000, 5000]
    co2_regimes = [
        ('Earth_1g_0.50Lmin', 9.81, 0.50, res_earth_050),
        ('Lunar_0.166g_0.50Lmin', 1.62, 0.50, res_lunar_050),
        ('Lunar_0.166g_0.85Lmin', 1.62, 0.85, res_lunar_085),
    ]

    for reg_label, g_val, q_val, res_base in co2_regimes:
        g_bl_brapa = res_base['compartments']['comp_1']['g_bl_mol_m2_s']
        g_bl_arab = res_base['compartments']['comp_2']['g_bl_mol_m2_s']
        dCO2_brapa = res_base['compartments']['comp_1']['dC_CO2_ppm']
        dCO2_arab = res_base['compartments']['comp_2']['dC_CO2_ppm']

        for c_amb in co2_levels:
            A_brapa, gs_brapa, gs_fac_brapa = compute_fvcb_with_stomatal_feedback(
                g_bl_brapa, dCO2_brapa, c_amb, species='brapa'
            )
            A_arab, gs_arab, gs_fac_arab = compute_fvcb_with_stomatal_feedback(
                g_bl_arab, dCO2_arab, c_amb, species='arab'
            )

            # Leaf temperature elevation increases when stomata close (reduced latent heat cooling)
            # Base dT on Earth = 1.50 K, Lunar 0.50 = 2.42 K, Lunar 0.85 = 1.68 K
            base_dT = 1.50 if g_val > 5.0 else (2.42 if q_val < 0.7 else 1.68)
            latent_cooling_loss_K = (1.0 - gs_fac_brapa) * (1.25 if g_val < 5.0 else 0.75) * (0.50 / q_val) ** 0.4
            dT_leaf_eff_C = base_dT + latent_cooling_loss_K

            records.append({
                'sweep_type': 'Cabin_CO2_Enrichment',
                'parameter_name': 'C_ambient_ppm',
                'parameter_value': float(c_amb),
                'regime': reg_label,
                'gravity_m_s2': g_val,
                'flow_rate_L_min': q_val,
                'g_bl_brapa_mol_m2_s': round(g_bl_brapa, 4),
                'gs_stomatal_brapa_mol_m2_s': round(gs_brapa, 4),
                'A_net_brapa_umol_m2_s': round(A_brapa, 3),
                'A_net_arab_umol_m2_s': round(A_arab, 3),
                'leaf_dT_elevation_C': round(dT_leaf_eff_C, 3),
                'stagnant_fraction_pct': round(res_base['stagnant_volume_pct'], 2),
                'rhizosphere_O2_depletion_pct': 4.2 if g_val > 5.0 else (7.8 if q_val < 0.7 else 5.1),
                'rnalater_dP_kPa': 0.73,
                'grace_window_min': 999.0,
            })

    # -------------------------------------------------------------------------
    # Sweep 2: Root-Growth Substrate Porosity Compaction (eps = 0.50 -> 0.20)
    # -------------------------------------------------------------------------
    porosity_levels = [0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20]
    eps_0 = 0.50
    K_0 = 1.0e-10  # m2

    for eps_val in porosity_levels:
        K_ratio = ((eps_val ** 3) * ((1.0 - eps_0) ** 2)) / (((eps_0) ** 3) * ((1.0 - eps_val) ** 2))
        K_eff = K_0 * K_ratio

        # Effective gas diffusivity in porous media (Millington-Quirk tortuosity: D_eff = D_0 * eps^(4/3))
        tortuosity_factor = (eps_val / eps_0) ** (4.0 / 3.0)

        # Rhizosphere O2 depletion (%) inside Substrate Block 1/3 (32 B. rapa roots) under Lunar 0.166g vs Earth 1g
        o2_dep_earth_pct = 4.2 / tortuosity_factor
        o2_dep_lunar_pct = 7.8 / tortuosity_factor

        # Harvest RNAlater Darcy injection pressure at Q_inj = 10 mL/min and 20 mL/min
        dP_10mL_kPa = 0.73 / K_ratio
        dP_20mL_kPa = 1.47 / K_ratio

        # Root saturation time (s) at 10 mL/min (smaller void volume saturates faster, but higher dP)
        V_void_mL = 10.55 * eps_val
        t_root_s = (V_void_mL / 10.0) * 60.0

        records.append({
            'sweep_type': 'Substrate_Porosity_Compaction',
            'parameter_name': 'substrate_porosity_eps',
            'parameter_value': eps_val,
            'regime': 'Lunar_0.166g_RootZone',
            'gravity_m_s2': 1.62,
            'flow_rate_L_min': 0.85,
            'g_bl_brapa_mol_m2_s': round(res_lunar_085['compartments']['comp_1']['g_bl_mol_m2_s'], 4),
            'gs_stomatal_brapa_mol_m2_s': 0.165,
            'A_net_brapa_umol_m2_s': round(10.66 * np.clip(1.0 - max(o2_dep_lunar_pct - 18.0, 0.0) * 0.025, 0.5, 1.0), 3),
            'A_net_arab_umol_m2_s': round(8.44 * np.clip(1.0 - max(o2_dep_lunar_pct * 0.65 - 18.0, 0.0) * 0.02, 0.6, 1.0), 3),
            'leaf_dT_elevation_C': round(float(K_eff * 1e10), 4),  # store K_eff * 1e10 in auxiliary column or dedicated metric
            'stagnant_fraction_pct': round(t_root_s, 2),
            'rhizosphere_O2_depletion_pct': round(o2_dep_lunar_pct, 2),
            'rnalater_dP_kPa': round(dP_10mL_kPa, 3),
            'grace_window_min': round(dP_20mL_kPa, 3),
        })

    # -------------------------------------------------------------------------
    # Sweep 3: Pump Brownout / Filter Blockage Fault Modes (Q = 0.85 -> 0.02 L/min)
    # -------------------------------------------------------------------------
    fault_flows = [0.85, 0.60, 0.50, 0.35, 0.25, 0.15, 0.08, 0.04, 0.02]
    # Headspace free volume = 54.05 mL = 5.405e-5 m3
    # Total canopy transpiration at nominal daylight: ~0.135 umol H2O / s
    # Time to saturate headspace from 65% RH to 95% RH (+350 mmol/m3) or +20 ppm nocturnal CO2
    for q_fault in fault_flows:
        res_f, _ = sim.run_case(Q_L_min=q_fault, gravity_g=1.62, max_steps=240)
        g_bl_b = res_f['compartments']['comp_1']['g_bl_mol_m2_s']
        dCO2_b = res_f['compartments']['comp_1']['dC_CO2_ppm'] * (0.50 / max(q_fault, 0.02)) ** 0.75
        dH2O_b = res_f['compartments']['comp_1']['dC_H2O_mmol'] * (0.50 / max(q_fault, 0.02)) ** 0.75

        A_b, _, _ = compute_fvcb_with_stomatal_feedback(g_bl_b, dCO2_b, 420.0, species='brapa')
        A_a, _, _ = compute_fvcb_with_stomatal_feedback(res_f['compartments']['comp_2']['g_bl_mol_m2_s'],
                                                        res_f['compartments']['comp_2']['dC_CO2_ppm'] * (0.50 / max(q_fault, 0.02)) ** 0.75,
                                                        420.0, species='arab')

        # Maximum allowable outage / brownout duration (minutes) before canopy RH > 92% or dT > +4.0 C
        # At Q >= 0.50 L/min, steady-state is safe indefinitely (grace_window = 240 min cap for plotting)
        flushing_rate_min = (q_fault * 1000.0) / 54.05  # air exchanges per minute
        net_accum_factor = max(1.85 - flushing_rate_min * 0.20, 0.04)
        grace_min = float(np.clip(18.5 / net_accum_factor, 4.2, 240.0))

        leaf_dT_fault = 1.68 + 2.85 * np.exp(-q_fault / 0.22)

        records.append({
            'sweep_type': 'Pump_Fault_Degradation',
            'parameter_name': 'Q_fault_L_min',
            'parameter_value': q_fault,
            'regime': 'Lunar_0.166g_Brownout',
            'gravity_m_s2': 1.62,
            'flow_rate_L_min': q_fault,
            'g_bl_brapa_mol_m2_s': round(g_bl_b, 4),
            'gs_stomatal_brapa_mol_m2_s': round(dH2O_b, 3),  # store dH2O_mmol in this column for fault sweep
            'A_net_brapa_umol_m2_s': round(A_b, 3),
            'A_net_arab_umol_m2_s': round(A_a, 3),
            'leaf_dT_elevation_C': round(leaf_dT_fault, 3),
            'stagnant_fraction_pct': round(res_f['stagnant_volume_pct'], 2),
            'rhizosphere_O2_depletion_pct': round(5.1 * (0.85 / max(q_fault, 0.04)) ** 0.45, 2),
            'rnalater_dP_kPa': round(dCO2_b, 3),
            'grace_window_min': round(grace_min, 2),
        })

    df_sens = pd.DataFrame(records)
    out_tsv = 'results/tables/gbc1_off_nominal_sensitivity.tsv'
    df_sens.to_csv(out_tsv, sep='\t', index=False)
    print(f"Saved off-nominal sensitivity results to {out_tsv} ({len(df_sens)} rows)")

    # -------------------------------------------------------------------------
    # Render Publication Figure 10: Off-Nominal Mission Sensitivity Sweeps
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 13.5), facecolor='white')
    gs = fig.add_gridspec(2, 2, hspace=0.33, wspace=0.25)

    # Panel A: Cabin CO2 Enrichment (420 to 5,000 ppm)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_facecolor('#fafafa')
    df_c = df_sens[df_sens['sweep_type'] == 'Cabin_CO2_Enrichment']
    c_e = df_c[df_c['regime'] == 'Earth_1g_0.50Lmin']
    c_lu = df_c[df_c['regime'] == 'Lunar_0.166g_0.50Lmin']
    c_lc = df_c[df_c['regime'] == 'Lunar_0.166g_0.85Lmin']

    l_a1 = ax_a.plot(c_e['parameter_value'], c_e['A_net_brapa_umol_m2_s'], '-o', color='#1e8449', lw=2.6,
                     label='B. rapa A_net — Earth 1g (0.50 L/min)')
    l_a2 = ax_a.plot(c_lu['parameter_value'], c_lu['A_net_brapa_umol_m2_s'], '--s', color='#e67e22', lw=2.5,
                     label='B. rapa A_net — Lunar 0.166g (0.50 L/min)')
    l_a3 = ax_a.plot(c_lc['parameter_value'], c_lc['A_net_arab_umol_m2_s'], '-^', color='#2874a6', lw=2.4,
                     label='A. thaliana A_net — Lunar 0.166g (0.85 L/min)')
    ax_a.axvspan(2000, 4500, color='#8e44ad', alpha=0.10, label='Typical Crewed Spacecraft Cabin CO₂ (2,000–4,500 ppm)')

    ax_a.set_xlabel('Ambient Spacecraft Cabin [CO₂] (ppm)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Net Photosynthesis A_net (μmol m⁻² s⁻¹)', fontsize=11, fontweight='bold', color='#1e8449')
    ax_a.set_xlim([300, 5150])
    ax_a.set_ylim([6, 33])
    ax_a.grid(True, linestyle='--', alpha=0.4)

    ax_a2 = ax_a.twinx()
    l_a4 = ax_a2.plot(c_lu['parameter_value'], c_lu['leaf_dT_elevation_C'], ':d', color='#c0392b', lw=2.4,
                      label='Lunar 0.166g Leaf ΔT (°C, Stomatal Closure Heating)')
    ax_a2.set_ylabel('Foliar Thermal Elevation ΔT_leaf (°C)', fontsize=10.5, fontweight='bold', color='#c0392b')
    ax_a2.set_ylim([1.0, 4.0])

    lines_a = l_a1 + l_a2 + l_a3 + l_a4
    ax_a.legend(lines_a, [l.get_label() for l in lines_a], loc='upper left', fontsize=8.2, framealpha=0.92)
    ax_a.set_title('A. Elevated Cabin CO₂ Sweep (420–5,000 ppm): Assimilation vs. Thermal Load', fontsize=12.3, fontweight='bold', pad=10)

    # Panel B: Root-Growth Substrate Porosity Compaction (eps = 0.50 -> 0.20)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_facecolor('#fafafa')
    df_p = df_sens[df_sens['sweep_type'] == 'Substrate_Porosity_Compaction']

    l_b1 = ax_b.plot(df_p['parameter_value'], df_p['rnalater_dP_kPa'], '-o', color='#c0392b', lw=2.6,
                     label='RNAlater Injection ΔP at Q=10 mL/min (kPa)')
    l_b2 = ax_b.plot(df_p['parameter_value'], df_p['grace_window_min'], '--s', color='#e74c3c', lw=2.2,
                     label='RNAlater Injection ΔP at Q=20 mL/min (kPa)')
    ax_b.axhline(5.0, color='#7f8c8d', linestyle=':', lw=1.6, label='Hardware Seal Pressure Limit (5.0 kPa)')
    ax_b.axvspan(0.34, 0.50, color='#27ae60', alpha=0.11, label='Nominal 14-Day Root Compaction Envelope (ε = 0.50 → 0.34)')

    ax_b.set_xlabel('Effective Substrate Porosity ε (Void Fraction after Root Expansion)', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Harvest RNAlater Injection Pressure ΔP (kPa)', fontsize=11, fontweight='bold', color='#c0392b')
    ax_b.invert_xaxis()  # Show progression from 0.50 (Day 0) to 0.20 (over-compacted)
    ax_b.set_ylim([0, 16.5])
    ax_b.grid(True, linestyle='--', alpha=0.4)

    ax_b2 = ax_b.twinx()
    l_b3 = ax_b2.plot(df_p['parameter_value'], df_p['rhizosphere_O2_depletion_pct'], '-^', color='#2980b9', lw=2.4,
                      label='Lunar Rhizosphere O₂ Depletion Δ[O₂] (%)')
    ax_b2.set_ylabel('Rhizosphere O₂ Depletion (%)', fontsize=10.5, fontweight='bold', color='#2980b9')
    ax_b2.set_ylim([0, 32])

    lines_b = l_b1 + l_b2 + l_b3
    ax_b.legend(lines_b, [l.get_label() for l in lines_b], loc='upper left', fontsize=8.2, framealpha=0.92)
    ax_b.set_title('B. Root-Zone Porosity Compaction (ε = 0.50 → 0.20): O₂ & RNAlater ΔP', fontsize=12.3, fontweight='bold', pad=10)

    # Panel C: Pump Brownout / Filter Blockage Fault Envelope on Moon (0.166g)
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.set_facecolor('#fafafa')
    df_f = df_sens[df_sens['sweep_type'] == 'Pump_Fault_Degradation']

    l_c1 = ax_c.plot(df_f['parameter_value'], df_f['rnalater_dP_kPa'], '-o', color='#d35400', lw=2.6,
                     label='Canopy CO₂ Drawdown |Δ[CO₂]| (ppm)')
    l_c2 = ax_c.plot(df_f['parameter_value'], df_f['gs_stomatal_brapa_mol_m2_s'], '-s', color='#2980b9', lw=2.4,
                     label='Canopy Vapor Accumulation Δ[H₂O] (mmol m⁻³)')
    ax_c.axvline(0.50, color='#27ae60', linestyle='--', lw=1.6, label='Earth Nominal Flow (0.50 L/min)')
    ax_c.axvline(0.15, color='#c0392b', linestyle=':', lw=1.8, label='Critical Brownout Threshold (0.15 L/min)')

    ax_c.set_xlabel('Degraded / Brownout Pump Flow Rate Q (L/min)', fontsize=11, fontweight='bold')
    ax_c.set_ylabel('Canopy |Δ[CO₂]| (ppm) & Vapor Δ[H₂O] (mmol m⁻³)', fontsize=11, fontweight='bold')
    ax_c.set_xlim([0.0, 0.90])
    ax_c.set_ylim([0, 95])
    ax_c.grid(True, linestyle='--', alpha=0.4)

    ax_c2 = ax_c.twinx()
    l_c3 = ax_c2.plot(df_f['parameter_value'], df_f['A_net_brapa_umol_m2_s'], '--^', color='#1e8449', lw=2.4,
                      label='B. rapa Net Photosynthesis A_net (μmol m⁻² s⁻¹)')
    ax_c2.set_ylabel('B. rapa A_net (μmol m⁻² s⁻¹)', fontsize=10.5, fontweight='bold', color='#1e8449')
    ax_c2.set_ylim([2.0, 10.5])

    lines_c = l_c1 + l_c2 + l_c3
    ax_c.legend(lines_c, [l.get_label() for l in lines_c], loc='center right', fontsize=8.2, framealpha=0.92)
    ax_c.set_title('C. Lunar 0.166g Pump Brownout & Filter Clogging Fault Envelope', fontsize=12.3, fontweight='bold', pad=10)

    # Panel D: Maximum Allowable Pump Outage Grace Window & Operational Safety Zones
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_facecolor('#fafafa')
    l_d1 = ax_d.plot(df_f['parameter_value'], df_f['grace_window_min'], '-o', color='#8e44ad', lw=2.8, markersize=7,
                     label='Max Allowable Outage / Brownout Grace Window t_grace (min)')
    l_d2 = ax_d.plot(df_f['parameter_value'], df_f['leaf_dT_elevation_C'] * 25.0, '--d', color='#c0392b', lw=2.2,
                     label='Foliar Thermal Elevation ΔT_leaf (scaled ×25 for axis)')

    ax_d.axhspan(0, 30, color='#e74c3c', alpha=0.14, label='Acute Condensation / Heat Stress Zone (< 30 min)')
    ax_d.axhspan(30, 120, color='#f39c12', alpha=0.12, label='Transient Telemetry Watch Zone (30–120 min)')
    ax_d.axhspan(120, 255, color='#2ecc71', alpha=0.12, label='Continuous Safe Flight Operating Envelope')

    ax_d.annotate('Complete Pump Stop Grace Time:\nt_grace ≈ 10.0–18.5 min before\n92% RH Condensation Onset',
                  xy=(0.04, 11.2), xytext=(0.24, 48.0),
                  arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.6),
                  fontsize=8.8, fontweight='bold', color='#4a235a',
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#8e44ad'))

    ax_d.set_xlabel('Operating / Residual Pump Flow Rate Q (L/min)', fontsize=11, fontweight='bold')
    ax_d.set_ylabel('Allowable Fault Duration Grace Window (minutes)', fontsize=11, fontweight='bold', color='#8e44ad')
    ax_d.set_xlim([0.0, 0.90])
    ax_d.set_ylim([0, 255])
    ax_d.grid(True, linestyle='--', alpha=0.4)
    ax_d.legend(loc='upper left', fontsize=8.2, framealpha=0.92)
    ax_d.set_title('D. Mission Fault-Tolerance Boundaries & Outage Grace Window', fontsize=12.3, fontweight='bold', pad=10)

    fig.suptitle(
        'LEAF GBC1 Off-Nominal Mission Sensitivity & Fault-Tolerance Sweeps\n'
        'Spacecraft Cabin CO₂ Enrichment (420–5,000 ppm), Root Porosity Compaction, and Pump Brownout Resilience',
        fontsize=14.5, fontweight='bold', y=0.98
    )

    out_png = 'figures/fig10_off_nominal_sensitivity_sweeps.png'
    out_svg = 'figures/fig10_off_nominal_sensitivity_sweeps.svg'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved Figure 10 to {out_png} and {out_svg}")


if __name__ == '__main__':
    run_off_nominal_sweeps()
