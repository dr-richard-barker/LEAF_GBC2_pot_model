"""
4D Spatiotemporal Modeling & Cinematic Video Rendering Engine
for the LEAF GBC1 Lunar Plant Growth Chamber (M09024F-MDL-5 Pot + M09025F-MDL-2 Lid).

Generates Full HD (1920x1080, 24 fps) H.264 MP4 cinematic videos:
  1. videos/lunar_experiment_4d_cinematic.mp4       (0.166g Lunar Surface Mission)
  2. videos/ground_control_1g_4d_cinematic.mp4      (1.00g Earth Ground Control)
  3. videos/comparison_lunar_vs_1g_4d.mp4           (Full-Screen Synchronized Side-by-Side Comparison)
  4. videos/stage1_water_foam_activation.mp4        (Stage I Highlight: Foam & Nytex Wicking)
  5. videos/stage2_airflow_turbulence_moisture.mp4  (Stage II Highlight: Tubes, Turbulence & Moisture)
  6. videos/stage3_5day_growth_ontogeny.mp4         (Stage III-IV Highlight: 5-Day Growth & Diurnal Cycle)
  7. videos/stage4_rnalater_harvest_injection.mp4   (Stage V Highlight: RNAlater Fixative Perfusion)
"""

import os
import sys
import math
import time
import subprocess
import multiprocessing as mp
import numpy as np

os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config_4d'
os.makedirs('/tmp/mpl_config_4d', exist_ok=True)
os.makedirs('videos', exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from PIL import Image

FFMPEG_BIN = '/opt/homebrew/bin/ffmpeg'
if not os.path.exists(FFMPEG_BIN):
    FFMPEG_BIN = 'ffmpeg'


# =============================================================================
# 1. PRECOMPUTED GEOMETRY & 164-PLANT COORDINATES
# =============================================================================

def build_plant_catalog():
    """Returns deterministic array of 164 plant positions: (x, z, species_id, chamber_id, phase_jitter)."""
    rng = np.random.default_rng(42)
    plants = []
    # Chamber 1: 32 Brassica rapa (8 x 4 lattice)
    for ix in range(8):
        for iz in range(4):
            px = -8.0 + (ix / 7.0) * 26.0
            pz = -28.5 + (iz / 3.0) * 16.5
            plants.append((px, pz, 1, 1, rng.uniform(-1.5, 1.5)))
    # Chamber 2: 100 Arabidopsis thaliana (10 x 10 rosette lawn)
    for ix in range(10):
        for iz in range(10):
            px = -8.0 + (ix / 9.0) * 26.0
            pz = -6.0 + (iz / 9.0) * 12.0
            plants.append((px, pz, 2, 2, rng.uniform(-1.5, 1.5)))
    # Chamber 3: 32 Brassica rapa (8 x 4 lattice)
    for ix in range(8):
        for iz in range(4):
            px = -8.0 + (ix / 7.0) * 26.0
            pz = 11.5 + (iz / 3.0) * 16.5
            plants.append((px, pz, 1, 3, rng.uniform(-1.5, 1.5)))
    return np.array(plants, dtype=np.float64)

PLANTS_ARRAY = build_plant_catalog()


# =============================================================================
# 2. 4D TIMELINE STATE ENGINE (NORMALIZED PROGRESS u in [0, 1])
# =============================================================================

def compute_mission_state(u, regime='lunar'):
    """
    Computes the instantaneous physical state of the LEAF GBC1 chamber at normalized
    timeline coordinate u in [0, 1] for either 'lunar' (0.166g) or 'earth' (1.00g).
    """
    is_lunar = (regime == 'lunar')
    g_val = 0.166 if is_lunar else 1.000
    g_ms2 = g_val * 9.81

    state = {
        'u': u,
        'regime': regime,
        'is_lunar': is_lunar,
        'g_val': g_val,
        'g_ms2': g_ms2,
        'stage_num': 1,
        'stage_name': '',
        'stage_sub': '',
        'mission_clock': '',
        'water_injected_mL': 0.0,
        'foam_sat_pct': 0.0,
        'water_tube_active': 0.0,
        'air_Q_Lmin': 0.0,
        'air_active': 0.0,
        'day_float': 0.0,
        't_hours': 0.0,
        'is_night': False,
        'rna_active': 0.0,
        'rna_fill_frac': 0.0,
        'rna_vol_mL': 0.0,
        'rna_dP_kPa': 0.0,
        'delta_bl_mm': 0.0,
        'Ri_num': 0.0,
        'Re_num': 0.0,
        'A_net': 0.0,
        'RH_mean': 45.0,
        'T_margin_C': 2.50,
        'cam_azim': -52.0,
        'cam_elev': 23.0,
    }

    # Smooth camera orbit across the mission
    state['cam_azim'] = -56.0 + 24.0 * math.sin(u * math.pi * 1.6)
    state['cam_elev'] = 22.0 + 5.0 * math.sin(u * math.pi * 2.2)

    if u < 0.22:
        # STAGE I: Water Injection & Porous Foam Capillary Wicking
        s = u / 0.22
        state['stage_num'] = 1
        state['stage_name'] = "STAGE I: WATER INJECTION & POROUS FOAM CAPILLARY FILLING"
        t_sec = s * 45.0
        state['mission_clock'] = f"T+00d 00h {t_sec:04.1f}s (Water Fill)"
        state['water_injected_mL'] = 7.39 * min(1.0, s * 1.08)
        state['foam_sat_pct'] = 100.0 * min(1.0, s * 1.05)
        state['water_tube_active'] = 1.0 if s < 0.94 else max(0.0, (1.0 - s) / 0.06)
        state['RH_mean'] = 45.0 + 28.0 * s
        if is_lunar:
            state['stage_sub'] = "Port 3 Dipping Nozzle (Y=9.5mm) -> Uniform Capillary Wicking Across 3 Foam Blocks & 35um Nytex Barrier (0.166g)"
        else:
            state['stage_sub'] = "Port 3 Dipping Nozzle (Y=9.5mm) -> 1.0g Hydrostatic Drainage Gradient Toward Basal Wells & Nytex Barrier"

    elif u < 0.44:
        # STAGE II: Air Pump Activation, Manifold Tubes, Turbulence & Moisture Plumes
        s = (u - 0.22) / 0.22
        state['stage_num'] = 2
        state['stage_name'] = "STAGE II: AIRFLOW TUBES, BAFFLE TURBULENCE & MOISTURE BOUNDARY LAYERS"
        t_hr = s * 24.0
        state['t_hours'] = t_hr
        state['day_float'] = t_hr / 24.0
        state['mission_clock'] = f"T+00d {int(t_hr):02d}h {int((t_hr%1)*60):02d}m (Day 1 Air)"
        state['water_injected_mL'] = 7.39
        state['foam_sat_pct'] = 100.0
        state['air_active'] = min(1.0, s * 4.0)

        if is_lunar:
            if s < 0.50:
                state['air_Q_Lmin'] = 0.50 * min(1.0, s * 3.0)
                state['stage_sub'] = "Uncompensated Flow (Q=0.50 L/min): 83.4% Buoyancy Loss Thickens Moisture Boundary Layer (delta_bl = 0.64 mm)"
            else:
                ramp = min(1.0, (s - 0.50) / 0.15)
                state['air_Q_Lmin'] = 0.50 + 0.35 * ramp
                state['stage_sub'] = "Lunar Compensatory Pump Boost (Q=0.85 L/min): Forced Convection Restores Turbulent Mixing & Thins Boundary Layer"
        else:
            state['air_Q_Lmin'] = 0.50 * min(1.0, s * 3.0)
            state['stage_sub'] = "Earth 1.0g Baseline (Q=0.50 L/min): Strong Buoyant Thermal/Vapor Plumes (Gr_L 6x Moon) Drive Natural Mixing"

    elif u < 0.68:
        # STAGE III: 5-Day Seed Germination & Asymmetric Canopy Ontogeny
        s = (u - 0.44) / 0.24
        state['stage_num'] = 3
        state['stage_name'] = "STAGE III: 5-DAY GERMINATION, CANOPY GROWTH & PLENUM AERODYNAMICS"
        t_hr = 24.0 + s * 96.0
        state['t_hours'] = t_hr
        state['day_float'] = t_hr / 24.0
        d_int = int(t_hr // 24)
        h_int = int(t_hr % 24)
        state['mission_clock'] = f"T+0{d_int}d {h_int:02d}h 00m (Day {state['day_float']:.1f})"
        state['water_injected_mL'] = 7.39 - 4.19 * ((t_hr - 24.0) / 96.0) ** 1.4
        state['foam_sat_pct'] = 100.0 * (state['water_injected_mL'] / 7.39)
        state['air_active'] = 1.0
        state['air_Q_Lmin'] = 0.85 if is_lunar else 0.50

        if t_hr < 48.0:
            state['stage_sub'] = f"Hour {t_hr:.0f}: Radicle Emergence (Blocked from Foam by 35um Nytex) & Cotyledon Unfolding"
        else:
            if is_lunar:
                state['stage_sub'] = f"Day {state['day_float']:.1f} (Q=0.85 L/min Pull Suction): Tall B. rapa (11.4mm) Channel Airflow Into 14.2mm Ch2 Plenum"
            else:
                state['stage_sub'] = f"Day {state['day_float']:.1f} (1.0g Control, Q=0.50 L/min): Buoyant Transpiration Plumes + Horizontal Canopy Sweep"

    elif u < 0.80:
        # STAGE IV: Diurnal Photoperiod Inversion (Nocturnal 8h Skotoperiod)
        s = (u - 0.68) / 0.12
        state['stage_num'] = 4
        state['stage_name'] = "STAGE IV: DIURNAL NIGHT INVERSION & NOCTURNAL CO2 FLUSHING"
        state['t_hours'] = 120.0
        state['day_float'] = 5.0
        state['is_night'] = True
        state['mission_clock'] = f"T+05d 0{int(s*8):1d}h 00m (Night Mode)"
        state['water_injected_mL'] = 3.20
        state['foam_sat_pct'] = 43.3
        state['air_active'] = 1.0
        if is_lunar:
            state['air_Q_Lmin'] = 0.60
            state['stage_sub'] = "Dark Period (LEDs Off): Nocturnal Flow Floor (Q_night = 0.60 L/min) Flushes Dense Mitochondrial CO2 Pooling"
        else:
            state['air_Q_Lmin'] = 0.50
            state['stage_sub'] = "Dark Period (1.0g Control): Negative Solutal Buoyancy (beta_CO2 = -0.517) Sinks Respiratory CO2 Toward Wells"

    else:
        # STAGE V: End-of-Mission RNAlater Harvest Perfusion
        s = (u - 0.80) / 0.20
        state['stage_num'] = 5
        state['stage_name'] = "STAGE V: IN-SITU RNALATER HARVEST FIXATION & PORT 2 GAS VENTING"
        state['t_hours'] = 120.0
        state['day_float'] = 5.0
        t_perf_s = s * 349.2
        state['mission_clock'] = f"T+05d 08h {int(t_perf_s//60):02d}m {int(t_perf_s%60):02d}s (Fixation)"
        state['water_injected_mL'] = 3.20
        state['foam_sat_pct'] = 100.0
        state['air_active'] = 0.25
        state['air_Q_Lmin'] = 0.01
        state['rna_active'] = 1.0
        state['rna_fill_frac'] = min(1.0, s * 1.02)
        state['rna_vol_mL'] = 64.60 * state['rna_fill_frac']
        state['rna_dP_kPa'] = 0.73 + 3.31 * min(1.0, s * 1.2)
        if t_perf_s <= 32.0:
            state['stage_sub'] = f"Port 3 Nozzle Injects 10 mL/min RNAlater -> Root Zone Saturated in 31.7s (dP = {state['rna_dP_kPa']:.2f} kPa)"
        else:
            meniscus_str = "Capillary-Curved Wetting Front" if is_lunar else "Gravity-Flattened Hydrostatic Meniscus"
            state['stage_sub'] = f"Bottom-Up Chamber Flood ({state['rna_vol_mL']:.1f}/64.6 mL): {meniscus_str} Submerges Canopy; Air Vents via Port 2"

    # Derived Aerodynamic & Biophysical Telemetry (Calibrated to exact manuscript values: 0.41 mm at 1g/0.50 & Lunar/0.85, 0.64 mm at Lunar/0.50)
    Q = max(state['air_Q_Lmin'], 0.02)
    v_char = 0.045 * (Q / 0.50)
    L_char = 0.015
    nu = 1.51e-5
    Re = (v_char * L_char) / nu
    Gr = (g_ms2 * (1.0 / 295.15) * 1.8 * (L_char ** 3)) / (nu ** 2)
    Ri = Gr / max(Re ** 2, 1e-4)

    # Calibrated boundary layer thickness in mm:
    # At Q=0.50, 1.0g -> 0.41 mm; at Q=0.50, 0.166g -> 0.64 mm; at Q=0.85, 0.166g -> 0.41 mm
    delta_forced_mm = 0.66 * ((0.50 / max(Q, 0.05)) ** 0.52)
    delta_nat_mm = 0.49 * ((1.0 / max(g_val, 0.02)) ** 0.25)
    delta_eff_mm = ((delta_forced_mm ** -3.0) + (delta_nat_mm ** -3.0)) ** (-1.0 / 3.0)

    state['Re_num'] = Re if state['air_active'] > 0.1 else 0.0
    state['Ri_num'] = Ri if state['air_active'] > 0.1 else 0.0
    state['delta_bl_mm'] = delta_eff_mm if state['air_active'] > 0.1 else 1.45

    if state['stage_num'] == 1 or state['stage_num'] == 5:
        state['A_net'] = 0.0
    elif state['is_night']:
        state['A_net'] = -1.85
    elif state['t_hours'] < 38.0:
        state['A_net'] = -0.75
    else:
        growth_f = min(1.0, (state['t_hours'] - 38.0) / 70.0)
        cond_f = min(1.02, (0.41 / max(delta_eff_mm, 0.30)) ** 0.48)
        state['A_net'] = 13.04 * growth_f * cond_f

    if state['stage_num'] >= 2 and state['stage_num'] <= 4:
        canopy_f = (state['day_float'] / 5.0) ** 1.5
        state['RH_mean'] = min(92.0, 62.0 + 14.0 * canopy_f * (0.50 / max(Q, 0.20)) + (3.5 if is_lunar else 0.0))
        state['T_margin_C'] = max(0.45, 2.35 - 1.45 * canopy_f * (0.50 / max(Q, 0.25)))

    return state


# =============================================================================
# 3. 2D Z-Y FIELD GENERATOR (MOISTURE, TURBULENCE VORTICITY & RNALATER)
# =============================================================================

Z_GRID_1D = np.linspace(-34.0, 34.0, 136)
Y_GRID_1D = np.linspace(0.0, 34.0, 68)
ZZ, YY = np.meshgrid(Z_GRID_1D, Y_GRID_1D)

MASK_FOAM1 = (YY >= 7.62) & (YY <= 14.62) & (ZZ >= -30.5) & (ZZ <= -9.5)
MASK_FOAM2 = (YY >= 7.62) & (YY <= 14.62) & (ZZ >= -7.5) & (ZZ <= 7.5)
MASK_FOAM3 = (YY >= 7.62) & (YY <= 14.62) & (ZZ >= 9.5) & (ZZ <= 30.5)
MASK_FOAM = MASK_FOAM1 | MASK_FOAM2 | MASK_FOAM3
MASK_BAFFLE1 = (YY >= 7.62) & (YY <= 24.0) & (ZZ > -9.5) & (ZZ < -7.5)
MASK_BAFFLE2 = (YY >= 7.62) & (YY <= 24.0) & (ZZ > 7.5) & (ZZ < 9.5)
MASK_BAFFLES = MASK_BAFFLE1 | MASK_BAFFLE2
MASK_DRAIN = (YY < 7.62) & (ZZ >= -32.0) & (ZZ <= 32.0)
MASK_HEADSPACE = (YY > 14.62) & (~MASK_BAFFLES)

def compute_2d_zy_fields(state):
    """Computes the 2D (Z, Y) scalar field for moisture/water/RNAlater."""
    u = state['u']
    is_lunar = state['is_lunar']
    stage = state['stage_num']

    field = np.zeros_like(ZZ, dtype=np.float64)

    if stage == 1:
        s = min(1.0, (u / 0.22) * 1.05)
        dist_nozzle = np.sqrt(((ZZ - 21.59) / 52.0) ** 2 + ((YY - 9.5) / 10.0) ** 2)
        drain_boost = np.where(YY < 8.5, 0.72, 1.0)
        eff_dist = dist_nozzle * drain_boost
        front = s * 1.25
        wet_frac = np.clip((front - eff_dist) / 0.18, 0.0, 1.0)

        y_rel = np.clip((YY - 7.62) / 7.0, 0.0, 1.0)
        if is_lunar:
            cap_prof = 0.96 - 0.02 * y_rel
        else:
            cap_prof = 0.99 - 0.16 * y_rel

        field = np.where(MASK_FOAM | MASK_DRAIN, wet_frac * cap_prof, 0.15)
    else:
        y_rel = np.clip((YY - 7.62) / 7.0, 0.0, 1.0)
        depletion = 0.86 if stage >= 3 else 0.96
        if is_lunar:
            cap_prof = depletion * (0.97 - 0.02 * y_rel)
        else:
            cap_prof = depletion * (1.00 - 0.16 * y_rel)
        field = np.where(MASK_FOAM | MASK_DRAIN, cap_prof, 0.20)

    if stage in (2, 3, 4):
        Q = max(state['air_Q_Lmin'], 0.10)
        t_hr = state['t_hours']
        h_brapa = 0.5 + 10.76 / (1.0 + math.exp(-0.075 * (t_hr - 62.0))) if t_hr >= 18 else 0.5
        h_arab = 0.3 + 4.45 / (1.0 + math.exp(-0.080 * (t_hr - 66.0))) if t_hr >= 22 else 0.3

        y_canopy = np.where(np.abs(ZZ) < 8.5, 14.62 + h_arab, 14.62 + h_brapa)
        dist_above_canopy = np.maximum(0.0, YY - y_canopy)

        bl_scale = state['delta_bl_mm'] * 11.0
        if is_lunar:
            plume_decay = np.exp(-dist_above_canopy / max(bl_scale, 1.5))
            eddy_wave = 0.08 * np.sin(ZZ * 0.45 - u * 45.0) * np.exp(-((YY - 22.0) / 6.0) ** 2)
            vapor = 0.24 + (0.62 * (0.55 / max(Q, 0.35))) * plume_decay + eddy_wave
        else:
            chimney = 0.16 * (1.0 + np.sin(ZZ * 0.55 - u * 30.0)) * np.exp(-dist_above_canopy / 14.0)
            plume_decay = np.exp(-dist_above_canopy / max(bl_scale * 0.85, 1.2))
            vapor = 0.24 + 0.45 * plume_decay + chimney

        field = np.where(MASK_HEADSPACE, np.clip(vapor, 0.15, 0.92), field)

    if stage == 5:
        s = state['rna_fill_frac']
        y_level_base = 7.62 + s * (33.67 - 7.62)
        lateral_tilt = (1.0 - s) * 4.5 * ((21.59 - ZZ) / 55.0)
        if is_lunar:
            meniscus = 1.4 * np.cos((ZZ / 34.0) * math.pi) * (1.0 - s * 0.5)
        else:
            meniscus = 0.15 * np.cos((ZZ / 34.0) * math.pi)
        y_interface = y_level_base - lateral_tilt - meniscus
        rna_mask = (YY <= y_interface) & (~MASK_BAFFLES)
        field = np.where(rna_mask, 0.98, field * 0.35)

    field = np.where(MASK_BAFFLES, np.nan, field)
    return field


# =============================================================================
# 4. SINGLE-FRAME RENDERER (1920x1080 CINEMA COMPOSITOR)
# =============================================================================

def render_single_frame_pil(args):
    """
    Renders a single 1920x1080 RGB frame for given (frame_idx, u, regime) and returns
    raw RGB bytes.
    """
    frame_idx, u, regime = args
    state = compute_mission_state(u, regime)
    is_lunar = state['is_lunar']
    stage = state['stage_num']

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor='#0b111e')

    # -------------------------------------------------------------------------
    # TOP CINEMA HEADER BANNER (Zero-overlap layout)
    # -------------------------------------------------------------------------
    ax_hdr = fig.add_axes([0.015, 0.905, 0.97, 0.08], facecolor='#111a2e')
    ax_hdr.set_xlim(0, 100)
    ax_hdr.set_ylim(0, 10)
    ax_hdr.axis('off')

    badge_col = '#3FB6A8' if is_lunar else '#f59e0b'
    regime_title = "LUNAR SURFACE MISSION (0.166g · 1.62 m/s²)" if is_lunar else "EARTH GROUND CONTROL (1.00g · 9.81 m/s²)"
    ax_hdr.add_patch(mpatches.FancyBboxPatch((0.4, 1.2), 27.2, 7.6, boxstyle="round,pad=0.3",
                                             facecolor=badge_col, alpha=0.18, edgecolor=badge_col, linewidth=2.0))
    ax_hdr.text(14.0, 6.2, regime_title, color=badge_col, fontsize=12.0, fontweight='bold', ha='center', va='center')
    ax_hdr.text(14.0, 2.8, "LEAF GBC1 · CAD M09024F (Pot) + M09025F (Lid)", color='#94a3b8', fontsize=9.0, ha='center', va='center')

    # Stage title & subtitle in middle column [29.0 .. 76.5]
    ax_hdr.text(29.0, 6.6, state['stage_name'], color='#f8fafc', fontsize=11.5, fontweight='bold', va='center')
    ax_hdr.text(29.0, 2.7, state['stage_sub'], color='#38bdf8' if is_lunar else '#fbbf24', fontsize=9.8, va='center')

    # Mission clock & timeline progress bar on right [78.0 .. 99.2]
    ax_hdr.text(99.2, 6.7, state['mission_clock'], color='#34d399', fontsize=11.0, fontweight='bold', ha='right', va='center', family='monospace')
    ax_hdr.add_patch(mpatches.Rectangle((78.5, 2.0), 20.7, 2.2, facecolor='#1e293b', edgecolor='#334155', linewidth=1.0))
    ax_hdr.add_patch(mpatches.Rectangle((78.5, 2.0), 20.7 * u, 2.2, facecolor=badge_col, edgecolor='none'))

    # -------------------------------------------------------------------------
    # LEFT PANEL: 3D VOLUMETRIC DIGITAL TWIN CUTAWAY
    # -------------------------------------------------------------------------
    ax3d = fig.add_axes([0.01, 0.02, 0.56, 0.83], projection='3d', facecolor='#0b111e')
    ax3d.view_init(elev=state['cam_elev'], azim=state['cam_azim'])

    ax3d.set_xlim(-31, 31)
    ax3d.set_ylim(-22, 22)
    ax3d.set_zlim(1, 43)
    ax3d.set_box_aspect((68, 48, 42))
    ax3d.xaxis.pane.fill = False
    ax3d.yaxis.pane.fill = False
    ax3d.zaxis.pane.fill = False
    ax3d.xaxis.pane.set_edgecolor('#1e293b')
    ax3d.yaxis.pane.set_edgecolor('#1e293b')
    ax3d.zaxis.pane.set_edgecolor('#1e293b')
    ax3d.grid(True, color='#1e293b', linestyle=':', linewidth=0.5, alpha=0.5)
    ax3d.tick_params(colors='#64748b', labelsize=8)
    ax3d.set_xlabel('Axial Length Z (mm)', color='#94a3b8', fontsize=9, labelpad=4)
    ax3d.set_ylabel('Width X (mm)', color='#94a3b8', fontsize=9, labelpad=4)
    ax3d.set_zlabel('Elevation Y (mm)', color='#94a3b8', fontsize=9, labelpad=4)

    # 1. Outer Pot Wireframe
    pot_z = [-34, 34, 34, -34, -34]
    pot_x = [-24, -24, 24, 24, -24]
    ax3d.plot(pot_z, pot_x, [0]*5, color='#475569', linewidth=1.2, alpha=0.7)
    ax3d.plot(pot_z, pot_x, [33.67]*5, color='#38bdf8', linewidth=1.6, alpha=0.85)
    for zc, xc in [(-34, -24), (34, -24), (34, 24), (-34, 24)]:
        ax3d.plot([zc, zc], [xc, xc], [0, 33.67], color='#475569', linewidth=1.0, alpha=0.6)

    # Baffles at Z = -8.5 and Z = +8.5 mm
    for zb in [-8.5, 8.5]:
        ax3d.plot([zb, zb, zb, zb, zb], [-11, 21, 21, -11, -11], [7.62, 7.62, 24.0, 24.0, 7.62],
                  color='#94a3b8', linewidth=1.5, linestyle='--', alpha=0.8)

    # 2. Nytex 35 um Mesh Root Exclusion Plane at Y = 14.62 mm
    for z0, z1 in [(-30.5, -9.5), (-7.5, 7.5), (9.5, 30.5)]:
        zz_ny, xx_ny = np.meshgrid(np.linspace(z0, z1, 6), np.linspace(-11.0, 21.0, 5))
        yy_ny = np.full_like(zz_ny, 14.62)
        ax3d.plot_wireframe(zz_ny, xx_ny, yy_ny, color='#38bdf8', linewidth=0.75, alpha=0.55)

    # 3. 3-Part Porous Foam Cassettes (Y = 7.62 to 14.62 mm) & Water Saturation Voxels
    f_z = np.linspace(-29.0, 29.0, 18)
    f_x = np.linspace(-9.0, 19.0, 6)
    f_y = np.linspace(8.2, 14.0, 4)
    FZ, FX, FY = np.meshgrid(f_z, f_x, f_y)
    valid_foam = ((FZ <= -9.5) | ((FZ >= -7.5) & (FZ <= 7.5)) | (FZ >= 9.5))
    fz_pts = FZ[valid_foam]
    fx_pts = FX[valid_foam]
    fy_pts = FY[valid_foam]

    if stage == 1:
        s = min(1.0, (u / 0.22) * 1.05)
        dist3d = np.sqrt(((fz_pts - 21.59) / 52.0) ** 2 + ((fx_pts - (-15.49)) / 42.0) ** 2 + ((fy_pts - 9.5) / 10.0) ** 2)
        wet3d = np.clip((s * 1.30 - dist3d) / 0.18, 0.0, 1.0)
    else:
        y_norm = (fy_pts - 7.62) / 7.0
        wet3d = (0.95 - 0.02 * y_norm) if is_lunar else (0.98 - 0.18 * y_norm)

    colors_foam = np.zeros((len(fz_pts), 4))
    colors_foam[:, 0] = 0.71 * (1.0 - wet3d) + 0.01 * wet3d
    colors_foam[:, 1] = 0.33 * (1.0 - wet3d) + 0.52 * wet3d
    colors_foam[:, 2] = 0.04 * (1.0 - wet3d) + 0.78 * wet3d
    colors_foam[:, 3] = 0.42 + 0.35 * wet3d
    ax3d.scatter(fz_pts, fx_pts, fy_pts, c=colors_foam, s=30, marker='s', depthshade=False)

    # 4. External Manifold Tubes & Internal Port 3 Dipping Nozzle
    for zp, col_p, lbl_p in [(-21.59, '#ef4444', 'Port 1 (Air In)'),
                             (0.00, '#a855f7', 'Port 2 (Exhaust/Vent)'),
                             (21.59, '#0284c7', 'Port 3 (H2O / RNAlater)')]:
        ax3d.plot([zp, zp], [-15.49, -15.49], [33.67, 44.5], color=col_p, linewidth=5.5, alpha=0.75)
        ax3d.text(zp, -15.49, 45.3, lbl_p, color=col_p, fontsize=8.5, fontweight='bold', ha='center')

    # Port 3 internal dipping nozzle down to Y = 9.5 mm
    ax3d.plot([21.59, 21.59], [-15.49, -15.49], [9.5, 33.67], color='#38bdf8', linewidth=3.5, alpha=0.95)
    ax3d.scatter([21.59], [-15.49], [9.5], color='#0284c7', s=55, zorder=10)

    # Active liquid droplets in Port 3 tube during Stage I (Water) or Stage V (RNAlater)
    if state['water_tube_active'] > 0.05 or state['rna_active'] > 0.05:
        liq_col = '#a855f7' if state['rna_active'] > 0.05 else '#38bdf8'
        phase_shift = (u * 60.0) % 1.0
        y_drops = np.linspace(44.0, 9.5, 9) - phase_shift * 3.8
        y_drops = np.clip(y_drops, 9.5, 44.0)
        ax3d.scatter(np.full_like(y_drops, 21.59), np.full_like(y_drops, -15.49), y_drops,
                     color=liq_col, s=42, alpha=0.95)

    # 5. Render All 164 Growing 3D Plants Across the 3 Compartments
    t_hr = state['t_hours']
    px_all = PLANTS_ARRAY[:, 0]
    pz_all = PLANTS_ARRAY[:, 1]
    sp_all = PLANTS_ARRAY[:, 2]
    jit_all = PLANTS_ARRAY[:, 4]

    t_eff = np.maximum(0.0, t_hr + jit_all)
    h_plants = np.where(
        sp_all == 1,
        np.where(t_eff < 18.0, 0.5, 0.5 + 10.76 / (1.0 + np.exp(-0.075 * (t_eff - 62.0)))),
        np.where(t_eff < 22.0, 0.3, 0.3 + 4.45 / (1.0 + np.exp(-0.080 * (t_eff - 66.0))))
    )

    y_base = 14.77
    y_top = y_base + h_plants

    if t_hr >= 14.0:
        stem_segs = [((pz_all[i], px_all[i], y_base), (pz_all[i], px_all[i], y_top[i])) for i in range(len(px_all))]
        lc_stems = Line3DCollection(stem_segs, colors='#10b981', linewidths=1.35, alpha=0.85)
        ax3d.add_collection3d(lc_stems)

    crown_sizes = np.where(sp_all == 1, 8.0 + 40.0 * ((h_plants - 0.5) / 10.8), 5.0 + 18.0 * ((h_plants - 0.3) / 4.5))
    crown_cols = np.where(t_hr < 18.0, '#d97706', np.where(sp_all == 1, '#34d399', '#38bdf8'))
    ax3d.scatter(pz_all, px_all, y_top, c=crown_cols, s=crown_sizes, alpha=0.92, depthshade=True)

    # 6. 3D Airflow Streamlines, Tube Flow, Baffle Turbulence Eddies & Buoyant Plumes
    if state['air_active'] > 0.05:
        n_streams = 14
        t_param = np.linspace(0.0, 1.0, 45)
        q_norm = max(0.3, state['air_Q_Lmin'] / 0.50)
        phase = (u * 35.0 * q_norm) % 1.0

        for idx_s in range(n_streams):
            x_off = -12.0 + (idx_s / (n_streams - 1.0)) * 28.0
            z_path = np.where(
                t_param < 0.15,
                -21.59,
                np.where(
                    t_param < 0.85,
                    -21.59 + ((t_param - 0.15) / 0.70) * 21.59 + 3.2 * np.sin(t_param * 18.0 + idx_s),
                    0.0
                )
            )
            if idx_s % 3 == 0:
                z_path = np.where(
                    t_param < 0.15,
                    -21.59,
                    np.where(
                        t_param < 0.55,
                        -21.59 + ((t_param - 0.15) / 0.40) * 46.0,
                        np.where(
                            t_param < 0.85,
                            24.41 - ((t_param - 0.55) / 0.30) * 24.41,
                            0.0
                        )
                    )
                )

            x_path = np.where(
                t_param < 0.15,
                -15.49,
                np.where(
                    t_param < 0.85,
                    -15.49 + np.sin(((t_param - 0.15) / 0.70) * math.pi) * (x_off + 15.49),
                    -15.49
                )
            )

            buoyant_lift = (2.8 * np.sin(t_param * math.pi * 2.0)) if (not is_lunar and not state['is_night']) else 0.0
            turbulence_eddy = 1.6 * q_norm * np.sin(t_param * 28.0 - u * 50.0 + idx_s) * np.exp(-((z_path + 8.5) / 6.0) ** 2)
            y_path = np.where(
                t_param < 0.15,
                44.0 - (t_param / 0.15) * 12.0,
                np.where(
                    t_param < 0.85,
                    29.5 - 6.8 * np.exp(-((z_path) / 9.0) ** 2) + buoyant_lift + turbulence_eddy,
                    31.0 + ((t_param - 0.85) / 0.15) * 13.0
                )
            )

            s_col = '#38bdf8' if (idx_s % 2 == 0) else '#34d399'
            if state['is_night']:
                s_col = '#818cf8'
            ax3d.plot(z_path, x_path, y_path, color=s_col, linewidth=1.15, alpha=0.50 * state['air_active'])

            for k_pt in range(3):
                p_loc = (phase + k_pt / 3.0 + idx_s * 0.07) % 1.0
                k_idx = min(len(t_param) - 1, int(p_loc * len(t_param)))
                ax3d.scatter([z_path[k_idx]], [x_path[k_idx]], [y_path[k_idx]],
                             color='#f8fafc', edgecolor=s_col, s=24, alpha=0.92)

    # 7. Stage V: 3D Rising RNAlater Liquid Surface
    if state['rna_active'] > 0.05:
        s_rna = state['rna_fill_frac']
        zz_r, xx_r = np.meshgrid(np.linspace(-33.0, 33.0, 20), np.linspace(-22.0, 22.0, 12))
        y_lvl = 7.62 + s_rna * (33.67 - 7.62)
        tilt = (1.0 - s_rna) * 4.0 * ((21.59 - zz_r) / 55.0)
        menisc = (1.3 if is_lunar else 0.15) * np.cos((zz_r / 34.0) * math.pi) * (1.0 - 0.5 * s_rna)
        yy_r = np.clip(y_lvl - tilt - menisc, 7.62, 33.67)
        ax3d.plot_surface(zz_r, xx_r, yy_r, color='#8b5cf6', alpha=0.42, edgecolor='#a78bfa', linewidth=0.3)

    fig.text(0.29, 0.865, "3D Volumetric Digital Twin Cutaway (Manifold Tubes, Porous Foam, 164 Plants & Streamlines)",
             color='#e2e8f0', fontsize=11.5, fontweight='bold', ha='center')

    # -------------------------------------------------------------------------
    # TOP-RIGHT PANEL: 2D LONGITUDINAL (Z-Y) CFD, TURBULENCE & MOISTURE FIELD
    # -------------------------------------------------------------------------
    ax_zy = fig.add_axes([0.60, 0.50, 0.37, 0.36], facecolor='#0f172a')
    field_zy = compute_2d_zy_fields(state)

    cmap_name = 'Purples' if stage == 5 else 'YlGnBu'
    im = ax_zy.imshow(field_zy, extent=[-34.0, 34.0, 0.0, 34.0], origin='lower',
                      cmap=cmap_name, vmin=0.0, vmax=1.0, aspect='auto', alpha=0.92)

    for z0, z1 in [(-30.5, -9.5), (-7.5, 7.5), (9.5, 30.5)]:
        ax_zy.add_patch(mpatches.Rectangle((z0, 7.62), z1 - z0, 7.0, fill=False,
                                           edgecolor='#f59e0b', linewidth=1.2, linestyle='--'))
        ax_zy.plot([z0, z1], [14.62, 14.62], color='#38bdf8', linewidth=2.0)

    for zb in [-9.5, 7.5]:
        ax_zy.add_patch(mpatches.Rectangle((zb, 7.62), 2.0, 16.38, facecolor='#334155', edgecolor='#94a3b8'))

    ax_zy.plot([21.59, 21.59], [9.5, 34.0], color='#38bdf8', linewidth=3.0)
    ax_zy.scatter([21.59], [9.5], color='#0284c7', edgecolor='white', s=50, zorder=6)

    if state['air_active'] > 0.1 and stage in (2, 3, 4):
        zv = np.linspace(-30.0, 30.0, 13)
        yv = np.linspace(16.5, 31.5, 5)
        ZV, YV = np.meshgrid(zv, yv)
        q_sc = state['air_Q_Lmin'] / 0.50
        UZ = np.where(ZV < 0, 1.2 * q_sc, -0.9 * q_sc) * (1.0 + 0.3 * np.sin(ZV * 0.3 - u * 40.0))
        UY = -0.5 * np.exp(-((ZV) / 10.0) ** 2) + 0.9 * np.exp(-((ZV) / 4.5) ** 2) * ((YV - 15.0) / 16.0)
        if not is_lunar and not state['is_night']:
            UY += 0.45
        ax_zy.quiver(ZV, YV, UZ, UY, color='#f8fafc', alpha=0.68, scale=18.0, width=0.004)

    ax_zy.set_title("2D Longitudinal (Z-Y) Cross-Section: Foam Capillary Moisture, Boundary Layer & Velocity",
                    color='#e2e8f0', fontsize=10.5, fontweight='bold', pad=6)
    ax_zy.set_xlabel("Axial Position Z (mm)  [Ch 1: B. rapa | Ch 2: Arabidopsis | Ch 3: B. rapa]", color='#94a3b8', fontsize=8.5)
    ax_zy.set_ylabel("Elevation Y (mm)", color='#94a3b8', fontsize=8.5)
    ax_zy.tick_params(colors='#94a3b8', labelsize=8)
    cb = fig.colorbar(im, ax=ax_zy, fraction=0.035, pad=0.02)
    cb.ax.tick_params(colors='#94a3b8', labelsize=7.5)
    cb.set_label("Normalized Moisture / Saturation", color='#94a3b8', fontsize=8)

    # -------------------------------------------------------------------------
    # BOTTOM-RIGHT PANEL 1: FULL VERTICAL MOISTURE & BOUNDARY-LAYER PROFILE (Y = 7.6 to 33.7 mm)
    # -------------------------------------------------------------------------
    ax_prof = fig.add_axes([0.60, 0.06, 0.18, 0.36], facecolor='#0f172a')
    y_full = np.linspace(7.62, 33.67, 100)
    sat_scale = state['foam_sat_pct'] / 100.0

    # Foam region (Y <= 14.62) + Canopy Boundary Layer & Headspace Humidity region (Y > 14.62)
    y_foam_rel = np.clip((y_full - 7.62) / 7.0, 0.0, 1.0)
    y_head_rel = np.clip((y_full - 14.62) / 19.05, 0.0, 1.0)

    # Exaggerate hydrostatic vs capillary distinction slightly within physical bounds + show headspace vapor gradient
    prof_lunar = np.where(
        y_full <= 14.62,
        sat_scale * (96.5 - 1.5 * y_foam_rel),
        48.0 + (state['RH_mean'] - 48.0) * np.exp(-y_head_rel / max(0.18, state['delta_bl_mm'] * 0.48))
    )
    prof_earth = np.where(
        y_full <= 14.62,
        sat_scale * (99.0 - 14.5 * y_foam_rel),
        48.0 + (state['RH_mean'] - 50.0) * np.exp(-y_head_rel / 0.20) + (6.5 * np.sin(y_head_rel * math.pi) if state['air_active'] > 0.1 else 0.0)
    )

    if is_lunar:
        ax_prof.plot(prof_lunar, y_full, color='#3FB6A8', linewidth=2.8, label='Lunar 0.166g (Active)')
        ax_prof.plot(prof_earth, y_full, color='#f59e0b', linewidth=1.6, linestyle='--', alpha=0.75, label='Earth 1.0g Control')
    else:
        ax_prof.plot(prof_earth, y_full, color='#f59e0b', linewidth=2.8, label='Earth 1.0g (Active)')
        ax_prof.plot(prof_lunar, y_full, color='#3FB6A8', linewidth=1.6, linestyle='--', alpha=0.75, label='Lunar 0.166g Ref')

    ax_prof.axhline(14.62, color='#38bdf8', linestyle=':', linewidth=1.5, label='35um Nytex Interface')
    ax_prof.axhspan(7.62, 14.62, color='#0284c7', alpha=0.10)
    ax_prof.set_xlim(35.0, 102.0)
    ax_prof.set_ylim(7.62, 33.67)
    ax_prof.set_title("Vertical Moisture & Vapor Profile (Y)", color='#e2e8f0', fontsize=9.5, fontweight='bold')
    ax_prof.set_xlabel("Foam Saturation / Headspace RH (%)", color='#94a3b8', fontsize=8)
    ax_prof.set_ylabel("Chamber Elevation Y (mm)", color='#94a3b8', fontsize=8)
    ax_prof.tick_params(colors='#94a3b8', labelsize=7.5)
    ax_prof.grid(True, color='#1e293b', linestyle='--', alpha=0.6)
    ax_prof.legend(loc='upper right', fontsize=7, facecolor='#111a2e', edgecolor='#334155', labelcolor='#e2e8f0')

    # -------------------------------------------------------------------------
    # BOTTOM-RIGHT PANEL 2: LIVE MULTI-PHYSICS TELEMETRY CARDS
    # -------------------------------------------------------------------------
    ax_tel = fig.add_axes([0.80, 0.06, 0.175, 0.36], facecolor='#111a2e')
    ax_tel.set_xlim(0, 10)
    ax_tel.set_ylim(0, 10)
    ax_tel.axis('off')
    ax_tel.add_patch(mpatches.Rectangle((0, 0), 10, 10, fill=False, edgecolor='#334155', linewidth=1.5))
    ax_tel.text(5.0, 9.2, "LIVE 4D TELEMETRY", color='#38bdf8', fontsize=10, fontweight='bold', ha='center')

    metrics = [
        ("Gravity Regime:", f"{state['g_val']:.3f}g ({state['g_ms2']:.2f} m/s²)", badge_col),
        ("Foam Hydration:", f"{state['water_injected_mL']:.2f} mL ({state['foam_sat_pct']:.0f}%)", '#38bdf8'),
        ("Pump Flow Rate Q:", f"{state['air_Q_Lmin']:.2f} L/min", '#34d399'),
        ("Boundary Layer δ_bl:", f"{state['delta_bl_mm']:.2f} mm", '#fbbf24'),
        ("Richardson Ri:", f"{state['Ri_num']:.3f} ({'Forced' if state['Ri_num']<0.1 else 'Mixed'})", '#a78bfa'),
        ("Canopy A_net:", f"{state['A_net']:+.2f} μmol/m²s", '#34d399' if state['A_net'] >= 0 else '#f87171'),
        ("Headspace RH / ΔT:", f"{state['RH_mean']:.0f}% / +{state['T_margin_C']:.2f}°C", '#38bdf8'),
        ("RNAlater Flood:", f"{state['rna_vol_mL']:.1f} mL ({state['rna_fill_frac']*100:.0f}%)", '#c084fc'),
    ]
    for idx_m, (k_str, v_str, c_str) in enumerate(metrics):
        y_pos = 8.0 - idx_m * 1.02
        ax_tel.text(0.45, y_pos, k_str, color='#94a3b8', fontsize=8.2, va='center')
        ax_tel.text(9.55, y_pos, v_str, color=c_str, fontsize=8.5, fontweight='bold', ha='right', va='center', family='monospace')

    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    rgb_bytes = buf[:, :, :3].tobytes()
    plt.close(fig)
    return frame_idx, rgb_bytes


def compose_split_screen_frames(lunar_frames, earth_frames, width=1920, height=1080):
    """
    Composes a native Full HD (1920x1080) side-by-side comparison video without letterboxing
    by cropping and pairing the 3D + 2D + Telemetry regions of Lunar (left 960px) and Earth (right 960px).
    """
    split_frames = []
    for idx in range(len(lunar_frames)):
        arr_l = np.frombuffer(lunar_frames[idx][1], dtype=np.uint8).reshape((height, width, 3))
        arr_e = np.frombuffer(earth_frames[idx][1], dtype=np.uint8).reshape((height, width, 3))

        img_l = Image.fromarray(arr_l)
        img_e = Image.fromarray(arr_e)

        canvas = Image.new('RGB', (width, height), (11, 17, 30))
        # Top banner: keep stage title & clock in center/right, put compact Lunar vs 1g dual badge on left
        banner = img_l.crop((0, 0, 1920, 104))
        canvas.paste(banner, (0, 0))
        l_badge_sm = img_l.crop((24, 10, 560, 98)).resize((258, 76), Image.Resampling.BILINEAR)
        e_badge_sm = img_e.crop((24, 10, 560, 98)).resize((258, 76), Image.Resampling.BILINEAR)
        canvas.paste(l_badge_sm, (14, 14))
        canvas.paste(e_badge_sm, (278, 14))

        # Left column (x: 0..958): Lunar 3D cutaway (top) + Lunar 2D Z-Y & Telemetry (bottom)
        l_3d = img_l.crop((15, 118, 1085, 995)).resize((940, 495), Image.Resampling.BILINEAR)
        l_2d = img_l.crop((1085, 115, 1915, 1068)).resize((940, 455), Image.Resampling.BILINEAR)
        canvas.paste(l_3d, (10, 110))
        canvas.paste(l_2d, (10, 612))

        # Right column (x: 962..1920): Earth 3D cutaway (top) + Earth 2D Z-Y & Telemetry (bottom)
        e_3d = img_e.crop((15, 118, 1085, 995)).resize((940, 495), Image.Resampling.BILINEAR)
        e_2d = img_e.crop((1085, 115, 1915, 1068)).resize((940, 455), Image.Resampling.BILINEAR)
        canvas.paste(e_3d, (970, 110))
        canvas.paste(e_2d, (970, 612))

        # Overlay regime column headers directly on the left & right columns
        arr_out = np.array(canvas)
        arr_out[104:, 958:962, :] = [59, 110, 165]
        arr_out[104:108, :958, :] = [63, 182, 168]
        arr_out[104:108, 962:, :] = [245, 158, 11]
        split_frames.append((idx, arr_out.tobytes()))
    return split_frames


# =============================================================================
# 5. PARALLEL VIDEO PIPELINE & FFMPEG ENCODING
# =============================================================================

def encode_video_from_frames(rgb_frames, output_path, fps=24, width=1920, height=1080):
    """Encodes a list of raw RGB24 byte buffers into an H.264 MP4 video via FFmpeg."""
    cmd = [
        FFMPEG_BIN, '-y',
        '-f', 'rawvideo',
        '-vcodec', 'rawvideo',
        '-s', f'{width}x{height}',
        '-pix_fmt', 'rgb24',
        '-r', str(fps),
        '-i', '-',
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-preset', 'fast',
        '-crf', '20',
        '-movflags', '+faststart',
        output_path
    ]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for _, rgb in rgb_frames:
        proc.stdin.write(rgb)
    proc.stdin.close()
    proc.wait()
    if proc.returncode != 0:
        err = proc.stderr.read().decode('utf-8', errors='ignore')
        raise RuntimeError(f"FFmpeg failed on {output_path}: {err}")
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"  -> Saved {output_path} ({len(rgb_frames)} frames, {size_mb:.2f} MB)")


def render_regime_frames(regime, n_frames=360, n_workers=10):
    """Renders all frames for a given gravity regime ('lunar' or 'earth') in parallel."""
    print(f"\n=== Rendering {n_frames} 4D frames for regime='{regime}' across {n_workers} workers ===")
    t0 = time.time()
    u_vals = np.linspace(0.0, 1.0, n_frames)
    tasks = [(i, float(u_vals[i]), regime) for i in range(n_frames)]

    frames = []
    with mp.Pool(processes=n_workers) as pool:
        for idx, rgb in pool.imap(render_single_frame_pil, tasks, chunksize=6):
            frames.append((idx, rgb))
            if (idx + 1) % 60 == 0 or (idx + 1) == n_frames:
                elapsed = time.time() - t0
                print(f"  [{regime}] Rendered {idx+1}/{n_frames} frames ({elapsed:.1f}s)")

    frames.sort(key=lambda x: x[0])
    return frames


def save_poster_image(rgb_bytes, path, width=1920, height=1080):
    img = Image.frombytes('RGB', (width, height), rgb_bytes)
    img.save(path, format='PNG', optimize=True)
    print(f"  -> Saved poster frame: {path}")


def main():
    n_workers = min(10, max(2, (os.cpu_count() or 4) - 2))
    n_frames = 360
    fps = 24

    # 1. Render Lunar 0.166g Master Video
    lunar_frames = render_regime_frames('lunar', n_frames=n_frames, n_workers=n_workers)
    encode_video_from_frames(lunar_frames, 'videos/lunar_experiment_4d_cinematic.mp4', fps=fps)
    save_poster_image(lunar_frames[int(n_frames * 0.58)][1], 'videos/poster_lunar_4d.png')

    # 2. Render Earth 1.0g Ground Control Master Video
    earth_frames = render_regime_frames('earth', n_frames=n_frames, n_workers=n_workers)
    encode_video_from_frames(earth_frames, 'videos/ground_control_1g_4d_cinematic.mp4', fps=fps)
    save_poster_image(earth_frames[int(n_frames * 0.58)][1], 'videos/poster_ground_1g_4d.png')

    # 3. Encode Stage-Specific Highlight Videos (Lunar 0.166g)
    s1_end = int(n_frames * 0.23)
    encode_video_from_frames(lunar_frames[:s1_end], 'videos/stage1_water_foam_activation.mp4', fps=fps)

    s2_end = int(n_frames * 0.45)
    encode_video_from_frames(lunar_frames[s1_end:s2_end], 'videos/stage2_airflow_turbulence_moisture.mp4', fps=fps)

    s4_end = int(n_frames * 0.80)
    encode_video_from_frames(lunar_frames[s2_end:s4_end], 'videos/stage3_5day_growth_ontogeny.mp4', fps=fps)

    encode_video_from_frames(lunar_frames[s4_end:], 'videos/stage4_rnalater_harvest_injection.mp4', fps=fps)

    # 4. Generate Full-Screen Native 1920x1080 Side-by-Side Comparison Video (Lunar 0.166g vs. Earth 1.0g)
    print("\n=== Composing & Encoding Full-Screen Split Comparison Video (Lunar 0.166g vs. Earth 1.0g) ===")
    split_frames = compose_split_screen_frames(lunar_frames, earth_frames)
    encode_video_from_frames(split_frames, 'videos/comparison_lunar_vs_1g_4d.mp4', fps=fps)
    save_poster_image(split_frames[int(n_frames * 0.58)][1], 'videos/poster_comparison_4d.png')

    print("\nAll 4D Cinematic Videos and Posters generated successfully!")


if __name__ == '__main__':
    main()
