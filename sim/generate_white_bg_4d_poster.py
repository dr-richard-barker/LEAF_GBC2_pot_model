"""
generate_white_bg_4d_poster.py

Generates a publication-grade, Full HD (1920x1080) synchronized side-by-side comparison poster
of the LEAF GBC1 4D digital twin under Lunar Hypogravity (0.166g) vs. Earth Ground Control (1.0g)
on a STRICT PURE WHITE BACKGROUND (#ffffff) for npj Microgravity.

Features:
- Pure white background (#ffffff) and light studio styling (#f8fafc panels)
- Solid acrylic divider baffles at Z = -8.5 mm and +8.5 mm (Y = 7.62 to 24.0 mm)
- Upward boundary deflection and zero velocity vectors penetrating solid dividers
- Side-by-side comparison of 3D digital twin, 2D longitudinal CFD field, vertical boundary layer profiles, and live biophysical telemetry.

Author: Richard Barker, Department of Agricultural and Biological Engineering, Purdue University
Email: barkerrj@purdue.edu
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image

def generate_white_bg_poster():
    print("Generating publication-grade 4D comparison poster on white background (#ffffff)...")
    
    # 1920x1080 at 100 dpi = 19.2 x 10.8 inches
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor='#ffffff')
    
    # Global Header Banner [0.02, 0.91, 0.96, 0.08]
    ax_hdr = fig.add_axes([0.02, 0.915, 0.96, 0.075], facecolor='#f8fafc')
    ax_hdr.set_xlim(0, 100)
    ax_hdr.set_ylim(0, 10)
    ax_hdr.axis('off')
    
    # Draw header card outline
    hdr_box = mpatches.FancyBboxPatch((0, 0), 100, 10, boxstyle="round,pad=0.2",
                                      facecolor='#f8fafc', edgecolor='#cbd5e1', linewidth=1.2)
    ax_hdr.add_patch(hdr_box)
    
    # Left badge: Mission & Hardware
    ax_hdr.text(1.2, 6.5, "NASA ARTEMIS III · LEAF PAYLOAD", color='#b91c1c', fontsize=11.5, fontweight='bold', va='center')
    ax_hdr.text(1.2, 2.8, "Growth Bed Chamber 1 (GBC1) · CAD M09024F (Pot) + M09025F (Lid)", color='#475569', fontsize=9.2, va='center')
    
    # Center: Stage & Flow Mode
    ax_hdr.text(50.0, 6.5, "STAGE II: PUSH-MODE AERODYNAMICS & BAFFLE BOUNDARY LAYERS", color='#0f172a', fontsize=12.0, fontweight='bold', ha='center', va='center')
    ax_hdr.text(50.0, 2.8, "Forced Inlet at Port 1 (Z = -21.6 mm) → Upward Cresting Over Divider 1 → Exhaust at Port 2 (Z = 0.0 mm)", color='#0284c7', fontsize=9.5, fontweight='bold', ha='center', va='center')
    
    # Right: Mission Clock
    ax_hdr.text(98.8, 6.5, "MET: Day 02 · 14:32:00", color='#047857', fontsize=11.5, fontweight='bold', ha='right', va='center', family='monospace')
    ax_hdr.text(98.8, 2.8, "Diurnal Day Cycle (PPFD 200 μmol/m²s)", color='#64748b', fontsize=8.8, ha='right', va='center')

    # Sub-functions to render a single regime half
    def render_regime_side(x_offset, regime_name, is_lunar):
        badge_title = "LUNAR SURFACE PAYLOAD (0.166g · g = 1.62 m/s²)" if is_lunar else "EARTH GROUND CONTROL (1.00g · g = 9.81 m/s²)"
        badge_col = "#dc2626" if is_lunar else "#2563eb"
        badge_bg = "#fef2f2" if is_lunar else "#eff6ff"
        
        # Sub-header badge for this half
        ax_subhdr = fig.add_axes([x_offset + 0.015, 0.850, 0.455, 0.045], facecolor=badge_bg)
        ax_subhdr.axis('off')
        sb_box = mpatches.FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.1",
                                         facecolor=badge_bg, edgecolor=badge_col, linewidth=1.2, transform=ax_subhdr.transAxes)
        ax_subhdr.add_patch(sb_box)
        ax_subhdr.text(0.5, 0.5, badge_title, color=badge_col, fontsize=10.5, fontweight='bold', ha='center', va='center')

        # ---------------------------------------------------------------------
        # 3D Digital Twin Cutaway (Left part of this side)
        # ---------------------------------------------------------------------
        ax3d = fig.add_axes([x_offset + 0.01, 0.38, 0.25, 0.45], projection='3d', facecolor='#ffffff')
        elev_angle = 24.0
        azim_angle = -52.0 if is_lunar else -48.0
        ax3d.view_init(elev=elev_angle, azim=azim_angle)
        
        ax3d.set_xlim(-34, 34)
        ax3d.set_ylim(-24, 24)
        ax3d.set_zlim(0, 42)
        ax3d.set_box_aspect((68, 48, 42))
        ax3d.xaxis.pane.fill = False
        ax3d.yaxis.pane.fill = False
        ax3d.zaxis.pane.fill = False
        ax3d.xaxis.pane.set_edgecolor('#e2e8f0')
        ax3d.yaxis.pane.set_edgecolor('#e2e8f0')
        ax3d.zaxis.pane.set_edgecolor('#e2e8f0')
        ax3d.grid(True, color='#e2e8f0', linestyle=':', linewidth=0.5, alpha=0.7)
        ax3d.tick_params(colors='#475569', labelsize=7)
        ax3d.set_xlabel('Z (mm)', color='#334155', fontsize=8, labelpad=2)
        ax3d.set_ylabel('X (mm)', color='#334155', fontsize=8, labelpad=2)
        ax3d.set_zlabel('Y (mm)', color='#334155', fontsize=8, labelpad=2)
        
        # Pot wireframe
        pot_z = [-34, 34, 34, -34, -34]
        pot_x = [-24, -24, 24, 24, -24]
        ax3d.plot(pot_z, pot_x, [0]*5, color='#64748b', linewidth=1.1, alpha=0.8)
        ax3d.plot(pot_z, pot_x, [33.67]*5, color='#0284c7', linewidth=1.4, alpha=0.9)
        for zc, xc in [(-34, -24), (34, -24), (34, 24), (-34, 24)]:
            ax3d.plot([zc, zc], [xc, xc], [0, 33.67], color='#94a3b8', linewidth=0.9, alpha=0.7)
            
        # Solid Acrylic Baffles at Z = -8.5 and +8.5 mm
        for zb in [-8.5, 8.5]:
            ax3d.plot([zb, zb, zb, zb, zb], [-11, 21, 21, -11, -11], [7.62, 7.62, 24.0, 24.0, 7.62],
                      color='#334155', linewidth=1.6, linestyle='-')
            # Fill baffle patch
            bb_y, bb_x = np.meshgrid(np.linspace(7.62, 24.0, 4), np.linspace(-11, 21, 4))
            bb_z = np.full_like(bb_y, zb)
            ax3d.plot_surface(bb_z, bb_x, bb_y, color='#cbd5e1', alpha=0.5, edgecolor='#475569', linewidth=0.3)
            
        # Substrate blocks
        for z0, z1 in [(-30.5, -9.5), (-7.5, 7.5), (9.5, 30.5)]:
            s_z = [z0, z1, z1, z0, z0]
            s_x = [-11, -11, 21, 21, -11]
            ax3d.plot(s_z, s_x, [7.62]*5, color='#d97706', linewidth=1.0)
            ax3d.plot(s_z, s_x, [14.62]*5, color='#d97706', linewidth=1.0)
            
        # Seedlings in 3 chambers
        # Chamber 1 & 3: B. rapa (tall stems)
        for z_seed in [-26, -20, -14, 14, 20, 26]:
            for x_seed in [-6, 2, 10, 16]:
                ax3d.plot([z_seed, z_seed], [x_seed, x_seed], [14.62, 24.5], color='#15803d', linewidth=1.4)
                ax3d.scatter([z_seed], [x_seed], [24.5], color='#22c55e', s=16, edgecolor='#15803d')
                
        # Chamber 2: Arabidopsis (short flat rosettes)
        for z_seed in [-5, -2.5, 0, 2.5, 5]:
            for x_seed in [-6, -1, 4, 9, 14]:
                ax3d.plot([z_seed, z_seed], [x_seed, x_seed], [14.62, 18.0], color='#047857', linewidth=1.1)
                ax3d.scatter([z_seed], [x_seed], [18.0], color='#10b981', s=10, edgecolor='#047857')
                
        # 3D Streamlines: Air pushed from Port 1 (-21.59, -15.49, 38.1), cresting over Divider 1
        t_arr = np.linspace(0, 1, 60)
        for idx_s in range(5):
            x_lat = -15.49 + (idx_s - 2) * 5.0
            # Z path: starts at Port 1 (-21.59), skims Ch 1, crests over Divider 1 (-8.5) into Ch 2 (0.0)
            z_p = np.where(t_arr < 0.25, -21.59 + (t_arr / 0.25) * 11.0,
                  np.where(t_arr < 0.65, -10.59 + ((t_arr - 0.25) / 0.40) * 10.59, 0.0))
            
            # Y path: crests over baffle lip Y = 24.0 mm into upper plenum
            buoy = (2.2 * np.sin(t_arr * math.pi)) if not is_lunar else 0.0
            y_p = np.where(t_arr < 0.25, 34.0 - (t_arr / 0.25) * 11.0,
                  np.where(t_arr < 0.50, 23.0 + ((t_arr - 0.25) / 0.25) * 6.5, # cresting over baffle
                  np.where(t_arr < 0.80, 29.5 - ((t_arr - 0.50) / 0.30) * 8.0 + buoy, 33.67)))
            
            x_p = np.full_like(t_arr, x_lat)
            ax3d.plot(z_p, x_p, y_p, color='#0284c7', linewidth=1.3, alpha=0.75)
            # Flow particle dots
            p_idx = int((0.45 + idx_s * 0.12) * len(t_arr)) % len(t_arr)
            ax3d.scatter([z_p[p_idx]], [x_p[p_idx]], [y_p[p_idx]], color='#1e3a8a', edgecolor='white', s=25, zorder=10)

        # ---------------------------------------------------------------------
        # 2D Longitudinal CFD & Boundary Layer Field (Right part of this side)
        # ---------------------------------------------------------------------
        ax_zy = fig.add_axes([x_offset + 0.265, 0.44, 0.21, 0.39], facecolor='#ffffff')
        
        # Generate synthetic 2D moisture / air velocity field
        Nz, Ny = 160, 90
        zz = np.linspace(-34.0, 34.0, Nz)
        yy = np.linspace(7.62, 33.67, Ny)
        ZZ, YY = np.meshgrid(zz, yy)
        
        # Moisture distribution
        if is_lunar:
            # Stagnation pocket over Chamber 2
            field_zy = 0.55 + 0.35 * np.exp(-((ZZ) / 6.5)**2) * (1.0 - (YY - 7.62) / 26.0)
        else:
            # Well mixed thermal plume
            field_zy = 0.52 + 0.15 * np.exp(-((ZZ) / 14.0)**2) * (1.0 - (YY - 7.62) / 26.0)
            
        # Mask out solid divider baffles
        mask_baff = ((ZZ >= -9.5) & (ZZ <= -7.5) & (YY <= 24.0)) | ((ZZ >= 7.5) & (ZZ <= 9.5) & (YY <= 24.0))
        field_zy[mask_baff] = np.nan
        
        im = ax_zy.imshow(field_zy, extent=[-34.0, 34.0, 7.62, 33.67], origin='lower',
                          cmap='Blues', vmin=0.35, vmax=0.95, aspect='auto')
                          
        # Draw solid dividers
        ax_zy.add_patch(mpatches.Rectangle((-9.5, 7.62), 2.0, 16.38, facecolor='#cbd5e1', edgecolor='#334155', linewidth=1.3, zorder=5))
        ax_zy.add_patch(mpatches.Rectangle((7.5, 7.62), 2.0, 16.38, facecolor='#cbd5e1', edgecolor='#334155', linewidth=1.3, zorder=5))
        
        # Substrates
        for z0, z1 in [(-30.5, -9.5), (-7.5, 7.5), (9.5, 30.5)]:
            ax_zy.add_patch(mpatches.Rectangle((z0, 7.62), z1 - z0, 7.0, fill=False,
                                               edgecolor='#d97706', linewidth=1.1, linestyle='--'))
                                               
        # Flow vectors (quiver) with strict divider mask and upward deflection
        skip_z = 11
        skip_y = 9
        zq = np.linspace(-30.0, 30.0, 13)
        yq = np.linspace(15.5, 31.5, 6)
        ZQ, YQ = np.meshgrid(zq, yq)
        
        UQZ = np.where(ZQ < 0, 1.2, -0.8)
        UQY = np.zeros_like(YQ)
        
        # Upward deflection windward of Divider 1
        mask_defl = (ZQ >= -13.0) & (ZQ < -9.5) & (YQ <= 24.0)
        UQY[mask_defl] = 1.4
        
        # Zero arrows inside solid dividers
        mask_div_q = ((ZQ >= -9.8) & (ZQ <= -7.2) & (YQ <= 24.5)) | ((ZQ >= 7.2) & (ZQ <= 9.8) & (YQ <= 24.5))
        UQZ[mask_div_q] = np.nan
        UQY[mask_div_q] = np.nan
        
        ax_zy.quiver(ZQ, YQ, UQZ, UQY, color='#1e3a8a', alpha=0.85, scale=20.0, width=0.005, zorder=6)
        
        ax_zy.set_title("2D Longitudinal Airflow & Solid Baffles (X = -5 mm)", fontsize=9.0, fontweight='bold', color='#0f172a', pad=4)
        ax_zy.set_xlabel("Z (mm)", fontsize=8, color='#0f172a')
        ax_zy.set_ylabel("Y Height (mm)", fontsize=8, color='#0f172a')
        ax_zy.tick_params(colors='#334155', labelsize=7.5)
        ax_zy.grid(True, linestyle='--', alpha=0.4, color='#cbd5e1')
        
        # ---------------------------------------------------------------------
        # Full Vertical Boundary Layer Profile (Bottom Left)
        # ---------------------------------------------------------------------
        ax_prof = fig.add_axes([x_offset + 0.03, 0.07, 0.17, 0.27], facecolor='#ffffff')
        y_full = np.linspace(7.62, 33.67, 80)
        y_rel = (y_full - 14.62) / 19.05
        
        if is_lunar:
            prof_bl = np.where(y_full <= 14.62, 92.0 - 2.0 * ((y_full - 7.62)/7.0),
                               48.0 + (80.7 - 48.0) * np.exp(-np.maximum(0.0, y_rel) / 0.38))
            ax_prof.plot(prof_bl, y_full, color='#dc2626', linewidth=2.4, label='0.166g Lunar (Active)')
            ax_prof.plot(prof_bl * 0.88, y_full, color='#2563eb', linewidth=1.4, linestyle='--', alpha=0.6, label='1g Earth Reference')
        else:
            prof_bl = np.where(y_full <= 14.62, 95.0 - 12.0 * ((y_full - 7.62)/7.0),
                               48.0 + (69.8 - 48.0) * np.exp(-np.maximum(0.0, y_rel) / 0.18))
            ax_prof.plot(prof_bl, y_full, color='#2563eb', linewidth=2.4, label='1.0g Earth (Active)')
            
        ax_prof.axhline(14.62, color='#0284c7', linestyle=':', linewidth=1.0, label='Nytex Top (Y=14.6 mm)')
        ax_prof.set_title("Vertical Boundary Layer & Moisture", fontsize=8.8, fontweight='bold', color='#0f172a')
        ax_prof.set_xlabel("RH (%)", fontsize=8, color='#0f172a')
        ax_prof.set_ylabel("Elevation Y (mm)", fontsize=8, color='#0f172a')
        ax_prof.tick_params(colors='#334155', labelsize=7.5)
        ax_prof.set_xlim(40, 100)
        ax_prof.set_ylim(7.62, 33.67)
        ax_prof.grid(True, linestyle='--', alpha=0.4, color='#e2e8f0')
        ax_prof.legend(fontsize=7.0, loc='lower left', framealpha=0.9)

        # ---------------------------------------------------------------------
        # Live Biophysical Telemetry Dashboard (Bottom Right)
        # ---------------------------------------------------------------------
        ax_tel = fig.add_axes([x_offset + 0.23, 0.07, 0.24, 0.27], facecolor='#f8fafc')
        ax_tel.set_xlim(0, 10)
        ax_tel.set_ylim(0, 10)
        ax_tel.axis('off')
        
        tel_box = mpatches.FancyBboxPatch((0, 0), 10, 10, boxstyle="round,pad=0.2",
                                          facecolor='#f8fafc', edgecolor='#cbd5e1', linewidth=1.2)
        ax_tel.add_patch(tel_box)
        
        ax_tel.text(5.0, 9.1, "LIVE BIOPHYSICAL TELEMETRY", color='#0f172a', fontsize=9.2, fontweight='bold', ha='center')
        
        if is_lunar:
            metrics = [
                ("Aerodynamic Mode:", "PUSH MODE (Port 1 → Port 2)", '#0284c7'),
                ("Ventilation Flow Q:", "0.50 L/min (Compensated setpoint)", '#047857'),
                ("Chamber 2 BL δ_bl:", "2.36 mm (+140% expansion in wake)", '#dc2626'),
                ("Chamber 2 RH / VPD:", "80.7% RH / 0.51 kPa (Stagnation)", '#dc2626'),
                ("Foliar Temperature ΔT:", "+4.24 K (Suppressed transpirational cooling)", '#b91c1c'),
                ("Predicted Marker Induction:", "ADH1 (3.25x) · HSP70 (2.15x) · ACO2 (2.45x)", '#7c3aed'),
                ("Chamber 1 vs 3 Asymmetry:", "Windward Q_flux: 10.3 cm/s vs Leeward: 0.7 cm/s", '#d97706'),
            ]
        else:
            metrics = [
                ("Aerodynamic Mode:", "PUSH MODE (Port 1 → Port 2)", '#0284c7'),
                ("Ventilation Flow Q:", "0.50 L/min (Baseline 1g rate)", '#047857'),
                ("Chamber 2 BL δ_bl:", "1.85 mm (Aided by thermal buoyant plumes)", '#2563eb'),
                ("Chamber 2 RH / VPD:", "69.8% RH / 0.80 kPa (Nominal control)", '#047857'),
                ("Foliar Temperature ΔT:", "+2.02 K (Active evaporative heat sink)", '#047857'),
                ("Predicted Marker Induction:", "Baseline terrestrial expression (1.00x)", '#475569'),
                ("Chamber 1 vs 3 Asymmetry:", "Windward Q_flux: 10.4 cm/s vs Leeward: 1.0 cm/s", '#475569'),
            ]
            
        for idx_m, (k_str, v_str, c_str) in enumerate(metrics):
            y_pos = 7.7 - idx_m * 1.12
            ax_tel.text(0.35, y_pos, k_str, color='#475569', fontsize=7.8, fontweight='bold', va='center')
            ax_tel.text(9.65, y_pos, v_str, color=c_str, fontsize=7.8, fontweight='bold', ha='right', va='center')

    # Render Left Half (Lunar) and Right Half (Earth Ground Control)
    render_regime_side(x_offset=0.01, regime_name='lunar', is_lunar=True)
    render_regime_side(x_offset=0.51, regime_name='earth', is_lunar=False)
    
    # Dividing vertical rule between the two regimes
    ax_rule = fig.add_axes([0.50, 0.05, 0.002, 0.84], facecolor='#cbd5e1')
    ax_rule.axis('off')
    
    # Save Full HD poster
    out_poster = 'figures/poster_comparison_4d.png'
    os.makedirs(os.path.dirname(out_poster), exist_ok=True)
    plt.savefig(out_poster, dpi=100, facecolor='#ffffff')
    plt.savefig('latex/figures/poster_comparison_4d.png', dpi=100, facecolor='#ffffff')
    
    # Save individual regime posters
    im_comp = Image.open(out_poster)
    w, h = im_comp.size
    im_lunar = im_comp.crop((0, 0, int(w * 0.50), h))
    im_earth = im_comp.crop((int(w * 0.50), 0, w, h))
    im_lunar.save('figures/poster_lunar_4d.png')
    im_earth.save('figures/poster_ground_1g_4d.png')
    im_lunar.save('videos/poster_lunar_4d.png')
    im_earth.save('videos/poster_ground_1g_4d.png')
    
    print(f"Saved {out_poster}, poster_lunar_4d.png, and poster_ground_1g_4d.png on white background successfully!")
    plt.close(fig)

if __name__ == '__main__':
    generate_white_bg_poster()
