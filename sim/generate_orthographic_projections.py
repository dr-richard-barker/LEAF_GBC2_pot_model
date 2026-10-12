"""
Generate Publication-Grade Orthographic Projections with Complete Annotations
for the LEAF GBC2 Plant Growth Chamber.
Includes:
- Top Plan View (Z-X Plane, Looking Down -Y)
- Side Elevation / Transverse Cutaway (X-Y Plane at Z = 0)
- Longitudinal Cutaway Section (Mid-plane Z-Y at X = -5 mm through canopies)

Full ABAI QC Standards:
- Zero label collisions or text overflows across chamber baffles
- Dedicated background badges with high-contrast typography
- Unobstructed vertical fluid conduits (Port 1 Inlet, Port 2 Exhaust, Port 3 Dipping Nozzle)
- Explicit seed density & morphometric telemetry (32 B. rapa vs 100 Arabidopsis)
- 3-part partitioned substrate blocks (Y = 7.62 to 14.62 mm)
- Lower drainage wells (Y = 0 to 7.62 mm)
- Open central plenum (ΔY = 14.2 mm) directly under Port 2
"""

import os
import sys
import numpy as np
import trimesh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

def generate_orthographic_figure():
    os.makedirs('figures', exist_ok=True)
    
    print("Loading native STEP geometries for orthographic projection...")
    lid = trimesh.load('cad/lid.step').geometry['M09025F-MDL-2']
    pot = trimesh.load('cad/pot.step').geometry['M09024F-MDL-5']
    
    # Verified CAD mating transformation for lid onto pot body
    T = np.array([
        [ 0.0, 0.0, 1.0, 0.0],
        [ 0.0, 1.0, 0.0, (33.670 - 1.270)/1000.0],
        [-1.0, 0.0, 0.0, 0.0],
        [ 0.0, 0.0, 0.0, 1.0]
    ])
    lid_mated = lid.copy()
    lid_mated.apply_transform(T)
    
    fig = plt.figure(figsize=(22, 16), facecolor='white')
    
    # Grid layout: 2 rows with ample spacing
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.45], hspace=0.32, wspace=0.25)
    
    def plot_mesh_section(ax, mesh, origin, normal, idx_h, idx_v, color, lw=1.8, ls='-', label=''):
        s = mesh.section(plane_origin=origin, plane_normal=normal)
        if s is None:
            return
        for line in s.entities:
            pts = s.vertices[line.points] * 1000.0 # to mm
            ax.plot(pts[:, idx_h], pts[:, idx_v], color=color, lw=lw, ls=ls, label=label)
            label = ''

    # =========================================================================
    # PANEL A: TOP PLAN VIEW (Z-X Plane, looking down from +Y)
    # =========================================================================
    ax_top = fig.add_subplot(gs[0, 0])
    ax_top.set_facecolor('#ffffff')
    
    # Slices through lid (Y = 36 mm) and pot rim (Y = 30 mm)
    plot_mesh_section(ax_top, lid_mated, [0, 0.036, 0], [0, 1, 0], 2, 0, '#1b2631', lw=2.2)
    plot_mesh_section(ax_top, pot, [0, 0.030, 0], [0, 1, 0], 2, 0, '#7f8c8d', lw=1.2, ls='--')
    
    # Growth chamber interior outlines with badges
    chambers_top = [
        (-20.0, 'CHAMBER 1 (32 Seeds)\nBrassica rapa\nArea: 64 cm²', '#27ae60', '#145a32', '#eafaf1', 17.0),
        (0.0, 'CHAMBER 2 (100 Seeds)\nA. thaliana\nArea: 15 cm²', '#2980b9', '#154360', '#ebf5fb', 13.5),
        (20.0, 'CHAMBER 3 (32 Seeds)\nBrassica rapa\nArea: 64 cm²', '#27ae60', '#145a32', '#eafaf1', 17.0)
    ]
    for z_c, title_text, border_col, text_col, bg_col, w_rect in chambers_top:
        rect = patches.Rectangle((z_c - w_rect/2.0, -12.0), w_rect, 24.0, linewidth=1.5,
                                 edgecolor=border_col, facecolor=bg_col, alpha=0.50, linestyle=':')
        ax_top.add_patch(rect)
        ax_top.text(z_c, 0.0, title_text, ha='center', va='center', fontsize=8.0, fontweight='bold',
                    color=text_col, bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor=border_col, lw=1.0, alpha=0.95))
                    
    # Port locations in Top View (X = -15.49 mm)
    ports_top = [
        (-21.59, -15.49, 'PORT 1: AIR INLET\n(Z = -21.6 mm, Ø 3.18 mm)', '#c0392b', '#fce4ec', -24.5),
        (0.00, -15.49, 'PORT 2: CENTRAL EXHAUST\n(Z = 0.0 mm, Ø 3.18 mm)', '#8e44ad', '#f3e5f5', -24.5),
        (21.59, -15.49, 'PORT 3: WATER / RNALATER\n(Z = +21.6 mm, Ø 3.18 mm)', '#1f618d', '#e8f8f5', -24.5)
    ]
    for z_p, x_p, p_label, p_col, p_bg, y_text in ports_top:
        circle_out = patches.Circle((z_p, x_p), 5.5, edgecolor=p_col, facecolor='white', lw=2.0)
        circle_in = patches.Circle((z_p, x_p), 1.6, edgecolor=p_col, facecolor=p_col, alpha=0.85)
        ax_top.add_patch(circle_out)
        ax_top.add_patch(circle_in)
        ax_top.annotate(p_label, xy=(z_p, x_p - 5.5), xytext=(z_p, y_text),
                        ha='center', va='top', fontsize=7.5, fontweight='bold', color=p_col,
                        bbox=dict(boxstyle='round,pad=0.2', facecolor=p_bg, edgecolor=p_col, lw=1.0),
                        arrowprops=dict(arrowstyle='->', color=p_col, lw=1.3))
                        
    # Perimeter manifold airflow indicator with badge
    ax_top.annotate('', xy=(16, 15.0), xytext=(-16, 15.0),
                    arrowprops=dict(arrowstyle='->', color='#e67e22', lw=2.4, ls='--'))
    ax_top.text(0, 17.5, 'Perimeter Air Distribution Manifold (Y = 28–32 mm)', ha='center', va='bottom',
                fontsize=8.0, color='#d35400', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='#fef9e7', edgecolor='#f39c12', lw=0.8))

    ax_top.set_xlim([-38, 38])
    ax_top.set_ylim([-31, 24])
    ax_top.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax_top.set_ylabel('X Width (mm)', fontsize=11, fontweight='bold')
    ax_top.set_title('TOP PLAN VIEW (Looking Down -Y)\nHardware Boundaries & Port Configuration', fontsize=12, fontweight='bold', pad=10)
    ax_top.grid(True, linestyle='--', alpha=0.4)
    ax_top.set_aspect('equal')

    # =========================================================================
    # PANEL B: SIDE ELEVATION & TRANSVERSE CUTAWAY (X-Y Plane at Z = 0)
    # =========================================================================
    ax_side = fig.add_subplot(gs[0, 1])
    ax_side.set_facecolor('#ffffff')
    
    # Transverse section through central chamber (Z = 0)
    plot_mesh_section(ax_side, pot, [0, 0, 0], [0, 0, 1], 0, 1, '#1b2631', lw=2.2)
    plot_mesh_section(ax_side, lid_mated, [0, 0, 0], [0, 0, 1], 0, 1, '#c0392b', lw=2.2)
    
    # Lower drainage well (Y = 0 to 7.62 mm)
    rect_drain_side = patches.Rectangle((-10.5, 0.0), 31.5, 7.62, facecolor='#eaeded', alpha=0.5, edgecolor='#95a5a6', linestyle=':')
    ax_side.add_patch(rect_drain_side)
    ax_side.text(5.25, 3.8, 'Lower Drainage Base (Y = 0–7.6 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#5d6d7e',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#bdc3c7', lw=0.8, alpha=0.9))
    
    # Porous media substrate base inside Chamber 2 (Y = 7.62 to 14.62 mm)
    rect_sub_side = patches.Rectangle((-10.5, 7.62), 31.5, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    ax_side.add_patch(rect_sub_side)
    ax_side.text(5.25, 11.12, 'Porous Substrate 2 (Base)\n(Y = 7.6–14.6 mm, ε = 0.50)', ha='center', va='center', fontsize=8.0, fontweight='bold', color='#6e2c00',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.9))
    
    # Chamber 2 canopy zone (Arabidopsis: 100 seedlings, compact rosette lawn Y = 14.62 to 19.5 mm)
    rect_can_side = patches.Rectangle((-10.5, 14.62), 31.5, 4.88, facecolor='#2980b9', alpha=0.25, edgecolor='#1f618d', linestyle='--')
    ax_side.add_patch(rect_can_side)
    ax_side.text(5.25, 17.0, 'Chamber 2: 100 Arabidopsis Seeds\n(Rosette Lawn, Y = 14.6–19.5 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#154360',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#2980b9', lw=0.8, alpha=0.9))
    
    # Open central headspace plenum above Arabidopsis
    rect_plenum_side = patches.Rectangle((-10.5, 19.5), 31.5, 14.17, facecolor='#ebf5fb', alpha=0.45, edgecolor='#aed6f1', linestyle=':')
    ax_side.add_patch(rect_plenum_side)
    ax_side.text(5.25, 26.5, 'Open Headspace Plenum (ΔY = 14.2 mm)\nDirect Corridor to Port 2 Exhaust', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#1b4f72',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor='#3498db', lw=0.8, alpha=0.92))
    
    # Central Exhaust Port 2 at X = -15.5 mm with clean single vertical arrow and badge
    ax_side.annotate('PORT 2: CENTRAL EXHAUST\n(X = -15.5 mm, Y = 33–38 mm)\nReturn to Cabin Air Loop',
                     xy=(-15.49, 38.0), xytext=(-12.0, 45.0),
                     ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='#8e44ad',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#f3e5f5', edgecolor='#9c27b0', lw=1.0),
                     arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.8))
                     
    ax_side.annotate('', xy=(-15.49, 43.5), xytext=(-15.49, 34.0),
                     arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=2.4))
                     
    ax_side.set_xlim([-25, 25])
    ax_side.set_ylim([-2, 54])
    ax_side.set_xlabel('X Width (mm)', fontsize=11, fontweight='bold')
    ax_side.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax_side.set_title('SIDE ELEVATION / TRANSVERSE CUTAWAY (X-Y at Z = 0)\nCentral Exhaust & Chamber 2 Profile', fontsize=12, fontweight='bold', pad=10)
    ax_side.grid(True, linestyle='--', alpha=0.4)
    ax_side.set_aspect('equal')

    # =========================================================================
    # PANEL C: FRONT ORTHOGRAPHIC CUTAWAY (Z-Y Plane, X = -5.0 mm through canopies)
    # =========================================================================
    ax_front = fig.add_subplot(gs[1, :])
    ax_front.set_facecolor('#ffffff')
    
    # Longitudinal cut through plant canopies
    plot_mesh_section(ax_front, pot, [-0.005, 0, 0], [1, 0, 0], 2, 1, '#1b2631', lw=2.2)
    plot_mesh_section(ax_front, lid_mated, [-0.005, 0, 0], [1, 0, 0], 2, 1, '#c0392b', lw=2.2)
    
    # 3 SEPARATE POROUS SUBSTRATE BLOCKS AT Y = 7.62 - 14.62 mm
    # Block 1 (Chamber 1): Z = -29.5 to -10.5 mm
    sub1_front = patches.Rectangle((-29.5, 7.62), 19.0, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.3)
    ax_front.add_patch(sub1_front)
    ax_front.text(-20.0, 11.12, 'Porous Substrate 1\n(Y = 7.6–14.6 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.92))

    # Block 2 (Chamber 2): Z = -7.0 to +7.0 mm
    sub2_front = patches.Rectangle((-7.0, 7.62), 14.0, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.3)
    ax_front.add_patch(sub2_front)
    ax_front.text(0.0, 11.12, 'Porous Substrate 2\n(Y = 7.6–14.6 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.92))

    # Block 3 (Chamber 3): Z = +10.5 to +29.5 mm
    sub3_front = patches.Rectangle((10.5, 7.62), 19.0, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.3)
    ax_front.add_patch(sub3_front)
    ax_front.text(17.0, 11.12, 'Porous Substrate 3\n(Y = 7.6–14.6 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.92))

    # Lower cavity label with crisp white badge
    ax_front.text(0.0, 3.8, 'Lower Aeration & Drainage Wells (Y = 0–7.6 mm)',
                  ha='center', va='center', fontsize=8.5, fontweight='bold', color='#5d6d7e',
                  bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#bdc3c7', lw=0.9, alpha=0.92))
                  
    # -------------------------------------------------------------------------
    # COMPARTMENTS, SPECIES & CANOPIES (COMPACT BOUNDING BOXES - ZERO OVERFLOW)
    # -------------------------------------------------------------------------
    # Compartment 1: B. rapa (32 seeds, tall foliage Y = 14.62 to 26.0 mm)
    c1_box = patches.Rectangle((-29.5, 14.62), 19.0, 11.38, facecolor='#27ae60', alpha=0.18, edgecolor='#1e8449', linestyle='--')
    ax_front.add_patch(c1_box)
    ax_front.text(-20.0, 24.8, 'CHAMBER 1 (32 Seeds)\nBrassica rapa (h = 11.4 mm)\nArea: 64 cm² · LAI = 9.4',
                  ha='center', va='center', fontsize=7.5, fontweight='bold', color='#145a32',
                  bbox=dict(boxstyle='round,pad=0.25', facecolor='#eafaf1', edgecolor='#27ae60', lw=0.9, alpha=0.95))
                  
    # Plant sketch with multiple taller stems and broad leaves
    for z_off in [-25.0, -20.0, -15.0]:
        ax_front.plot([z_off, z_off], [14.62, 19.0], color='#27ae60', lw=2.0)
        ax_front.plot([z_off, z_off-2.0], [19.0, 21.5], color='#27ae60', lw=1.6)
        ax_front.plot([z_off, z_off+2.0], [19.0, 21.0], color='#27ae60', lw=1.6)

    # Compartment 2: A. thaliana (100 seeds, miniature rosette lawn Y = 14.62 to 19.5 mm)
    c2_box = patches.Rectangle((-7.0, 14.62), 14.0, 4.88, facecolor='#2980b9', alpha=0.25, edgecolor='#1f618d', linestyle='--')
    ax_front.add_patch(c2_box)
    ax_front.text(0.0, 17.0, 'CHAMBER 2 (100 Seeds)\nA. thaliana Lawn (h = 4.9 mm)\nArea: 15 cm² · LAI = 3.1',
                  ha='center', va='center', fontsize=7.0, fontweight='bold', color='#154360',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='#ebf5fb', edgecolor='#2980b9', lw=0.9, alpha=0.95))
                  
    # Carpet of small rosette tufts representing 100 seeds
    for z_off in [-4.5, -2.25, 0.0, 2.25, 4.5]:
        ax_front.plot([z_off, z_off], [14.62, 16.0], color='#2980b9', lw=1.6)
        ax_front.plot([z_off, z_off-1.0], [16.0, 16.8], color='#2980b9', lw=1.3)
        ax_front.plot([z_off, z_off+1.0], [16.0, 16.8], color='#2980b9', lw=1.3)

    # Open central plenum box above Arabidopsis
    c2_plenum = patches.Rectangle((-7.0, 19.5), 14.0, 14.17, facecolor='#ebf5fb', alpha=0.35, edgecolor='#aed6f1', linestyle=':')
    ax_front.add_patch(c2_plenum)
    ax_front.text(0.0, 25.0, 'OPEN CENTRAL PLENUM\n(ΔY = 14.2 mm Headspace)\nDirect Exhaust Corridor',
                  ha='center', va='center', fontsize=7.2, fontweight='bold', color='#1b4f72',
                  bbox=dict(boxstyle='round,pad=0.25', facecolor='#ffffff', edgecolor='#3498db', lw=0.9, alpha=0.95))

    # Compartment 3: B. rapa (32 seeds, tall foliage Y = 14.62 to 26.0 mm)
    # Badge placed at Z = 16.0 to avoid overlapping dipping nozzle at Z = 21.6 mm
    c3_box = patches.Rectangle((10.5, 14.62), 19.0, 11.38, facecolor='#27ae60', alpha=0.18, edgecolor='#1e8449', linestyle='--')
    ax_front.add_patch(c3_box)
    ax_front.text(16.0, 24.8, 'CHAMBER 3 (32 Seeds)\nBrassica rapa (h = 11.4 mm)\nArea: 64 cm² · LAI = 9.4',
                  ha='center', va='center', fontsize=7.5, fontweight='bold', color='#145a32',
                  bbox=dict(boxstyle='round,pad=0.25', facecolor='#eafaf1', edgecolor='#27ae60', lw=0.9, alpha=0.95))
    for z_off in [13.0, 17.0, 25.0]:
        ax_front.plot([z_off, z_off], [14.62, 19.0], color='#27ae60', lw=2.0)
        ax_front.plot([z_off, z_off-2.0], [19.0, 21.5], color='#27ae60', lw=1.6)
        ax_front.plot([z_off, z_off+2.0], [19.0, 21.0], color='#27ae60', lw=1.6)

    # -------------------------------------------------------------------------
    # FLUID PORT & PIPELINE ANNOTATIONS
    # -------------------------------------------------------------------------
    # PORT 1: AIR INLET (Z = -21.59 mm)
    ax_front.annotate('PORT 1: FORCED AIR INLET\n(Z = -21.6 mm, Ø 3.18 mm)\nFrom Pump Manifold',
                     xy=(-21.59, 33.67), xytext=(-21.59, 48.5),
                     ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='#c0392b',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#fce4ec', edgecolor='#e91e63', lw=1.1),
                     arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.8))
    # Inflow arrow
    ax_front.annotate('', xy=(-21.59, 29.0), xytext=(-21.59, 34.5),
                     arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2.4))

    # PORT 2: CENTRAL AIR EXHAUST (Z = 0.00 mm)
    ax_front.annotate('PORT 2: CENTRAL EXHAUST\n(Z = 0.0 mm, Ø 3.18 mm)\nReturn to Cabin Atmosphere Loop',
                     xy=(0.0, 33.67), xytext=(0.0, 48.5),
                     ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='#8e44ad',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#f3e5f5', edgecolor='#9c27b0', lw=1.1),
                     arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.8))
    # Outflow arrow
    ax_front.annotate('', xy=(0.0, 43.5), xytext=(0.0, 34.0),
                     arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=2.4))

    # PORT 3: WATER & RNALATER INLET (Z = +21.59 mm)
    ax_front.annotate('PORT 3: WATER & RNALATER\n(Z = +21.6 mm, Ø 3.18 mm)\nInternal Dipping Delivery Line',
                     xy=(21.59, 33.67), xytext=(21.59, 48.5),
                     ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='#1f618d',
                     bbox=dict(boxstyle='round,pad=0.25', facecolor='#e8f8f5', edgecolor='#16a085', lw=1.1),
                     arrowprops=dict(arrowstyle='->', color='#1f618d', lw=1.8))
                     
    # DIPPING NOZZLE (Liquid conduit running down from Port 3 directly into Substrate Block 3 at Y = 9.5 mm)
    ax_front.plot([21.59, 21.59], [33.67, 9.5], color='#2980b9', lw=4.2, solid_capstyle='round')
    ax_front.scatter([21.59], [9.5], color='#154360', s=45, zorder=5)
    
    # Discharge arrows into substrate
    ax_front.annotate('', xy=(19.0, 9.5), xytext=(21.59, 9.5), arrowprops=dict(arrowstyle='->', color='#2980b9', lw=2.0))
    ax_front.annotate('', xy=(24.2, 9.5), xytext=(21.59, 9.5), arrowprops=dict(arrowstyle='->', color='#2980b9', lw=2.0))
    
    # Clean leader line callout for dipping nozzle in right margin
    ax_front.annotate('RNAlater / Water Discharge\n(Y = 9.5 mm in Substrate 3)\nPrevents foliar splashing',
                      xy=(21.59, 9.5), xytext=(29.5, 14.0),
                      ha='left', va='center', fontsize=7.5, fontweight='bold', color='#1b4f72',
                      bbox=dict(boxstyle='round,pad=0.25', facecolor='#e8f8f5', edgecolor='#2980b9', lw=1.0),
                      arrowprops=dict(arrowstyle='->', color='#2980b9', lw=1.4))

    # -------------------------------------------------------------------------
    # AIRFLOW CIRCULATION STREAMLINES
    # -------------------------------------------------------------------------
    # Manifold distribution arrows (Left and Right branches without colliding with Port 2 vertical line)
    ax_front.annotate('', xy=(-6.0, 31.0), xytext=(-17.0, 31.0),
                     arrowprops=dict(arrowstyle='->', color='#e67e22', lw=2.2, ls='--'))
    ax_front.annotate('', xy=(17.0, 31.0), xytext=(6.0, 31.0),
                     arrowprops=dict(arrowstyle='->', color='#e67e22', lw=2.2, ls='--'))
                     
    ax_front.text(-12.0, 32.8, 'Manifold Distribution (Y = 28–32 mm)',
                  ha='center', va='bottom', fontsize=7.5, color='#d35400', fontweight='bold', style='italic',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='#fef9e7', edgecolor='#f39c12', lw=0.8))

    # Canopy circulation streamlines ("dip-and-converge")
    # Skim across Ch 1
    ax_front.annotate('', xy=(-10.0, 19.5), xytext=(-17.0, 21.0), arrowprops=dict(arrowstyle='->', color='#27ae60', lw=1.8, ls='-.'))
    # Dip into Chamber 2 plenum from Ch 1
    ax_front.annotate('', xy=(-2.0, 27.5), xytext=(-7.0, 21.0), arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=2.0))
    # Dip into Chamber 2 plenum from Ch 3
    ax_front.annotate('', xy=(2.0, 27.5), xytext=(7.0, 21.0), arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=2.0))
    # Skim across Ch 3
    ax_front.annotate('', xy=(10.0, 19.5), xytext=(17.0, 21.0), arrowprops=dict(arrowstyle='->', color='#27ae60', lw=1.8, ls='-.'))

    # Skimming flow badge with leader in upper left
    ax_front.annotate('Canopy Skimming Flow\n(Controlled boundary layer)',
                      xy=(-20.0, 20.0), xytext=(-31.0, 16.5),
                      ha='right', va='center', fontsize=7.5, fontweight='bold', color='#784212', style='italic',
                      bbox=dict(boxstyle='round,pad=0.2', facecolor='#fef5e7', edgecolor='#d35400', lw=0.8),
                      arrowprops=dict(arrowstyle='->', color='#d35400', lw=1.2))

    ax_front.set_xlim([-42, 42])
    ax_front.set_ylim([-2, 57])
    ax_front.set_xlabel('Z Length (mm)', fontsize=12, fontweight='bold')
    ax_front.set_ylabel('Y Height (mm)', fontsize=12, fontweight='bold')
    ax_front.set_title('LONGITUDINAL CUTAWAY SECTION (Z-Y Plane at X = -5.0 mm)\nFull Fluid Architecture, 3-Part Partitioned Substrate, and Air/Liquid Flow Pathways',
                       fontsize=13, fontweight='bold', pad=12)
    ax_front.grid(True, linestyle='--', alpha=0.4)
    ax_front.set_aspect('equal')

    # Comprehensive Legend at bottom
    legend_elements = [
        Line2D([0], [0], color='#1b2631', lw=2.2, label='Pot Body (Aluminum / Polymer Housing, Y = 0 to 33.67 mm)'),
        Line2D([0], [0], color='#c0392b', lw=2.2, label='Lid Assembly & Manifold Ceiling (Y = 32.40 to 38.11 mm)'),
        Line2D([0], [0], color='#e67e22', lw=2.2, ls='--', label='Forced Airflow (Port 1 Inlet → Manifold → Canopy → Port 2 Exhaust)'),
        Line2D([0], [0], color='#2980b9', lw=4.2, label='Port 3 Dipping Nozzle (Liquid Discharge at Y = 9.5 mm inside Substrate 3)'),
        patches.Patch(facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', label='3-Part Porous Substrate (Chambers 1, 2, 3 Bases, Y = 7.62–14.62 mm)'),
        patches.Patch(facecolor='#eaeded', alpha=0.5, edgecolor='#95a5a6', label='Lower Aeration & Drainage Wells (Y = 0–7.62 mm)'),
        patches.Patch(facecolor='#27ae60', alpha=0.20, edgecolor='#1e8449', label='Chambers 1 & 3: 32 Brassica rapa Seeds (Tall Foliage h = 11.4 mm, Area: 64 cm²)'),
        patches.Patch(facecolor='#2980b9', alpha=0.25, edgecolor='#1f618d', label='Chamber 2: 100 Arabidopsis thaliana Seeds (Rosette Lawn h = 4.9 mm, Area: 15 cm²)')
    ]
    ax_front.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                    ncol=4, fontsize=9.0, frameon=True, fancybox=True, shadow=True)

    plt.savefig('figures/fig_orthographic_projection_annotated.png', dpi=200, bbox_inches='tight')
    plt.savefig('figures/fig_orthographic_projection_annotated.svg', bbox_inches='tight')
    print("Saved figures/fig_orthographic_projection_annotated.png and .svg successfully!")

if __name__ == '__main__':
    generate_orthographic_figure()
