"""
Generate CAD Mated Assembly Cross-Sections for LEAF GBC2 Plant Growth Chamber.
Extracts true 3D cross-section contours directly from native STEP files
without unaligned 2D local projection distortions.

Includes 3 orthogonal / cutaway sections:
1. Slice at X = -5.0 mm (Longitudinal through the 3 plant growth compartments)
2. Slice at X = -15.49 mm (Longitudinal passing directly through the 3 fluid ports)
3. Transverse Section at Z = 0.0 mm (Across central Chamber 2 and Port 2 Exhaust)

Full ABAI QC Standards:
- Zero label overlap, high contrast badges, clear spacing
- Exact STEP part mating: Pot (Y = 0.0 to 33.67 mm), Lid (Y = 32.40 to 38.11 mm)
- 3-part partitioned substrate blocks (Y = 7.62 to 14.62 mm)
- Lower drainage wells (Y = 0.0 to 7.62 mm)
- 32 B. rapa (tall, h = 11.4 mm) vs 100 Arabidopsis (lawn, h = 4.9 mm)
- Port 3 dipping nozzle down to Y = 9.5 mm inside Substrate Block 3
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

def generate_assembly_cross_sections():
    os.makedirs('figures', exist_ok=True)
    
    print("Loading native STEP geometries for assembly cross sections...")
    lid = trimesh.load('cad/lid.step').geometry['M09025F-MDL-2']
    pot = trimesh.load('cad/pot.step').geometry['M09024F-MDL-5']
    
    # Precise CAD mating transformation for lid onto pot body
    T = np.array([
        [ 0.0, 0.0, 1.0, 0.0],
        [ 0.0, 1.0, 0.0, (33.670 - 1.270)/1000.0],
        [-1.0, 0.0, 0.0, 0.0],
        [ 0.0, 0.0, 0.0, 1.0]
    ])
    lid_mated = lid.copy()
    lid_mated.apply_transform(T)
    
    # =========================================================================
    # FIGURE 1: ANNOTATED ASSEMBLY CROSS SECTIONS
    # =========================================================================
    fig, axes = plt.subplots(1, 3, figsize=(22, 7.8), facecolor='white')
    plt.subplots_adjust(left=0.05, right=0.97, bottom=0.16, top=0.88, wspace=0.22)
    
    def plot_mesh_section(ax, mesh, origin, normal, idx_h, idx_v, color, lw=2.0, ls='-', label=''):
        s = mesh.section(plane_origin=origin, plane_normal=normal)
        if s is None:
            return
        for line in s.entities:
            pts = s.vertices[line.points] * 1000.0 # to mm
            ax.plot(pts[:, idx_h], pts[:, idx_v], color=color, lw=lw, ls=ls, label=label)
            label = ''

    # -------------------------------------------------------------------------
    # PANEL 1: LONGITUDINAL SECTION AT X = -5.0 mm (GROWTH COMPARTMENTS)
    # -------------------------------------------------------------------------
    ax1 = axes[0]
    ax1.set_facecolor('#ffffff')
    
    plot_mesh_section(ax1, pot, [-0.005, 0, 0], [1, 0, 0], 2, 1, '#1b2631', lw=2.2)
    plot_mesh_section(ax1, lid_mated, [-0.005, 0, 0], [1, 0, 0], 2, 1, '#c0392b', lw=2.2)
    
    # 3-Part Partitioned Substrate (Y = 7.62 to 14.62 mm)
    sub1 = patches.Rectangle((-29.5, 7.62), 19.0, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    sub2 = patches.Rectangle((-7.0, 7.62), 14.0, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    sub3 = patches.Rectangle((10.5, 7.62), 19.0, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    ax1.add_patch(sub1)
    ax1.add_patch(sub2)
    ax1.add_patch(sub3)
    
    # Lower drainage wells (Y = 0 to 7.62 mm)
    drain = patches.Rectangle((-29.5, 0.0), 59.0, 7.62, facecolor='#eaeded', alpha=0.45, edgecolor='#95a5a6', linestyle=':')
    ax1.add_patch(drain)
    
    # Foliage canopies
    c1_canopy = patches.Rectangle((-29.5, 14.62), 19.0, 11.38, facecolor='#27ae60', alpha=0.18, edgecolor='#1e8449', linestyle='--')
    c2_canopy = patches.Rectangle((-7.0, 14.62), 14.0, 4.88, facecolor='#2980b9', alpha=0.25, edgecolor='#1f618d', linestyle='--')
    c2_plenum = patches.Rectangle((-7.0, 19.5), 14.0, 14.17, facecolor='#ebf5fb', alpha=0.35, edgecolor='#aed6f1', linestyle=':')
    c3_canopy = patches.Rectangle((10.5, 14.62), 19.0, 11.38, facecolor='#27ae60', alpha=0.18, edgecolor='#1e8449', linestyle='--')
    ax1.add_patch(c1_canopy)
    ax1.add_patch(c2_canopy)
    ax1.add_patch(c2_plenum)
    ax1.add_patch(c3_canopy)
    
    # Crisp descriptive labels with badges
    ax1.text(-20.0, 20.5, 'Chamber 1:\n32 B. rapa Seeds\n(h = 11.4 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#145a32',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='#eafaf1', edgecolor='#27ae60', lw=0.8, alpha=0.92))
    ax1.text(0.0, 17.0, 'Chamber 2:\n100 Arabidopsis\n(h = 4.9 mm)', ha='center', va='center', fontsize=7.0, fontweight='bold', color='#154360',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='#ebf5fb', edgecolor='#2980b9', lw=0.8, alpha=0.92))
    ax1.text(0.0, 26.5, 'Open Headspace\n(ΔY = 14.2 mm Plenum)', ha='center', va='center', fontsize=7.2, fontweight='bold', color='#1b4f72',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor='#3498db', lw=0.8, alpha=0.92))
    ax1.text(20.0, 20.5, 'Chamber 3:\n32 B. rapa Seeds\n(h = 11.4 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#145a32',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='#eafaf1', edgecolor='#27ae60', lw=0.8, alpha=0.92))
    
    ax1.text(-20.0, 11.12, 'Substrate 1', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.9))
    ax1.text(0.0, 11.12, 'Substrate 2', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.9))
    ax1.text(20.0, 11.12, 'Substrate 3', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.9))
             
    ax1.text(0.0, 3.8, 'Lower Aeration & Drainage Wells (Y = 0–7.6 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#5d6d7e',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#bdc3c7', lw=0.8, alpha=0.9))
    
    ax1.set_xlim([-37, 37])
    ax1.set_ylim([-2, 46])
    ax1.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax1.set_title('Cross-Section at X = -5.0 mm\n(Through 3 Plant Growth Compartments)', fontsize=12, fontweight='bold', pad=10)
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.set_aspect('equal')

    # -------------------------------------------------------------------------
    # PANEL 2: LONGITUDINAL SECTION AT X = -15.49 mm (3 FLUID PORTS)
    # -------------------------------------------------------------------------
    ax2 = axes[1]
    ax2.set_facecolor('#ffffff')
    
    plot_mesh_section(ax2, pot, [-0.015494, 0, 0], [1, 0, 0], 2, 1, '#1b2631', lw=2.2)
    plot_mesh_section(ax2, lid_mated, [-0.015494, 0, 0], [1, 0, 0], 2, 1, '#c0392b', lw=2.2)
    
    # 3-part substrate profile at port X plane
    sub1_p2 = patches.Rectangle((-29.5, 7.62), 19.0, 7.0, facecolor='#e59866', alpha=0.35, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    sub2_p2 = patches.Rectangle((-7.0, 7.62), 14.0, 7.0, facecolor='#e59866', alpha=0.35, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    sub3_p2 = patches.Rectangle((10.5, 7.62), 19.0, 7.0, facecolor='#e59866', alpha=0.35, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    ax2.add_patch(sub1_p2)
    ax2.add_patch(sub2_p2)
    ax2.add_patch(sub3_p2)
    
    # Port 3 Dipping Nozzle tube extending down to Y = 9.5 mm inside Substrate 3
    ax2.plot([21.59, 21.59], [33.67, 9.5], color='#2980b9', lw=4.2, solid_capstyle='round')
    ax2.scatter([21.59], [9.5], color='#154360', s=45, zorder=5)
    ax2.annotate('Nozzle Discharge\n(Y = 9.5 mm)', xy=(21.59, 9.5), xytext=(27.5, 4.0),
                 ha='center', va='top', fontsize=7.5, fontweight='bold', color='#1f618d',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='#e8f8f5', edgecolor='#1f618d', lw=0.8),
                 arrowprops=dict(arrowstyle='->', color='#1f618d', lw=1.3))
                 
    # Port annotations with clear horizontal spacing
    ax2.annotate('PORT 1: AIR INLET\n(Z = -21.6 mm, Ø 3.18 mm)', xy=(-21.59, 38.1), xytext=(-21.59, 44.5),
                 ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='#c0392b',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#fce4ec', edgecolor='#e91e63', lw=1.0),
                 arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5))
                 
    ax2.annotate('PORT 2: CENTRAL EXHAUST\n(Z = 0.0 mm, Ø 3.18 mm)', xy=(0.0, 38.1), xytext=(0.0, 44.5),
                 ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='#8e44ad',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#f3e5f5', edgecolor='#9c27b0', lw=1.0),
                 arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.5))
                 
    ax2.annotate('PORT 3: WATER / RNALATER\n(Z = +21.6 mm, Ø 3.18 mm)', xy=(21.59, 38.1), xytext=(21.59, 44.5),
                 ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='#1f618d',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#e8f8f5', edgecolor='#16a085', lw=1.0),
                 arrowprops=dict(arrowstyle='->', color='#1f618d', lw=1.5))
                 
    # Upper perimeter manifold channelling indication
    ax2.annotate('', xy=(15.0, 30.5), xytext=(-15.0, 30.5),
                 arrowprops=dict(arrowstyle='->', color='#e67e22', lw=2.2, ls='--'))
    ax2.text(0.0, 31.8, 'Perimeter Air Distribution Manifold', ha='center', va='bottom',
             fontsize=7.5, color='#d35400', fontweight='bold', style='italic',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='#fef9e7', edgecolor='#f39c12', lw=0.8))

    ax2.set_xlim([-37, 37])
    ax2.set_ylim([-2, 53])
    ax2.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax2.set_title('Cross-Section at X = -15.49 mm\n(Passing Directly Through the 3 Fluid Ports)', fontsize=12, fontweight='bold', pad=10)
    ax2.grid(True, linestyle='--', alpha=0.4)
    ax2.set_aspect('equal')

    # -------------------------------------------------------------------------
    # PANEL 3: TRANSVERSE SECTION AT Z = 0.0 mm (ACROSS CENTER COMPARTMENT)
    # -------------------------------------------------------------------------
    ax3 = axes[2]
    ax3.set_facecolor('#ffffff')
    
    plot_mesh_section(ax3, pot, [0, 0, 0], [0, 0, 1], 0, 1, '#1b2631', lw=2.2)
    plot_mesh_section(ax3, lid_mated, [0, 0, 0], [0, 0, 1], 0, 1, '#c0392b', lw=2.2)
    
    # Substrate Block 2 inside Chamber 2 (Y = 7.62 to 14.62 mm)
    sub2_trans = patches.Rectangle((-10.5, 7.62), 31.5, 7.0, facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', lw=1.2)
    ax3.add_patch(sub2_trans)
    
    # Arabidopsis canopy (Y = 14.62 to 19.5 mm)
    can2_trans = patches.Rectangle((-10.5, 14.62), 31.5, 4.88, facecolor='#2980b9', alpha=0.25, edgecolor='#1f618d', linestyle='--')
    ax3.add_patch(can2_trans)
    
    # Open plenum (Y = 19.5 to 33.67 mm)
    plen_trans = patches.Rectangle((-10.5, 19.5), 31.5, 14.17, facecolor='#ebf5fb', alpha=0.35, edgecolor='#aed6f1', linestyle=':')
    ax3.add_patch(plen_trans)
    
    # Lower drainage well
    drain_trans = patches.Rectangle((-10.5, 0.0), 31.5, 7.62, facecolor='#eaeded', alpha=0.45, edgecolor='#95a5a6', linestyle=':')
    ax3.add_patch(drain_trans)
    
    # Badges
    ax3.text(5.25, 11.12, 'Porous Substrate 2 (Base)\n(Y = 7.6–14.6 mm, ε = 0.50)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#6e2c00',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ba4a00', lw=0.8, alpha=0.9))
    ax3.text(5.25, 17.0, 'Chamber 2: 100 Arabidopsis Seeds\n(Rosette Lawn, Y = 14.6–19.5 mm)', ha='center', va='center', fontsize=7.2, fontweight='bold', color='#154360',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#2980b9', lw=0.8, alpha=0.9))
    ax3.text(5.25, 26.5, 'Open Headspace Plenum (ΔY = 14.2 mm)\nDirect Corridor to Port 2 Exhaust', ha='center', va='center', fontsize=7.2, fontweight='bold', color='#1b4f72',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='#ffffff', edgecolor='#3498db', lw=0.8, alpha=0.92))
    ax3.text(5.25, 3.8, 'Lower Drainage Base (Y = 0–7.6 mm)', ha='center', va='center', fontsize=7.5, fontweight='bold', color='#5d6d7e',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#bdc3c7', lw=0.8, alpha=0.9))
    
    # Port 2 exhaust stack annotation
    ax3.annotate('PORT 2: CENTRAL EXHAUST\n(X = -15.5 mm, Y = 33–38 mm)', xy=(-15.49, 38.1), xytext=(-12.0, 44.5),
                 ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='#8e44ad',
                 bbox=dict(boxstyle='round,pad=0.25', facecolor='#f3e5f5', edgecolor='#9c27b0', lw=1.0),
                 arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.5))
                 
    # Exit flow arrow
    ax3.annotate('', xy=(-15.49, 43.5), xytext=(-15.49, 34.0),
                 arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=2.2))

    ax3.set_xlim([-25, 25])
    ax3.set_ylim([-2, 53])
    ax3.set_xlabel('X Width (mm)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax3.set_title('Transverse Section at Z = 0.0 mm\n(Across Center Compartment & Port 2)', fontsize=12, fontweight='bold', pad=10)
    ax3.grid(True, linestyle='--', alpha=0.4)
    ax3.set_aspect('equal')

    # Overall Legend
    legend_elements = [
        Line2D([0], [0], color='#1b2631', lw=2.2, label='Pot Body (STEP M09024F-MDL-5, Y = 0 to 33.67 mm)'),
        Line2D([0], [0], color='#c0392b', lw=2.2, label='Lid Assembly (STEP M09025F-MDL-2, Y = 32.40 to 38.11 mm)'),
        patches.Patch(facecolor='#e59866', alpha=0.45, hatch='xxx', edgecolor='#ba4a00', label='3-Part Porous Substrate (Chambers 1, 2, 3 Bases, Y = 7.62–14.62 mm)'),
        patches.Patch(facecolor='#eaeded', alpha=0.45, edgecolor='#95a5a6', label='Lower Aeration & Drainage Wells (Y = 0–7.62 mm)'),
        patches.Patch(facecolor='#27ae60', alpha=0.18, edgecolor='#1e8449', label='Chambers 1 & 3: 32 B. rapa Seeds (Tall Foliage h = 11.4 mm, Area: 64 cm²)'),
        patches.Patch(facecolor='#2980b9', alpha=0.25, edgecolor='#1f618d', label='Chamber 2: 100 Arabidopsis Seeds (Rosette Lawn h = 4.9 mm, Area: 15 cm²)'),
        Line2D([0], [0], color='#2980b9', lw=4.2, label='Port 3 Dipping Nozzle (Discharge at Y = 9.5 mm in Substrate 3)')
    ]
    fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.02),
               ncol=4, fontsize=9.0, frameon=True, fancybox=True, shadow=True)

    plt.savefig('figures/assembly_cross_sections_annotated.png', dpi=200, bbox_inches='tight')
    print("Saved figures/assembly_cross_sections_annotated.png")

    # =========================================================================
    # FIGURE 2: PURE CAD ENGINEERING CROSS SECTIONS WITH KEY DIMENSIONS (EXACT)
    # =========================================================================
    fig2, axes2 = plt.subplots(1, 3, figsize=(22, 7.5), facecolor='white')
    plt.subplots_adjust(left=0.05, right=0.97, bottom=0.12, top=0.88, wspace=0.22)
    
    # Slice 1: X = -5.0 mm
    ax_e1 = axes2[0]
    ax_e1.set_facecolor('#ffffff')
    plot_mesh_section(ax_e1, pot, [-0.005, 0, 0], [1, 0, 0], 2, 1, '#1b2631', lw=2.2, label='Pot Body')
    plot_mesh_section(ax_e1, lid_mated, [-0.005, 0, 0], [1, 0, 0], 2, 1, '#c0392b', lw=2.2, label='Lid Assembly')
    ax_e1.axhline(33.67, color='#7f8c8d', linestyle=':', lw=1.2)
    ax_e1.text(34.0, 33.67, 'Mating Plane (Y = 33.67 mm)', fontsize=7.5, color='#7f8c8d', va='center')
    ax_e1.axhline(7.62, color='#b9770e', linestyle='--', lw=1.0)
    ax_e1.text(34.0, 7.62, 'Substrate Base Ledge (Y = 7.62 mm)', fontsize=7.5, color='#b9770e', va='center')
    ax_e1.set_xlim([-37, 48])
    ax_e1.set_ylim([-2, 44])
    ax_e1.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax_e1.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax_e1.set_title('CAD Mated Geometry: Slice at X = -5.0 mm\n(Longitudinal Section through Compartments)', fontsize=12, fontweight='bold', pad=10)
    ax_e1.grid(True, linestyle='--', alpha=0.4)
    ax_e1.set_aspect('equal')
    ax_e1.legend(loc='lower left', fontsize=8.5)

    # Slice 2: X = -15.49 mm (Ports)
    ax_e2 = axes2[1]
    ax_e2.set_facecolor('#ffffff')
    plot_mesh_section(ax_e2, pot, [-0.015494, 0, 0], [1, 0, 0], 2, 1, '#1b2631', lw=2.2, label='Pot Body')
    plot_mesh_section(ax_e2, lid_mated, [-0.015494, 0, 0], [1, 0, 0], 2, 1, '#c0392b', lw=2.2, label='Lid Assembly')
    ax_e2.annotate('Port 1 (Ø 3.18)', xy=(-21.59, 38.1), xytext=(-21.59, 43.0), ha='center', va='bottom', fontsize=8, fontweight='bold', color='#c0392b', arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.4))
    ax_e2.annotate('Port 2 (Ø 3.18)', xy=(0.0, 38.1), xytext=(0.0, 43.0), ha='center', va='bottom', fontsize=8, fontweight='bold', color='#8e44ad', arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.4))
    ax_e2.annotate('Port 3 (Ø 3.18)', xy=(21.59, 38.1), xytext=(21.59, 43.0), ha='center', va='bottom', fontsize=8, fontweight='bold', color='#1f618d', arrowprops=dict(arrowstyle='->', color='#1f618d', lw=1.4))
    ax_e2.set_xlim([-37, 37])
    ax_e2.set_ylim([-2, 48])
    ax_e2.set_xlabel('Z Length (mm)', fontsize=11, fontweight='bold')
    ax_e2.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax_e2.set_title('CAD Mated Geometry: Slice at X = -15.49 mm\n(Passing Directly Through Fluid Ports)', fontsize=12, fontweight='bold', pad=10)
    ax_e2.grid(True, linestyle='--', alpha=0.4)
    ax_e2.set_aspect('equal')
    ax_e2.legend(loc='lower left', fontsize=8.5)

    # Slice 3: Z = 0.0 mm
    ax_e3 = axes2[2]
    ax_e3.set_facecolor('#ffffff')
    plot_mesh_section(ax_e3, pot, [0, 0, 0], [0, 0, 1], 0, 1, '#1b2631', lw=2.2, label='Pot Body')
    plot_mesh_section(ax_e3, lid_mated, [0, 0, 0], [0, 0, 1], 0, 1, '#c0392b', lw=2.2, label='Lid Assembly')
    ax_e3.annotate('Port 2 Centerline\n(X = -15.49 mm)', xy=(-15.49, 38.1), xytext=(-15.49, 44.0), ha='center', va='bottom', fontsize=8, fontweight='bold', color='#8e44ad', arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=1.4))
    ax_e3.set_xlim([-25, 25])
    ax_e3.set_ylim([-2, 50])
    ax_e3.set_xlabel('X Width (mm)', fontsize=11, fontweight='bold')
    ax_e3.set_ylabel('Y Height (mm)', fontsize=11, fontweight='bold')
    ax_e3.set_title('CAD Mated Geometry: Transverse at Z = 0.0 mm\n(Across Central Compartment)', fontsize=12, fontweight='bold', pad=10)
    ax_e3.grid(True, linestyle='--', alpha=0.4)
    ax_e3.set_aspect('equal')
    ax_e3.legend(loc='lower left', fontsize=8.5)

    plt.savefig('figures/assembly_cross_sections_exact.png', dpi=200, bbox_inches='tight')
    print("Saved figures/assembly_cross_sections_exact.png")

if __name__ == '__main__':
    generate_assembly_cross_sections()
