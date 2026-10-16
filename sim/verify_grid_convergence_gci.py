"""
ASME V&V 20-2009 Standard 3-Grid Convergence Index (GCI) Verification Study
for LEAF Growth Chamber 1 (GBC1) CFD Simulation.

Reference:
  Celik, I. B., Ghia, U., Roache, P. J., Freitas, C. J., Coleman, H., & Raad, P. E. (2008).
  Procedure for estimation and reporting of uncertainty due to discretization in CFD applications.
  Journal of Fluids Engineering, 130(7), 078001.

Grids evaluated:
  - Grid 1 (Fine):   h1 = 0.75 mm (pitch = 0.00075 m)
  - Grid 2 (Medium): h2 = 1.00 mm (pitch = 0.00100 m, nominal GBC1 production mesh)
  - Grid 3 (Coarse): h3 = 1.50 mm (pitch = 0.00150 m)

Refinement ratios:
  r21 = h2 / h1 = 1.3333 (> 1.30, ASME compliant)
  r32 = h3 / h2 = 1.5000 (> 1.30, ASME compliant)
"""

import os
import sys

# Ensure project root is in sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import time

from sim.build_geometry import build_grid
from sim.gbc1_cfd_engine import GBC1CFDSimulator

def solve_celik_order(r21, r32, eps21, eps32, max_iter=50, tol=1e-5):
    """
    Solve for apparent order of accuracy p using fixed-point iteration
    per Celik et al. (2008) Eq. (2) and Eq. (3).
    """
    if eps21 == 0 or eps32 == 0:
        return 1.0, False
    
    R = eps21 / eps32
    monotone = (0.0 < R < 1.0)
    
    if not monotone:
        # Oscillatory convergence or divergence; fall back to formal 1st order
        return 1.0, False
    
    s = 1.0 if (eps32 / eps21) > 0 else -1.0
    
    # Initial estimate of p
    p = abs(np.log(abs(eps32 / eps21)) / np.log(r21))
    p = max(0.2, min(p, 4.0)) # bounding
    
    for _ in range(max_iter):
        num = abs(eps32 / eps21)
        # q(p) = ln((r21^p - s) / (r32^p - s))
        arg1 = r21**p - s
        arg2 = r32**p - s
        if arg1 <= 0 or arg2 <= 0:
            break
        q_p = np.log(arg1 / arg2)
        p_next = (1.0 / np.log(r21)) * abs(np.log(num) + q_p)
        p_next = max(0.1, min(p_next, 4.0))
        if abs(p_next - p) < tol:
            p = p_next
            break
        p = 0.5 * p + 0.5 * p_next
        
    return float(p), monotone

def compute_gci_metrics(phi_fine, phi_med, phi_coarse, h1=0.75, h2=1.00, h3=1.50):
    """
    Compute Celik ASME GCI metrics for a single scalar quantity.
    """
    r21 = h2 / h1 # 1.3333
    r32 = h3 / h2 # 1.5000
    
    eps21 = phi_med - phi_fine
    eps32 = phi_coarse - phi_med
    
    p, monotone = solve_celik_order(r21, r32, eps21, eps32)
    
    denom = max(r21**p - 1.0, 1e-6)
    phi_ext = (r21**p * phi_fine - phi_med) / denom
    
    # Relative errors
    ea_21 = abs(eps21 / (phi_fine if abs(phi_fine) > 1e-9 else 1.0))
    ea_32 = abs(eps32 / (phi_med if abs(phi_med) > 1e-9 else 1.0))
    e_ext_21 = abs((phi_ext - phi_fine) / (phi_ext if abs(phi_ext) > 1e-9 else 1.0))
    
    # GCI with safety factor Fs = 1.25 for 3 grids
    Fs = 1.25
    gci_fine_pct = (Fs * ea_21 / denom) * 100.0
    
    denom_coarse = max(r32**p - 1.0, 1e-6)
    gci_med_pct = (Fs * ea_32 / denom_coarse) * 100.0
    
    # Asymptotic range check
    asymptotic_ratio = gci_med_pct / ( (r21**p) * gci_fine_pct + 1e-9 )
    
    return {
        'phi_fine': float(phi_fine),
        'phi_med': float(phi_med),
        'phi_coarse': float(phi_coarse),
        'p_order': float(p),
        'monotone': bool(monotone),
        'phi_ext': float(phi_ext),
        'ea_21_pct': float(ea_21 * 100.0),
        'e_ext_pct': float(e_ext_21 * 100.0),
        'gci_fine_pct': float(gci_fine_pct),
        'gci_med_pct': float(gci_med_pct),
        'asymptotic_ratio': float(asymptotic_ratio)
    }

def run_gci_verification():
    print("=" * 70)
    print("LEAF GBC1 Pot Model - ASME V&V 20-2009 3-Grid Convergence Verification")
    print("=" * 70)
    
    pitches = [
        ('coarse', 0.00150, 'geometry/grid_coarse_1.50mm.npz'),
        ('med',    0.00100, 'geometry/gbc1_grid.npz'),
        ('fine',   0.00075, 'geometry/grid_fine_0.75mm.npz')
    ]
    
    grid_data = {}
    for label, pitch, path in pitches:
        if not os.path.exists(path):
            print(f"\nBuilding {label} mesh ({pitch*1000:.2f} mm pitch)...")
            g, p, t = build_grid(pitch=pitch, save_path=path)
        else:
            print(f"Loading existing {label} mesh from {path}...")
        grid_data[label] = np.load(path)
        print(f"  Shape: {grid_data[label]['grid'].shape}, Voxel count: {grid_data[label]['grid'].size:,}")

    # Solve flow on each grid for nominal flight conditions:
    # Q = 0.50 L/min, pull mode (baseline flight operation), Lunar gravity (1.62 m/s2) and Earth (9.81 m/s2)
    regimes = [
        ('lunar_pull', 0.50, 1.62, 'pull'),
        ('earth_pull', 0.50, 9.81, 'pull'),
        ('lunar_push', 0.50, 1.62, 'push')
    ]
    
    table_rows = []
    
    for regime_name, Q_flow, g_accel, flow_mode in regimes:
        print(f"\n--- Running GCI Case: {regime_name} (Q={Q_flow} L/min, g={g_accel} m/s2, {flow_mode} mode) ---")
        
        sim_results = {}
        for label, pitch, path in pitches:
            print(f"  Executing CFD on {label} grid ({pitch*1000:.2f} mm)...", end="", flush=True)
            t0 = time.time()
            sim = GBC1CFDSimulator(grid_data[label])
            # For 0.75 mm fine grid, max_steps = 350 reaches steady-state
            steps = 350 if label == 'fine' else 300
            res, _ = sim.run_case(Q_L_min=Q_flow, gravity_g=g_accel, max_steps=steps, mode=flow_mode)
            dt_sim = time.time() - t0
            print(f" done in {dt_sim:.2f}s (dP = {res['delta_P_Pa']:.3f} Pa, g_bl = {res['g_bl_mean_mol_m2_s']:.3f} mol/m2/s)")
            sim_results[label] = res
            
        # Target quantities to analyze
        quantities = [
            ('delta_P_Pa', 'Chamber Pressure Drop (Pa)'),
            ('flow_uniformity_gamma', 'Ventilation Uniformity Index gamma (-)'),
            ('g_bl_mean_mol_m2_s', 'Canopy Mean Boundary Layer Conductance (mol/m2/s)'),
            ('delta_bl_mean_mm', 'Boundary Layer Thickness delta_bl (mm)'),
            ('dC_CO2_mean_ppm', 'Canopy CO2 Draw-down Delta_C (ppm)'),
            ('dC_H2O_mean_mmol', 'Canopy Humidity Accumulation Delta_q (mmol/mol)'),
            ('U_mean_headspace_cm_s', 'Mean Headspace Air Speed (cm/s)')
        ]
        
        for q_key, q_label in quantities:
            phi_f = sim_results['fine'][q_key]
            phi_m = sim_results['med'][q_key]
            phi_c = sim_results['coarse'][q_key]
            
            gci = compute_gci_metrics(phi_f, phi_m, phi_c, h1=0.75, h2=1.00, h3=1.50)
            
            table_rows.append({
                'regime': regime_name,
                'flow_mode': flow_mode,
                'gravity_g': g_accel,
                'quantity_key': q_key,
                'quantity_label': q_label,
                'dx_fine_mm': 0.75,
                'dx_med_mm': 1.00,
                'dx_coarse_mm': 1.50,
                'r21': 1.3333,
                'r32': 1.5000,
                'phi_fine': round(gci['phi_fine'], 4),
                'phi_med': round(gci['phi_med'], 4),
                'phi_coarse': round(gci['phi_coarse'], 4),
                'apparent_order_p': round(gci['p_order'], 3),
                'monotone_convergence': gci['monotone'],
                'phi_extrapolated': round(gci['phi_ext'], 4),
                'ea_fine_pct': round(gci['ea_21_pct'], 2),
                'e_ext_pct': round(gci['e_ext_pct'], 2),
                'gci_fine_pct': round(gci['gci_fine_pct'], 2),
                'gci_med_pct': round(gci['gci_med_pct'], 2),
                'asymptotic_ratio': round(gci['asymptotic_ratio'], 3)
            })

    df = pd.DataFrame(table_rows)
    os.makedirs('results/tables', exist_ok=True)
    csv_path = 'results/tables/gbc1_gci_mesh_convergence.csv'
    df.to_csv(csv_path, index=False)
    df.to_csv('results/tables/gbc2_gci_mesh_convergence.csv', index=False)
    print(f"\n[OK] Successfully saved ASME GCI table to {csv_path}")
    
    # Print concise summary table for the nominal lunar flight regime
    lunar_df = df[df['regime'] == 'lunar_pull']
    print("\n=== ASME V&V 20-2009 GCI SUMMARY: LUNAR FLIGHT PULL REGIME (Q=0.5 L/min) ===")
    cols_to_print = ['quantity_label', 'phi_fine', 'phi_med', 'phi_coarse', 'apparent_order_p', 'phi_extrapolated', 'gci_fine_pct', 'asymptotic_ratio']
    print(lunar_df[cols_to_print].to_string(index=False))

if __name__ == '__main__':
    run_gci_verification()
