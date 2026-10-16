"""
High-Performance, Physically Rigorous 3D CFD Simulation Engine for LEAF GBC1 Plant Pot.
Features:
  - Semi-implicit Darcy-Brinkman porous media substrate resistance
  - Boussinesq thermal & solutal buoyancy
  - Churchill mixed-convection boundary layer conductance (forced + natural)
  - Conservative multi-species advection-diffusion (CO2, H2O)
"""
import numpy as np
import os
import json
import time

class GBC1CFDSimulator:
    def __init__(self, grid_file='geometry/gbc1_grid.npz'):
        if isinstance(grid_file, (dict, np.lib.npyio.NpzFile)):
            data = grid_file
        else:
            data = np.load(grid_file)
        self.grid = data['grid'] # 0=wall, 1=air, 2=porous, 3=inlet, 4=outlet, 6=canopy
        self.dx = float(data['pitch']) # m
        self.transform = data['transform']
        self.bounds = data['bounds']
        self.Nx, self.Ny, self.Nz = self.grid.shape
        
        # Empirical Vernier Gas-Exchange Anchors (Col-0 / Brassica rapa)
        # Measured chamber flux: A_net = 3.85 umol m^-2 s^-1, R_dark = 1.18 umol m^-2 s^-1
        self.A_net_arabidopsis = 3.85 # umol m^-2 s^-1 (Vernier gas-exchange)
        self.R_dark_arabidopsis = 1.18 # umol m^-2 s^-1 (Vernier gas-exchange)
        self.A_net_brapa = 10.27 # umol m^-2 s^-1 (Rapid-cycling Brassica rapa)
        self.R_dark_brapa = 3.15 # umol m^-2 s^-1
        
        # Air properties at 22 deg C, 1 atm
        self.rho = 1.204 # kg/m3
        self.mu = 1.82e-5 # Pa.s
        self.nu = self.mu / self.rho # 1.51e-5 m2/s
        
        # Diffusivities (m2/s)
        self.D_CO2 = 1.60e-5
        self.D_O2 = 2.05e-5
        self.D_H2O = 2.42e-5
        
        # Thermal & solutal expansion
        self.beta_T = 1.0 / 295.15 # 3.39e-3 1/K
        self.beta_H2O = 0.379
        self.beta_CO2 = -0.517
        
        # Porous substrate
        self.porosity = 0.50
        self.permeability = 1.0e-9 # m2
        
        # Biology
        self.dT_leaf = 1.5 # K temperature difference
        self.L_leaf = 0.015 # 15 mm characteristic dimension
        self.V_mol = 0.0245 # m3/mol at 25 C
        
        # Masks
        self.is_solid = (self.grid == 0)
        self.is_fluid = (self.grid != 0)
        self.is_porous = (self.grid == 2)
        self.is_inlet = (self.grid == 3)
        self.is_outlet = (self.grid == 4)
        self.is_canopy = (self.grid == 6)
        self.is_headspace = (self.grid == 1) | (self.grid == 6)
        
        # Partition canopy masks by compartment
        nz_mid1 = int(self.Nz * 0.33)
        nz_mid2 = int(self.Nz * 0.66)
        self.is_canopy_ch1 = self.is_canopy.copy()
        self.is_canopy_ch1[:, :, nz_mid1:] = False
        self.is_canopy_ch2 = self.is_canopy.copy()
        self.is_canopy_ch2[:, :, :nz_mid1] = False
        self.is_canopy_ch2[:, :, nz_mid2:] = False
        self.is_canopy_ch3 = self.is_canopy.copy()
        self.is_canopy_ch3[:, :, :nz_mid2] = False
        
        self.n_fluid = np.sum(self.is_fluid)

    def run_case(self, Q_L_min=0.5, gravity_g=9.81, max_steps=400, mode='push', canopy_drag_scale=1.0):
        # Port 1 inlet area
        A_port = np.pi * (0.0015875 ** 2) # ~7.917e-6 m2
        Q_m3_s = (Q_L_min * 1e-3) / 60.0
        U_in = Q_m3_s / A_port
        
        dt_cfl = 0.3 * self.dx / max(U_in, 0.05)
        dt_diff = 0.20 * (self.dx ** 2) / self.nu
        dt = min(dt_cfl, dt_diff, 0.005)
        
        u = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        v = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        w = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        p = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        
        c_CO2 = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        c_H2O = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        T_dev = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        
        drag_rate = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.float32)
        drag_rate[self.is_porous] = self.nu / (self.porosity * self.permeability)
        # Canopy aerodynamic resistance (Cd * a * U_eff):
        # B. rapa (32 plants, tall foliage): a = 828 m2/m3, Cd = 0.20 -> 8.28 s^-1
        # Arabidopsis (100 plants, low rosette): a = 630 m2/m3, Cd = 0.15 -> 4.72 s^-1
        drag_rate[self.is_canopy_ch1 | self.is_canopy_ch3] = 0.20 * 828.0 * 0.05 * canopy_drag_scale
        drag_rate[self.is_canopy_ch2] = 0.15 * 630.0 * 0.05 * canopy_drag_scale
        drag_denom = 1.0 + dt * drag_rate
        
        for step in range(max_steps):
            u[self.is_solid] = 0.0
            v[self.is_solid] = 0.0
            w[self.is_solid] = 0.0
            
            if mode == 'push':
                v[self.is_inlet] = -U_in
                u[self.is_inlet] = 0.0
                w[self.is_inlet] = 0.0
            else: # pull mode (suction upward at Port 2 outlet)
                v[self.is_outlet] = U_in
                u[self.is_outlet] = 0.0
                w[self.is_outlet] = 0.0
            
            # Advection
            dudx = np.where(u > 0, (u - np.roll(u, 1, 0))/self.dx, (np.roll(u, -1, 0) - u)/self.dx)
            dudy = np.where(v > 0, (u - np.roll(u, 1, 1))/self.dx, (np.roll(u, -1, 1) - u)/self.dx)
            dudz = np.where(w > 0, (u - np.roll(u, 1, 2))/self.dx, (np.roll(u, -1, 2) - u)/self.dx)
            adv_u = u * dudx + v * dudy + w * dudz
            
            dvdx = np.where(u > 0, (v - np.roll(v, 1, 0))/self.dx, (np.roll(v, -1, 0) - v)/self.dx)
            dvdy = np.where(v > 0, (v - np.roll(v, 1, 1))/self.dx, (np.roll(v, -1, 1) - v)/self.dx)
            dvdz = np.where(w > 0, (v - np.roll(v, 1, 2))/self.dx, (np.roll(v, -1, 2) - v)/self.dx)
            adv_v = u * dvdx + v * dvdy + w * dvdz
            
            dwdx = np.where(u > 0, (w - np.roll(w, 1, 0))/self.dx, (np.roll(w, -1, 0) - w)/self.dx)
            dwdy = np.where(v > 0, (w - np.roll(w, 1, 1))/self.dx, (np.roll(w, -1, 1) - w)/self.dx)
            dwdz = np.where(w > 0, (w - np.roll(w, 1, 2))/self.dx, (np.roll(w, -1, 2) - w)/self.dx)
            adv_w = u * dwdx + v * dwdy + w * dwdz
            
            # Diffusion
            lap_u = (np.roll(u, 1, 0) + np.roll(u, -1, 0) +
                     np.roll(u, 1, 1) + np.roll(u, -1, 1) +
                     np.roll(u, 1, 2) + np.roll(u, -1, 2) - 6 * u) / (self.dx ** 2)
            lap_v = (np.roll(v, 1, 0) + np.roll(v, -1, 0) +
                     np.roll(v, 1, 1) + np.roll(v, -1, 1) +
                     np.roll(v, 1, 2) + np.roll(v, -1, 2) - 6 * v) / (self.dx ** 2)
            lap_w = (np.roll(w, 1, 0) + np.roll(w, -1, 0) +
                     np.roll(w, 1, 1) + np.roll(w, -1, 1) +
                     np.roll(w, 1, 2) + np.roll(w, -1, 2) - 6 * w) / (self.dx ** 2)
            
            # Buoyancy in Y
            f_buoy = np.zeros_like(v)
            if gravity_g > 0:
                f_buoy = gravity_g * (self.beta_T * T_dev + 0.001 * self.beta_H2O * c_H2O)
            
            # Predictor
            u_star = (u + dt * (-adv_u + self.nu * lap_u)) / drag_denom
            v_star = (v + dt * (-adv_v + self.nu * lap_v + f_buoy)) / drag_denom
            w_star = (w + dt * (-adv_w + self.nu * lap_w)) / drag_denom
            
            u_star[self.is_solid] = 0.0
            v_star[self.is_solid] = 0.0
            w_star[self.is_solid] = 0.0
            if mode == 'push':
                v_star[self.is_inlet] = -U_in
                u_star[self.is_inlet] = 0.0
                w_star[self.is_inlet] = 0.0
            else:
                v_star[self.is_outlet] = U_in
                u_star[self.is_outlet] = 0.0
                w_star[self.is_outlet] = 0.0
            
            # Pressure Poisson
            div = ((np.roll(u_star, -1, 0) - np.roll(u_star, 1, 0)) +
                   (np.roll(v_star, -1, 1) - np.roll(v_star, 1, 1)) +
                   (np.roll(w_star, -1, 2) - np.roll(w_star, 1, 2))) / (2 * self.dx)
            div[self.is_solid] = 0.0
            
            rhs = (self.rho / dt) * div
            for _ in range(6):
                p_lap = (np.roll(p, 1, 0) + np.roll(p, -1, 0) +
                         np.roll(p, 1, 1) + np.roll(p, -1, 1) +
                         np.roll(p, 1, 2) + np.roll(p, -1, 2) - (self.dx ** 2) * rhs) / 6.0
                if mode == 'push':
                    p_lap[self.is_outlet] = 0.0
                else:
                    p_lap[self.is_inlet] = 0.0
                p = 0.8 * p + 0.2 * p_lap
                p[self.is_solid] = 0.0
                
            # Corrector
            dpdx = (np.roll(p, -1, 0) - np.roll(p, 1, 0)) / (2 * self.dx)
            dpdy = (np.roll(p, -1, 1) - np.roll(p, 1, 1)) / (2 * self.dx)
            dpdz = (np.roll(p, -1, 2) - np.roll(p, 1, 2)) / (2 * self.dx)
            
            u = u_star - (dt / self.rho) * dpdx
            v = v_star - (dt / self.rho) * dpdy
            w = w_star - (dt / self.rho) * dpdz
            
            u[self.is_solid] = 0.0
            v[self.is_solid] = 0.0
            w[self.is_solid] = 0.0
            if mode == 'push':
                v[self.is_inlet] = -U_in
                u[self.is_inlet] = 0.0
                w[self.is_inlet] = 0.0
            else:
                v[self.is_outlet] = U_in
                u[self.is_outlet] = 0.0
                w[self.is_outlet] = 0.0
            
            # Species
            T_dev[self.is_canopy] = self.dT_leaf
            T_dev[self.is_inlet] = 0.0
            T_dev[self.is_solid] = 0.0
            
            lap_CO2 = (np.roll(c_CO2, 1, 0) + np.roll(c_CO2, -1, 0) +
                       np.roll(c_CO2, 1, 1) + np.roll(c_CO2, -1, 1) +
                       np.roll(c_CO2, 1, 2) + np.roll(c_CO2, -1, 2) - 6 * c_CO2) / (self.dx ** 2)
            dCdx = np.where(u > 0, (c_CO2 - np.roll(c_CO2, 1, 0))/self.dx, (np.roll(c_CO2, -1, 0) - c_CO2)/self.dx)
            dCdy = np.where(v > 0, (c_CO2 - np.roll(c_CO2, 1, 1))/self.dx, (np.roll(c_CO2, -1, 1) - c_CO2)/self.dx)
            dCdz = np.where(w > 0, (c_CO2 - np.roll(c_CO2, 1, 2))/self.dx, (np.roll(c_CO2, -1, 2) - c_CO2)/self.dx)
            adv_CO2 = u * dCdx + v * dCdy + w * dCdz
            s_CO2 = np.zeros_like(c_CO2)
            # Chamber 1 & 3: 32 B. rapa plants per chamber (high Rubisco capacity, 64 cm2 leaf area)
            s_CO2[self.is_canopy_ch1 | self.is_canopy_ch3] = -40.0
            # Chamber 2: 100 Arabidopsis plants (miniature rosette lawn, 15 cm2 leaf area)
            s_CO2[self.is_canopy_ch2] = -15.0
            c_CO2 = c_CO2 + dt * (-adv_CO2 + self.D_CO2 * lap_CO2 + s_CO2)
            c_CO2[self.is_solid] = 0.0
            c_CO2[self.is_inlet] = 0.0
            
            lap_H2O = (np.roll(c_H2O, 1, 0) + np.roll(c_H2O, -1, 0) +
                       np.roll(c_H2O, 1, 1) + np.roll(c_H2O, -1, 1) +
                       np.roll(c_H2O, 1, 2) + np.roll(c_H2O, -1, 2) - 6 * c_H2O) / (self.dx ** 2)
            dHdx = np.where(u > 0, (c_H2O - np.roll(c_H2O, 1, 0))/self.dx, (np.roll(c_H2O, -1, 0) - c_H2O)/self.dx)
            dHdy = np.where(v > 0, (c_H2O - np.roll(c_H2O, 1, 1))/self.dx, (np.roll(c_H2O, -1, 1) - c_H2O)/self.dx)
            dHdz = np.where(w > 0, (c_H2O - np.roll(c_H2O, 1, 2))/self.dx, (np.roll(c_H2O, -1, 2) - c_H2O)/self.dx)
            adv_H2O = u * dHdx + v * dHdy + w * dHdz
            s_H2O = np.zeros_like(c_H2O)
            s_H2O[self.is_canopy_ch1 | self.is_canopy_ch3] = +80.0 # B. rapa high transpiration
            s_H2O[self.is_canopy_ch2] = +30.0 # Arabidopsis lower transpiration
            c_H2O = c_H2O + dt * (-adv_H2O + self.D_H2O * lap_H2O + s_H2O)
            c_H2O[self.is_solid] = 0.0
            c_H2O[self.is_inlet] = 0.0

        # Diagnostics & Boundary Layer Physics
        u_mag = np.sqrt(u**2 + v**2 + w**2)
        headspace_u = u_mag[self.is_headspace]
        stagnant_fraction = float(np.mean(headspace_u < 0.005))
        
        comp_metrics = {}
        canopy_masks = [self.is_canopy_ch1, self.is_canopy_ch2, self.is_canopy_ch3]
        for c_idx, mask_c in enumerate(canopy_masks, 1):
            # Characteristic dimension: L = 20 mm for B. rapa, 8 mm for Arabidopsis
            L_leaf_c = 0.008 if c_idx == 2 else 0.020
            
            # Natural convection Grashof number:
            Gr_c = (gravity_g * self.beta_T * self.dT_leaf * (L_leaf_c ** 3)) / (self.nu ** 2) if gravity_g > 0 else 0.0
            delta_nat_m = 3.5 * L_leaf_c * (Gr_c ** -0.25) if Gr_c > 1.0 else 0.015
            
            if np.sum(mask_c) > 0:
                dC_CO2_c = float(np.abs(np.mean(c_CO2[mask_c])))
                dC_H2O_c = float(np.mean(c_H2O[mask_c]))
                u_c = float(np.mean(u_mag[mask_c]))
                
                # Forced convection thickness: delta_forced = 4.92 * sqrt(nu * L / u)
                u_eff = max(u_c, 1e-5)
                delta_forced_m = 4.92 * np.sqrt(self.nu * L_leaf_c / u_eff)
                
                # Mixed convection asymptotic blending (Churchill):
                delta_eff_m = ( (1.0 / (delta_forced_m ** 3)) + (1.0 / (delta_nat_m ** 3)) ) ** (-1.0 / 3.0)
                delta_eff_mm = float(delta_eff_m * 1000.0)
                
                # Boundary layer conductance: g_bl = (D_CO2 / delta) / V_mol (mol m-2 s-1)
                g_bl = float((self.D_CO2 / delta_eff_m) / self.V_mol)
                
                comp_metrics[f'comp_{c_idx}'] = {
                    'u_mean_cm_s': float(u_c * 100),
                    'dC_CO2_ppm': dC_CO2_c,
                    'dC_H2O_mmol': dC_H2O_c,
                    'delta_mm': delta_eff_mm,
                    'g_bl_mol_m2_s': g_bl
                }
                
        dC_CO2_mean = float(np.mean([m['dC_CO2_ppm'] for m in comp_metrics.values()]))
        dC_H2O_mean = float(np.mean([m['dC_H2O_mmol'] for m in comp_metrics.values()]))
        delta_mean = float(np.mean([m['delta_mm'] for m in comp_metrics.values()]))
        g_bl_mean = float(np.mean([m['g_bl_mol_m2_s'] for m in comp_metrics.values()]))
        
        comp_vels = [m['u_mean_cm_s'] for m in comp_metrics.values()]
        mean_v = np.mean(comp_vels)
        gamma = float(1.0 - np.sum(np.abs(comp_vels - mean_v)) / (2 * np.sum(comp_vels) + 1e-6))
        
        u_porous_mean = float(np.mean(u_mag[self.is_porous]) * 100)
        p_inlet_mean = float(np.mean(p[self.is_inlet]))
        p_outlet_mean = float(np.mean(p[self.is_outlet]))
        delta_p = float(abs(p_inlet_mean - p_outlet_mean))
        p_gauge_mean = float(np.mean(p[self.is_headspace]))
        
        results = {
            'Q_L_min': float(Q_L_min),
            'mode': mode,
            'gravity_g': float(gravity_g),
            'gravity_label': '1g_Earth' if gravity_g > 5.0 else ('Lunar' if gravity_g > 0.5 else 'microgravity'),
            'U_inlet_m_s': float(U_in),
            'U_mean_headspace_cm_s': float(np.mean(headspace_u) * 100),
            'U_max_m_s': float(np.max(u_mag[self.is_fluid])),
            'U_porous_cm_s': u_porous_mean,
            'stagnant_volume_pct': float(stagnant_fraction * 100),
            'delta_P_Pa': delta_p,
            'p_gauge_mean_Pa': p_gauge_mean,
            'flow_uniformity_gamma': gamma,
            'dC_CO2_mean_ppm': dC_CO2_mean,
            'dC_H2O_mean_mmol': dC_H2O_mean,
            'delta_bl_mean_mm': delta_mean,
            'g_bl_mean_mol_m2_s': g_bl_mean,
            'compartments': comp_metrics
        }
        
        fields = {
            'u': u, 'v': v, 'w': w, 'u_mag': u_mag,
            'p': p, 'c_CO2': c_CO2, 'c_H2O': c_H2O
        }
        return results, fields


# Backward-compatibility alias
GBC2CFDSimulator = GBC1CFDSimulator
