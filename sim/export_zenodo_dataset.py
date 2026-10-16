"""
Zenodo Open-Science Dataset Packager for LEAF GBC1 'Model Lunar LEAF' 5-Day Spaceflight Experiment.
Complies with FAIR data principles and Frictionless Data standards.
"""

import os
import json
import hashlib
import shutil
import pandas as pd

def compute_sha256(filepath):
    hasher = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192):
            hasher.update(chunk)
    return hasher.hexdigest()

def export_zenodo_package():
    print("\n--- Exporting Zenodo Open-Science Data Package ---")
    zenodo_dir = 'zenodo'
    os.makedirs(zenodo_dir, exist_ok=True)
    os.makedirs(os.path.join(zenodo_dir, 'data'), exist_ok=True)
    os.makedirs(os.path.join(zenodo_dir, 'figures'), exist_ok=True)
    
    # 1. Copy CSV datasets
    csv_files = [
        'results/tables/zenodo_5day_hourly_telemetry.csv',
        'results/tables/zenodo_push_vs_pull_aerodynamics.csv',
        'results/tables/zenodo_gravity_temp_response_surface.csv',
        'results/tables/zenodo_water_activation_foam_nytex.csv',
        'results/tables/gbc1_gci_mesh_convergence.csv',
        'results/tables/gbc1_foliar_surface_distribution.csv',
        'results/tables/gbc1_164plants_foliar_telemetry.csv',
        'results/tables/gbc1_spaceflight_hardware_comprehensive_benchmark.csv'
    ]
    
    for f in csv_files:
        if os.path.exists(f):
            shutil.copy2(f, os.path.join(zenodo_dir, 'data', os.path.basename(f)))
            
    # 2. Copy Figures 12-18
    fig_files = [
        'figures/fig12_water_activation_foam_and_nytex.png',
        'figures/fig12_water_activation_foam_and_nytex.svg',
        'figures/fig13_aerodynamic_push_vs_pull_regimes.png',
        'figures/fig13_aerodynamic_push_vs_pull_regimes.svg',
        'figures/fig14_model_lunar_leaf_5day_ontogeny.png',
        'figures/fig14_model_lunar_leaf_5day_ontogeny.svg',
        'figures/fig15_gravity_temperature_response_surfaces.png',
        'figures/fig15_gravity_temperature_response_surfaces.svg',
        'figures/fig16_optical_lid_condensation_envelope.png',
        'figures/fig16_optical_lid_condensation_envelope.svg',
        'figures/fig17_foliar_surface_distribution_and_hotspots.png',
        'figures/fig17_foliar_surface_distribution_and_hotspots.svg',
        'figures/fig18_spaceflight_hardware_comprehensive_benchmark.png',
        'figures/fig18_spaceflight_hardware_comprehensive_benchmark.svg'
    ]
    
    for f in fig_files:
        if os.path.exists(f):
            shutil.copy2(f, os.path.join(zenodo_dir, 'figures', os.path.basename(f)))
            
    # 3. Compute Checksums
    checksums = {}
    for root, _, files in os.walk(zenodo_dir):
        for file in files:
            if file == 'checksums.txt' or file == 'datapackage.json':
                continue
            p = os.path.join(root, file)
            rel = os.path.relpath(p, zenodo_dir)
            checksums[rel] = compute_sha256(p)
            
    with open(os.path.join(zenodo_dir, 'checksums.txt'), 'w') as f:
        for path, digest in sorted(checksums.items()):
            f.write(f"{digest}  {path}\n")
            
    # 4. Generate Frictionless Data Package Descriptor (datapackage.json)
    datapackage = {
        "name": "leaf-gbc1-model-lunar-leaf-5day-simulation",
        "title": "LEAF GBC1 'Model Lunar LEAF' 5-Day Spaceflight Simulation Dataset",
        "description": "Comprehensive 5-day (120-hour) computational multi-physics dataset for the LEAF GBC1 spaceflight plant growth chamber. Encompasses subsurface water activation in porous foam, Nytex mesh root exclusion, push vs pull aerodynamic CFD, developmental ontogeny, mixed-convection gravity (0g-1g), temperature sweeps (18-30 C), and optical lid condensation defense.",
        "version": "1.0.0",
        "homepage": "https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/",
        "license": "CC-BY-4.0",
        "keywords": [
            "space-biology",
            "computational-fluid-dynamics",
            "microgravity",
            "lunar-gravity",
            "plant-growth-chamber",
            "porous-media",
            "ontogeny",
            "leaf-gbc1",
            "artemis"
        ],
        "resources": [
            {
                "name": "5day_hourly_telemetry",
                "path": "data/zenodo_5day_hourly_telemetry.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Hourly 120-hour multi-species state matrix covering canopy height, leaf area, drag factor, CO2 exchange, transpiration, and condensation margin."
            },
            {
                "name": "push_vs_pull_aerodynamics",
                "path": "data/zenodo_push_vs_pull_aerodynamics.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Comparative aerodynamic metrics across 10 pump speeds for Push (positive gauge) vs Pull (suction negative gauge) flow regimes."
            },
            {
                "name": "gravity_temp_response_surface",
                "path": "data/zenodo_gravity_temp_response_surface.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Full factorial response surface spanning 4 gravity levels (0g, 0.166g, 0.38g, 1.0g), 5 air temperatures (18-30 C), and 5 ventilation rates."
            },
            {
                "name": "water_activation_foam_nytex",
                "path": "data/zenodo_water_activation_foam_nytex.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Vertical moisture saturation profiles in the 7.0 mm porous foam cassette and Nytex mesh Young-Laplace bubble point curve."
            },
            {
                "name": "gbc1_gci_mesh_convergence",
                "path": "data/gbc1_gci_mesh_convergence.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "ASME V&V 20-2009 standard 3-grid convergence index verification study across 1.50, 1.00, and 0.75 mm pitches."
            },
            {
                "name": "gbc1_foliar_surface_distribution",
                "path": "data/gbc1_foliar_surface_distribution.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Canopy microclimate percentiles (p10, p50, p90), spatial CoV, and dynamic range across flight regimes."
            },
            {
                "name": "gbc1_164plants_foliar_telemetry",
                "path": "data/gbc1_164plants_foliar_telemetry.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Individual scalar microenvironment telemetry for all 164 plants (32 B. rapa + 100 Arabidopsis + 32 B. rapa)."
            },
            {
                "name": "gbc1_spaceflight_hardware_comprehensive_benchmark",
                "path": "data/gbc1_spaceflight_hardware_comprehensive_benchmark.csv",
                "format": "csv",
                "mediatype": "text/csv",
                "encoding": "utf-8",
                "description": "Multi-system spaceflight hardware comparative benchmark across NASA BRIC, CARA, VEGGIE, and LEAF GBC1."
            }
        ]
    }
    
    with open(os.path.join(zenodo_dir, 'datapackage.json'), 'w') as f:
        json.dump(datapackage, f, indent=2)
        
    # 5. Generate Zenodo README.md Data Dictionary
    readme_content = """# LEAF GBC1 "Model Lunar LEAF" 5-Day Spaceflight Simulation Dataset

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.leaf--gbc1--lunar-blue)](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This archive contains the open-science research data, CFD metrics, and developmental biophysical time series for the **Model Lunar LEAF** 5-day (120-hour) spaceflight simulation conducted in the **LEAF GBC1** plant growth hardware (`M09024F-MDL-5` pot and `M09025F-MDL-2` lid assembly).

## Repository & Interactive Portal
- **Interactive Research Summary Portal**: https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/
- **Standalone 3D WebGL Digital Twin**: https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/viewer_3d.html
- **Full Manuscript & Equations**: `MANUSCRIPT_LEAF_GBC1_CFD.md`

---

## Dataset Inventory (`data/`)

### 1. `zenodo_5day_hourly_telemetry.csv`
Contains 121 hourly records (Hours 0 to 120) tracking multi-species developmental ontogeny under 16h photoperiod / 8h skotoperiod:
- `mission_hour`: Elapsed hours from water activation (0 to 120 h).
- `mission_day`: Mission day fraction (0.0 to 5.0).
- `is_daylight`: Boolean indicator for photoperiod LED state (True = Day, False = Night).
- `h_brapa_mm`: *Brassica rapa* canopy height above substrate (mm).
- `h_arab_mm`: *Arabidopsis thaliana* rosette height above substrate (mm).
- `A_brapa_cm2`: Single-chamber *B. rapa* foliar surface area (cm² across 32 plants).
- `A_arab_cm2`: Central-chamber *Arabidopsis* foliar surface area (cm² across 100 plants).
- `A_total_foliar_cm2`: Total whole-chamber foliar area (2x B. rapa + 1x Arabidopsis).
- `drag_factor_brapa_s_1`: Canopy aerodynamic Brinkman drag coefficient for *B. rapa* (s⁻¹).
- `drag_factor_arab_s_1`: Canopy aerodynamic Brinkman drag coefficient for *Arabidopsis* (s⁻¹).
- `A_net_chamber_umol_s`: Whole-chamber net CO₂ assimilation rate (µmol CO₂ / s).
- `delta_CO2_ppm`: Inter-chamber CO₂ concentration drawdown / elevation relative to inlet (ppm).
- `transpiration_rate_mL_h`: Hourly whole-canopy moisture transpiration rate (mL / h).
- `cum_water_transpired_mL`: Cumulative water transpired from water activation (mL).
- `RH_chamber_pct`: Steady-state chamber relative humidity (%).
- `T_dew_C`: Chamber dew point temperature (°C).
- `T_lid_C`: Acrylic optical lid inner surface temperature (°C).
- `delta_T_margin_C`: Optical lid anti-fogging safety margin (°C, must be > +0.50 °C).

### 2. `zenodo_push_vs_pull_aerodynamics.csv`
Comparative aerodynamic data evaluating Push Mode (forced inlet at Port 1) vs. Pull Mode (suction at Port 2) across 10 flow rates:
- `Q_L_min`: Volumetric ventilation rate (0.05 to 1.50 L/min).
- `U_port_m_s`: Linear gas velocity through 1/8" ID fluid port (m/s).
- `tau_residence_s`: Mean headspace flush residence time (s).
- `push_delta_P_Pa` / `pull_delta_P_Pa`: Total aerodynamic pressure drop across loop (Pa).
- `push_p_gauge_Pa` / `pull_p_gauge_Pa`: Internal static gauge pressure (Pa).
- `push_g_bl_mean_mol_m2_s` / `pull_g_bl_mean_mol_m2_s`: Mean boundary layer conductance (mol/m² s).
- `push_uniformity_gamma` / `pull_uniformity_gamma`: Inter-chamber flow uniformity coefficient γ.

### 3. `zenodo_gravity_temp_response_surface.csv`
100-case full factorial matrix crossing 4 gravity regimes (0g, 0.166g, 0.38g, 1.0g), 5 temperatures (18, 21, 24, 27, 30 °C), and 5 ventilation rates:
- `gravity_g_rel`: Relative gravity (0.0, 0.166, 0.38, 1.0).
- `T_air_C`: Chamber bulk air temperature (°C).
- `Q_L_min`: Ventilation rate (L/min).
- `Re`: Canopy Reynolds number.
- `Gr`: Thermal Grashof number.
- `Ri`: Richardson number (Ri = Gr / Re²).
- `convection_regime`: Classification (Forced, Mixed, or Natural Convection).
- `g_bl_mol_m2_s`: Effective boundary layer conductance (mol/m² s).
- `r_b_s_m`: Boundary layer resistance (s/m).
- `A_net_brapa_umol_m2_s` / `A_net_arab_umol_m2_s`: Net canopy assimilation rate (µmol/m² s).
- `delta_T_margin_C`: Optical lid anti-fogging margin (°C).
- `condensation_risk`: Boolean flag for window fogging hazard.

### 4. `zenodo_water_activation_foam_nytex.csv`
Vertical saturation and aeration profile across the 7.0 mm porous foam cassette (Y = 7.62 to 14.62 mm):
- `elevation_in_foam_mm`: Local elevation from foam base (0.0 to 7.0 mm).
- `elevation_CAD_Y_mm`: Absolute CAD coordinate Y (7.62 to 14.62 mm).
- `theta_w_Earth_1g` / `theta_w_Lunar_0.166g`: Volumetric water content (m³/m³).
- `eps_air_Earth_1g` / `eps_air_Lunar_0.166g`: Air-filled porosity (m³/m³).
- `D_eff_O2_Earth_m2_s` / `D_eff_O2_Lunar_m2_s`: Effective O₂ diffusivity in substrate (m²/s).

---

## Publication Figures (`figures/`)
- `fig12_water_activation_foam_and_nytex.png / .svg`: Water injection, Richards saturation curve, Nytex root barrier.
- `fig13_aerodynamic_push_vs_pull_regimes.png / .svg`: Push vs Pull CFD velocity fields, pressure drop, containment.
- `fig14_model_lunar_leaf_5day_ontogeny.png / .svg`: 0-120h developmental ontogeny, height, LAI, CO₂ balance.
- `fig15_gravity_temperature_response_surfaces.png / .svg`: Gravity x Temperature response surfaces, Richardson number.
- `fig16_optical_lid_condensation_envelope.png / .svg`: 5-day lid temperature, dew point, anti-fogging envelope.

---

## Citation & Reproducibility
Data generated with the LEAF GBC1 computational suite (`sim/simulate_model_lunar_leaf_5day.py`).
License: Creative Commons Attribution 4.0 International (CC-BY 4.0).
"""
    with open(os.path.join(zenodo_dir, 'README.md'), 'w') as f:
        f.write(readme_content)
        
    print(f"  -> Successfully compiled Zenodo open-science archive in '{zenodo_dir}/'.")
    print(f"  -> Generated checksums for {len(checksums)} files.")
    print(f"  -> Created Frictionless 'datapackage.json' and data dictionary 'README.md'.")

if __name__ == '__main__':
    export_zenodo_package()
