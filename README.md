# LEAF GBC2 Lunar Plant Growth Chamber — 3D CFD & Multi-Species Biophysics Model

[![3D Digital Twin](https://img.shields.io/badge/Interactive-3D_Digital_Twin-27ae60?style=for-the-badge)](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/viewer_3d.html)
[![GitHub Pages Portal](https://img.shields.io/badge/Live_Portal-GitHub_Pages-2980b9?style=for-the-badge)](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/)
[![Technical Manuscript](https://img.shields.io/badge/Manuscript-Full_Text_Report-8e44ad?style=for-the-badge)](./MANUSCRIPT_LEAF_GBC2_CFD.md)

Computational Fluid Dynamics (CFD), multi-species gas exchange, 14-day developmental canopy ontogeny, diurnal photoperiod inversion, off-nominal spacecraft fault-tolerance, and harvest-phase RNAlater fixative perfusion modeling for the **LEAF GBC2 spaceflight plant growth pot (`M09024F-MDL-5.STEP`) and lid (`M09025F-MDL-2.STEP`) assembly** under **Earth (1.0g, 9.81 m/s²)** and **Lunar (0.166g, 1.62 m/s²)** gravity regimes.

---

## Quick Links

- **Live GitHub Pages Research Portal**: [**`https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/`**](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/) ([`index.html`](./index.html))
- **Live Interactive 3D Chamber Digital Twin & CFD Simulator**: [**`https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/viewer_3d.html`**](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/viewer_3d.html) ([`viewer_3d.html`](./viewer_3d.html))
- **Full Technical Manuscript (Methods, Results & Figure Legends)**: [`MANUSCRIPT_LEAF_GBC2_CFD.md`](./MANUSCRIPT_LEAF_GBC2_CFD.md)
- **Orthographic Blueprint (ABAI QC Annotated)**: [`figures/fig_orthographic_projection_annotated.png`](./figures/fig_orthographic_projection_annotated.png) | [`.svg`](./figures/fig_orthographic_projection_annotated.svg)
- **Mated CAD Assembly Cross-Sections**: [`figures/assembly_cross_sections_annotated.png`](./figures/assembly_cross_sections_annotated.png) | [`exact`](./figures/assembly_cross_sections_exact.png)
- **Quantitative Telemetry Tables**: [`results/tables/`](./results/tables/)

---

## 1. Hardware Architecture & CAD Geometry

![Orthographic Projections Annotated](figures/fig_orthographic_projection_annotated.png)

### 1.1 Mated STEP CAD Envelope
- **Pot Body (`cad/pot.step` / `M09024F-MDL-5`)**: 47.63 mm (X width) × 33.67 mm (Y height) × 67.95 mm (Z length).
- **Lid Assembly (`cad/lid.step` / `M09025F-MDL-2`)**: 45.09 mm (X width) × 5.71 mm (Y height) × 65.41 mm (Z length), seated on the pot rim at _Y_ = 33.67 mm (_Y_ = 32.40–38.11 mm).
- **Total Internal Fluid Volume**: 64.60 cm³ (headspace air volume = 54.05 cm³; porous media bulk volume = 10.55 cm³, void volume = 5.28 mL at porosity ε = 0.50).
- **Lower Aeration & Drainage Wells (_Y_ = 0.00–7.62 mm)**: Molded basal channels (_V_drain = 8.74 cm³) beneath the chamber retainer ledges providing root aeration and capillary buffering.
- **3-Part Partitioned Porous Substrate (_Y_ = 7.62–14.62 mm)**: Three discrete porous media blocks (Δ_Y_ = 7.00 mm, permeability _K_ = 1.0 × 10⁻¹⁰ m², porosity ε = 0.50) seated inside each growth compartment and separated by solid vertical baffles (_Z_ = -9.5 to -7.5 mm and _Z_ = +7.5 to +9.5 mm).

### 1.2 Fluid Ports (_X_ = -15.49 mm, ⌀ 3.18 mm / 1/8″ ID)
1. **Port 1 (_Z_ = -21.59 mm) — Forced Air Inlet**: Injects conditioned air from the pump manifold into the upper perimeter distribution channel (_Y_ = 28–32 mm).
2. **Port 2 (_Z_ = 0.00 mm) — Central Air Exhaust**: Positioned directly above central Chamber 2 to collect converging airflow from both lateral chambers and return it to the cabin loop (also serves as the displaced gas vent during harvest fixation).
3. **Port 3 (_Z_ = +21.59 mm) — Water & RNAlater Delivery Port**: Equipped with an internal vertical dipping nozzle (_L_ = 24.2 mm) extending down to **_Y_ = 9.5 mm inside Substrate Block 3** to ensure bottom-up capillary wetting and eliminate foliar splash damage.

![Mated CAD Cross Sections](figures/assembly_cross_sections_annotated.png)

---

## 2. Multi-Species Planting Densities & Canopy Morphometrics

| Parameter | Chamber 1 (*Brassica rapa*) | Chamber 2 (*Arabidopsis thaliana*) | Chamber 3 (*Brassica rapa*) |
| :--- | :--- | :--- | :--- |
| **Axial Span (_Z_)** | -30.5 to -9.5 mm | -7.5 to +7.5 mm (Central) | +9.5 to +30.5 mm |
| **Seed Count & Density** | **32 seeds** (47.2 seeds/dm²) | **100 seeds** (206.2 seeds/dm²) | **32 seeds** (47.2 seeds/dm²) |
| **Single-Plant Leaf Area** | 2.00 cm² | 0.15 cm² | 2.00 cm² |
| **Total Chamber Foliar Area** | **64.0 cm²** (0.0064 m²) | **15.0 cm²** (0.0015 m²) | **64.0 cm²** (0.0064 m²) |
| **Canopy Vertical Span (_Y_)** | 15.0–26.0 mm (_h_ = 11.4 mm) | 15.0–19.5 mm (_h_ = 4.9 mm) | 15.0–26.0 mm (_h_ = 11.4 mm) |
| **Leaf Area Index (LAI)** | 9.44 | 3.09 | 9.44 |
| **Leaf Area Density (_a_)** | 828 m²/m³ | 630 m²/m³ | 828 m²/m³ |
| **Canopy Drag Factor (_D_)** | 8.28 s⁻¹ | 4.72 s⁻¹ | 8.28 s⁻¹ |
| **Characteristic Leaf Length (_L_)** | 20 mm | 8 mm | 20 mm |
| **Substrate Bulk / Void Volume** | 4.02 cm³ / 2.01 mL | 2.51 cm³ / 1.26 mL | 4.02 cm³ / 2.01 mL |
| **Upper Headspace Plenum** | Δ_Y_ = 7.7 mm | **Δ_Y_ = 14.2 mm (Open Plenum)** | Δ_Y_ = 7.7 mm |

---

## 3. Key Simulation Findings & Flight Operational Guidelines

### 3.1 "Dip-and-Converge" Aerodynamic Channelling & Age-of-Air Uniformity
- Forced air entering Port 1 (_Z_ = -21.6 mm) encounters tall foliage drag (_D_ = 8.28 s⁻¹) in Chamber 1, skimming across the upper manifold until reaching the step-down above the prostrate *Arabidopsis* rosette carpet in Chamber 2.
- The **Δ_Y_ = 14.2 mm open central plenum** directly beneath Port 2 acts as a low-resistance convergence corridor where lateral streams from Chambers 1 and 3 dip downward, flush the *Arabidopsis* lawn, and exhaust upward through Port 2.
- **Sandberg Conservation Verification**: At nominal flow _Q_ = 0.50 L/min, the numerical flux-weighted outlet mean age of air is θ_out = 7.27 s vs. theoretical τ = _V_fluid / _Q_ = 7.75 s (ratio = 0.938). Inter-chamber mean age uniformity is preserved within **< 1.0%** (Ch 1: 6.05 s, Ch 2: 6.10 s, Ch 3: 6.06 s).

### 3.2 Asymmetric Multi-Species Gas Exchange (Earth 1.0g vs. Lunar 0.166g)

| Regime | Pump Speed (_Q_) | Gravity (_g_) | *B. rapa* (Ch 1 & 3, 32 Seeds) _A_net | *B. rapa* Δ[CO₂] / Δ[H₂O] | *A. thaliana* (Ch 2, 100 Seeds) _A_net | *A. thaliana* Δ[CO₂] / Δ[H₂O] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Earth 1.0g Baseline** | 0.50 L/min | 9.81 m/s² | **13.04 μmol m⁻² s⁻¹** | -3.76 ppm / +7.11 mmol/m³ | **10.06 μmol m⁻² s⁻¹** | -1.31 ppm / +2.40 mmol/m³ |
| **Lunar 0.166g Uncompensated** | 0.50 L/min | 1.62 m/s² | **10.63 μmol m⁻² s⁻¹** (81.5%) | -3.77 ppm / +7.12 mmol/m³ | **8.43 μmol m⁻² s⁻¹** (83.8%) | -1.32 ppm / +2.42 mmol/m³ |
| **Lunar 0.166g Compensated** | **0.85 L/min** | 1.62 m/s² | **10.66 μmol m⁻² s⁻¹** | **-2.35 ppm / +4.50 mmol/m³** (-37%) | **8.44 μmol m⁻² s⁻¹** | **-0.85 ppm / +1.59 mmol/m³** (-37%) |

### 3.3 Diurnal Photoperiod Inversion & Lunar Nocturnal Pump Setback
- During the **8h dark skotoperiod**, LED arrays extinguish, leaves radiatively cool (Δ_T_ = -0.3 K), and mitochondrial dark respiration emits dense CO₂ (_S_CO2 = +6.0 μmol m⁻³ s⁻¹).
- Under Lunar 0.166g, low-power sleep setbacks (_Q_ ≤ 0.10 L/min) allow dense respiratory CO₂ to pool in the lower wells (Δ[CO₂] > +10 ppm).
- **Flight Rule**: Maintain a **nocturnal lunar pump floor of _Q_night ≥ 0.60 L/min** (Ri ≤ 0.05) to flush respiratory CO₂.

### 3.4 Harvest-Phase RNAlater Fixative Perfusion
- **Nominal Injection (_Q_inj = 10 mL/min)**: Root substrate saturation (5.28 mL void volume) completes in **_t_root = 31.7 s**; full chamber immersion (64.60 mL) completes in **_t_flood = 387.6 s (6.46 min)** with total pressure drop **Δ_P_ = 0.73 kPa**.
- **Rapid Injection (_Q_inj = 20 mL/min)**: Root saturation in **15.8 s**, full chamber flooding in **3.23 min**, with **Δ_P_ = 1.47 kPa** (well below seal safety margins).

### 3.5 14-Day Developmental Growth Trajectory (Days 0–14)
- **Imbibition Setback (Days 0–2)**: A gentle _Q_ = 0.25 L/min germination flow prevents seed-coat desiccation on the top surface of Substrate Blocks 1–3 (_Y_ = 14.62 mm) while clearing dark seed respiration (-0.45 to -0.72 nmol CO₂/s).
- **Autotrophic Transition**: Net daily carbon fixation turns positive at **Day 2.8** for *B. rapa* and **Day 3.2** for *A. thaliana*.
- **Adaptive Lunar Pump Schedule**: Stepping _Q_day from 0.25 → 0.50 → 0.65 → 0.85 L/min across Days 0–14 reduces *B. rapa* boundary-layer drawdown by **41%** (Δ[CO₂] from -3.98 ppm to -2.35 ppm at Day 14) while root biomass compaction (ε: 0.50 → 0.34) raises Day 14 RNAlater injection Δ_P_ from 0.73 kPa to 4.04 kPa (safely below the 5.0 kPa seal limit).

### 3.6 Off-Nominal Spacecraft Sensitivity & Fault-Tolerance Sweeps
- **Elevated Spacecraft Cabin CO₂ (420–5,000 ppm)**: Cabin CO₂ ≥ 2,000 ppm saturates FvCB carboxylation at the TPU plateau (_A_net ≈ 29.9–31.6 μmol m⁻² s⁻¹ in *B. rapa*), masking Lunar carbon limitation, **but** induces 50% stomatal closure (_g_s: 0.165 → 0.082 mol m⁻² s⁻¹ at 4,000 ppm), reducing evaporative cooling and elevating Lunar leaf temperature to Δ_T_leaf = +3.05 °C at 0.50 L/min vs. +2.19 °C at 0.85 L/min.
- **Root-Induced Substrate Compaction (ε = 0.50 → 0.20)**: If root density compresses substrate porosity below ε = 0.32, stepping RNAlater injection down from 10 mL/min to 5 mL/min maintains Δ_P_ < 3.5 kPa.
- **Pump Brownout Grace Window**: Complete pump outage (_Q_ ≤ 0.04 L/min) has a **grace window of _t_grace = 10.4–18.5 min** before transpired humidity reaches 92% RH condensation onset on the lid window.

### 3.7 Porous Substrate Capillarity & Optical Lid Condensation Defense
- **Lunar Capillary Retention (Richards / van Genuchten)**: In Lunar 0.166g, capillary forces (_P_c = 2σ cos θ_c / _r_pore) exceed gravity drainage by >250×, producing a uniform hydration profile (θ_w = 0.314 m³/m³) across the 7.0 mm substrate cassettes. Air-filled porosity is sustained at ε_air = 0.266 m³/m³ (above the hypoxia ceiling ε_crit = 0.15), with active drainage through the lower aeration wells (_Y_ = 0–7.62 mm).
- **Optical Lid Condensation Margin (Δ_T_margin)**: Transpiration across 164 plants injects 4.4 mL/day of moisture. Nominal flow (_Q_ = 0.50 L/min) maintains inner lid temperature above dew point (Δ_T_margin = +0.83 °C daytime, +2.15 °C at 0.85 L/min), preventing camera window fogging unless ventilation drops below _Q_crit ≈ 0.08 L/min.
- **Two-Phase RNAlater Venting**: Injection at 10 mL/min via Port 3 dipping nozzle (_Y_ = 9.5 mm) operates at Ca = 1.78 × 10⁻³ and Bo_lunar = 6.31 × 10⁻⁵ ≪ 1. Root substrate wets in 26.4 s, with 100% displaced air vented via Port 2 at Δ_P_vent < 0.01 kPa with zero trapped bubbles.

---

## 4. Simulation Figure Suite

| Figure | Description | Files |
| :--- | :--- | :--- |
| **Orthographic Blueprint** | 3-view CAD projection (Top Plan, Side Elevation, Longitudinal Cutaway) with full ABAI QC badges | [PNG](figures/fig_orthographic_projection_annotated.png) · [SVG](figures/fig_orthographic_projection_annotated.svg) |
| **Assembly Cutaways** | Mated STEP cross-sections across compartments (_X_ = -5 mm), ports (_X_ = -15.49 mm), and center (_Z_ = 0) | [Annotated](figures/assembly_cross_sections_annotated.png) · [Exact CAD](figures/assembly_cross_sections_exact.png) |
| **Figure 1** | Parametric pump sweep (0.05–4.0 L/min): velocity, boundary-layer conductance, Δ[CO₂], _A_net | [PNG](figures/fig1_pump_speed_vs_aerodynamics.png) · [SVG](figures/fig1_pump_speed_vs_aerodynamics.svg) |
| **Figure 2** | 2D/3D velocity vector fields & scalar distributions comparing Earth 1.0g vs. Lunar 0.166g | [PNG](figures/fig2_velocity_fields_earth_vs_lunar.png) · [SVG](figures/fig2_velocity_fields_earth_vs_lunar.svg) |
| **Figure 3** | Lunar compensatory pump sizing curve & regime classification (Ri, Pe, Gr/Re²) | [PNG](figures/fig3_lunar_compensatory_pump_curve.png) · [SVG](figures/fig3_lunar_compensatory_pump_curve.svg) |
| **Figure 4** | 3D spatial Mean Age of Air (θ), chamber uniformity, and Sandberg identity verification | [PNG](figures/fig4_mean_age_of_air_and_chamber_uniformity.png) · [SVG](figures/fig4_mean_age_of_air_and_chamber_uniformity.svg) |
| **Figure 5** | Flight hardware benchmark comparing LEAF GBC2 against NASA BRIC, CARA, and Veggie | [PNG](figures/fig5_flight_hardware_benchmark_gbc2_vs_nasa.png) · [SVG](figures/fig5_flight_hardware_benchmark_gbc2_vs_nasa.svg) |
| **Figure 6** | Predicted transcriptomic stress & phytohormone signaling phenotypes (Auxin, ABA, Ethylene) | [PNG](figures/fig6_predicted_hormone_carbon_phenotypes.png) · [SVG](figures/fig6_predicted_hormone_carbon_phenotypes.svg) |
| **Figure 7** | Transient RNAlater perfusion front progression & multi-species asymmetric gas exchange | [PNG](figures/fig7_rnalater_perfusion_and_species_dynamics.png) · [SVG](figures/fig7_rnalater_perfusion_and_species_dynamics.svg) |
| **Figure 8** | Diurnal day/night photoperiod inversion, respiratory CO₂ pooling, and nocturnal setback envelope | [PNG](figures/fig8_diurnal_day_night_photoperiod_inversion.png) · [SVG](figures/fig8_diurnal_day_night_photoperiod_inversion.svg) |
| **Figure 9** | 14-day developmental growth time-series: logistic canopy ontogeny, channelling, cumulative C, & root compaction | [PNG](figures/fig9_developmental_growth_timeseries.png) · [SVG](figures/fig9_developmental_growth_timeseries.svg) |
| **Figure 10** | Off-nominal spacecraft sensitivity sweeps: cabin CO₂ (420–5,000 ppm), root compaction, & pump brownout | [PNG](figures/fig10_off_nominal_sensitivity_sweeps.png) · [SVG](figures/fig10_off_nominal_sensitivity_sweeps.svg) |
| **Figure 11** | Porous substrate unsaturated capillarity, evapotranspiration, optical lid condensation defense, & two-phase venting | [PNG](figures/fig11_water_capillarity_evapotranspiration_condensation.png) · [SVG](figures/fig11_water_capillarity_evapotranspiration_condensation.svg) |

---

## 5. Repository Structure & Reproducibility

```text
LEAF_GBC2_pot_model/
├── index.html                           # GitHub Pages Interactive Research Summary Portal (CoSE Theme)
├── viewer_3d.html                       # Interactive 3D WebGL Digital Twin (Three.js & OrbitControls)
├── MANUSCRIPT_LEAF_GBC2_CFD.md          # Complete Technical Manuscript (Methods, Results, Legends)
├── README.md                            # Project documentation & summary
├── requirements.txt                     # Python environment dependencies
├── cad/
│   ├── pot.step                         # Native STEP CAD for GBC2 Pot Body (M09024F-MDL-5)
│   └── lid.step                         # Native STEP CAD for GBC2 Lid Assembly (M09025F-MDL-2)
├── geometry/
│   ├── gbc2_assembly.stl                # Mated 3D assembly mesh
│   └── gbc2_grid.npz                    # 3D voxelized domain (60x48x80) with species & substrate tags
├── sim/
│   ├── build_geometry.py                # CAD mating & 3D voxelization pipeline
│   ├── gbc2_cfd_engine.py               # 3D Navier-Stokes, Brinkman-Darcy & FvCB solver
│   ├── advanced_analytics.py            # Mean Age of Air, Sandberg conservation & hormone indices
│   ├── generate_orthographic_projections.py  # 3-view orthographic projection generator
│   ├── generate_assembly_cross_sections.py   # Mated CAD cross-section generator
│   ├── run_parametric_sweep.py          # Parametric pump speed sweep (Figs 1-3)
│   ├── run_advanced_pipeline.py         # Age of air, hardware benchmark & phenotypes (Figs 4-6)
│   ├── simulate_rnalater_perfusion.py   # RNAlater fixation & species asymmetry (Fig 7)
│   ├── simulate_diurnal_cycle.py        # Day/night photoperiod inversion (Fig 8)
│   ├── simulate_developmental_growth_timeseries.py # 14-day developmental growth trajectory (Fig 9)
│   ├── simulate_off_nominal_sensitivity.py         # Off-nominal mission sensitivity sweeps (Fig 10)
│   └── simulate_water_capillarity_and_humidity.py  # Capillary retention & lid condensation (Fig 11)
├── figures/                             # All publication figures (.png and .svg)
└── results/
    └── tables/                          # All quantitative simulation tables (.tsv)
```

### Running the Full Simulation & Visualization Pipeline

```bash
# 1. Create virtual environment and install dependencies
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# 2. Rebuild 3D voxel geometry from STEP CAD files
.venv/bin/python sim/build_geometry.py

# 3. Generate Orthographic Projections & CAD Cross-Sections
.venv/bin/python sim/generate_orthographic_projections.py
.venv/bin/python sim/generate_assembly_cross_sections.py

# 4. Run CFD & Biophysics Simulation Suites (Figures 1-10 & TSV tables)
.venv/bin/python sim/run_parametric_sweep.py
.venv/bin/python sim/run_advanced_pipeline.py
.venv/bin/python sim/simulate_rnalater_perfusion.py
.venv/bin/python sim/simulate_diurnal_cycle.py
.venv/bin/python sim/simulate_developmental_growth_timeseries.py
.venv/bin/python sim/simulate_off_nominal_sensitivity.py
```

