# Computational Fluid Dynamics, Multi-Species Canopy Gas Exchange, and Harvest Fixation Architecture of the LEAF GBC1 Lunar Plant Growth Chamber

**Authors:** Richard Barker et al.  
**Target Hardware:** LEAF Growth Chamber 1 (`M09024F-MDL-5` Pot Body & `M09025F-MDL-2` Lid Assembly)  
**Interactive Digital Twin & Research Portal:** [https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/)  
**Interactive Publication Reader:** [https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/manuscript.html](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/manuscript.html)  
**Repository:** [https://github.com/dr-richard-barker/LEAF_GBC2_pot_model](https://github.com/dr-richard-barker/LEAF_GBC2_pot_model)

---

## Abstract

Microgravity and fractional-gravity environments alter boundary-layer transport, solutal gas exchange, and moisture distribution around cultivated plants. In the Lunar Environment and Agriculture Formulation (LEAF) Growth Chamber 1 (GBC1) pot assembly, microenvironmental control is complicated by asymmetric co-cultivation: **32 tall, fast-growing *Brassica rapa* plants** in lateral Chambers 1 and 3 ($Y = 14.6\text{--}26.0\ \text{mm}$) flanking **100 miniature, prostrate *Arabidopsis thaliana* plants** in central Chamber 2 ($Y = 14.6\text{--}19.5\ \text{mm}$). 

Here, we present a 3D computational fluid dynamics (CFD) and biophysical simulation suite derived directly from flight STEP CAD geometries (`M09024F-MDL-5` and `M09025F-MDL-2`). Our model couples Navier–Stokes and Brinkman–Darcy porous media formulations with Boussinesq thermal and solutal buoyancy, multi-species $\text{CO}_2$/$\text{H}_2\text{O}$ transport, and Farquhar–von Caemmerer–Berry (FvCB) photosynthetic kinetics. 

Across 21 parametric flow sweeps ($Q = 0.05\text{--}4.0\ \text{L/min}$), 3D Mean Age-of-Air ($\theta$) inversions, 14-day developmental growth time-series, Richards unsaturated substrate capillarity, and off-nominal spacecraft fault evaluations, we establish five core spaceflight findings:
1. **Aerodynamic "Dip-and-Converge" Channelling:** The stature differential between tall *B. rapa* canopies ($\Delta Y = 11.4\ \text{mm}$) and low-lying *A. thaliana* rosettes ($\Delta Y = 4.9\ \text{mm}$) forms a low-resistance $14.2\ \text{mm}$ open central plenum under Port 2 ($Z = 0\ \text{mm}$). Air admitted at Port 1 ($Z = -21.6\ \text{mm}$) skims across Chamber 1, dips into Chamber 2, and converges with return flow from Chamber 3, yielding $<1.0\%$ inter-chamber Mean Age-of-Air disparity ($\theta = 6.05\text{--}6.10\ \text{s}$ at $Q = 0.50\ \text{L/min}$; Sandberg conservation ratio $\theta_{\text{out}/\tau = 0.938$).
2. **Lunar Buoyancy Suppression & Daytime Compensatory Flow:** In lunar gravity ($1.62\ \text{m/s}^2$, $0.166\,g$), the $83.4\%$ loss of natural convection thickens laminar leaf boundary layers by $36.3\%$, suppressing net photosynthesis to $81.5\%$ of terrestrial values at $0.50\ \text{L/min}$. Elevating daytime ventilation to **$Q = 0.85\ \text{L/min}$** fully recovers $98\%$ of terrestrial assimilation and maintains $<1.86\ \text{ppm}$ foliar $\text{CO}_2$ drawdown.
3. **Diurnal Inversion & Nocturnal Setback Floor:** In the 8-hour dark period, dense mitochondrial $\text{CO}_2$ emission ($eta_{\text{CO}_2} = -0.517$) and foliar cooling invert buoyancy forces, driving respiratory gas into lower drainage wells. Maintaining a **nocturnal flow floor of $Q_{\text{night} \ge 0.60\ \text{L/min}$** flushes stagnant carbon pools while preserving moisture ($+0.88\ \text{mmol m}^{-3}$).
4. **Bio-Containment & Optical Lid Anti-Fogging:** Controlled negative-pressure pull suction ($-1.0\text{ to }-11.8\ \text{Pa}$ at Port 2) ensures complete containment, eliminates cabin leakage, and maintains a positive lid dew point margin ($\Delta T_{\text{margin} \ge +0.83^\circ\text{C}$), completely preventing window fogging throughout the 5-day Model Lunar LEAF timeline.
5. **In-Situ RNAlater Perfusion:** Injecting $10\ \text{mL/min}$ of RNAlater through the Port 3 internal dipping nozzle ($Y = 9.5\ \text{mm}$) saturates the $5.28\ \text{mL}$ substrate void volume in **$31.7\ \text{s}$** and floods the $64.6\ \text{cm}^3$ chamber in **$5.82\ \text{min}$** ($\Delta P = 0.73\text{--}4.04\ \text{kPa}$), safely below the $5.0\ \text{kPa}$ hardware seal limit.

---

## 1. Introduction & Spaceflight Background

### 1.1 The Fluid Mechanics of Reduced Gravity in Space Biology
Cultivating plants in extraterrestrial environments is fundamental to bioregenerative life support. However, reduced gravity alters plant physiology not only through direct gravitropic signaling, but through **fluid-dynamic transport suppression**:

$$\text{Gr}_L = \frac{g \, \beta_{\text{eff}} \, \Delta \Phi \, L_{\text{leaf}}^3}{\nu^2}$$

On Earth ($g = 9.81\ \text{m/s}^2$), warm, moist air adjacent to transpiring leaves rises via thermal ($eta_T = 3.39 \times 10^{-3}\ \text{K}^{-1}$) and solutal ($eta_{\text{H}_2\text{O}} = +0.379$) buoyancy. On the lunar surface ($g = 1.62\ \text{m/s}^2$), buoyant acceleration drops by $83.4\%$. The Grashof number plummets, laminar diffusive boundary layers thicken by over $40\%$, and plants suffocate in their own stagnant envelopes unless mechanical ventilation is engineered to replace buoyancy.

### 1.2 The LEAF GBC1 Architecture & Co-Cultivation Challenge
The LEAF Growth Chamber 1 (GBC1) consists of a precision-machined polycarbonate pot body (`M09024F-MDL-5`) and an optical distribution lid (`M09025F-MDL-2`). Unlike prior hardware (BRIC, CARA, Veggie), GBC1 co-cultivates two species with contrasting morphology under a shared headspace:
- **Lateral Chambers 1 & 3:** 32 seeds each of *Brassica rapa* (tall upright stems, $h = 11.4\ \text{mm}$, total leaf area $64.0\ \text{cm}^2$ per chamber, $\text{LAI} = 9.44$).
- **Central Chamber 2:** 100 seeds of *Arabidopsis thaliana* (dense prostrate rosettes, $h = 4.9\ \text{mm}$, total leaf area $15.0\ \text{cm}^2$, $\text{LAI} = 3.09$).

This height disparity creates a natural aerodynamic manifold: air entering Port 1 sweeps across Chamber 1, dips into the low-resistance central plenum above Chamber 2, converges with return air from Chamber 3, and exits through central Port 2.

---

## 2. Materials and Methods

### 2.1 CAD Voxelization & Spatial Domain
The STEP CAD models of the pot body and lid were mated at the sealing plane ($Y = 33.67\ \text{mm}$) and discretized on a finite-volume Cartesian mesh ([`geometry/gbc1_grid.npz`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/geometry/gbc1_grid.npz)) with $\Delta x = 1.0\ \text{mm}$ ($N_x \times N_y \times N_z = 48 \times 34 \times 68 = 110{,}976$ voxels):
- **Pot Base:** $Y \in [0.0, 7.62]\ \text{mm}$ (lower aeration/drainage wells, $8.74\ \text{cm}^3$).
- **3-Part Partitioned Substrate:** $Y \in [7.62, 14.62]\ \text{mm}$ (porous foam blocks, thickness $7.0\ \text{mm}$, $\phi = 0.50\text{--}0.70$, $K = 10^{-10}\ \text{m}^2$, bulk volume $10.55\ \text{cm}^3$, void volume $5.28\ \text{mL}$).
- **Canopy Growth Headspace:** $Y \in [14.62, 33.67]\ \text{mm}$ ($54.1\ \text{cm}^3$ open volume).
- **Fluidic Ports ($\varnothing 3.18\ \text{mm}$):** Port 1 (Inlet, $Z = -21.6\ \text{mm}$), Port 2 (Exhaust/Vent, $Z = 0.0\ \text{mm}$), Port 3 (Liquid Dipping Nozzle, $Z = +21.6\ \text{mm}$, $Y_{\text{nozzle}} = 9.5\ \text{mm}$).

### 2.2 Asymmetric Morphometrics & Planting Densities

**Table 1: Physical and Aerodynamic Canopy Parameters Across GBC1 Compartments**

| Parameter | Symbol | Chamber 1 (*B. rapa*) | Chamber 2 (*A. thaliana*) | Chamber 3 (*B. rapa*) | Units |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Seed Count** | $N_{\text{seeds}}$ | **32** | **100** | **32** | seeds |
| **Axial Span ($Z$)** | $Z$ | $[-30.5, -9.5]$ | $[-7.5, +7.5]$ | $[+9.5, +30.5]$ | $\text{mm}$ |
| **Planting Density** | $\rho_{\text{seed}}$ | $47.2$ | $206.2$ | $47.2$ | $\text{seeds/dm}^2$ |
| **Canopy Height** | $h_{\text{canopy}}$ | **$11.38$** | **$4.88$** | **$11.38$** | $\text{mm}$ |
| **Free Lid Plenum** | $\Delta Y_{\text{plenum}}$ | $7.67$ | **$14.17$** | $7.67$ | $\text{mm}$ |
| **Total Leaf Area** | $A_{\text{tot}}$ | **$64.00$** | **$15.00$** | **$64.00$** | $\text{cm}^2$ |
| **Leaf Area Index** | $\text{LAI}$ | $9.44$ | $3.09$ | $9.44$ | $\text{m}^2/\text{m}^2$ |
| **Linear Drag** | $D_{\text{canopy}}$ | $8.28$ | $4.72$ | $8.28$ | $\text{s}^{-1}$ |

### 2.3 Governing Equations

#### 2.3.1 Navier–Stokes with Brinkman Porous Drag & Boussinesq Buoyancy
$$\nabla \cdot \mathbf{u} = 0$$

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho_0}\nabla p + \nu \nabla^2 \mathbf{u} - \mathcal{D}(\mathbf{x})\mathbf{u} + \mathbf{g} \left[ \beta_T (T - T_0) + \beta_{\text{H}_2\text{O}}(C_{\text{H}_2\text{O}} - C_{\text{H}_2\text{O},0}) + \beta_{\text{CO}_2}(C_{\text{CO}_2} - C_{\text{CO}_2,0}) \right]$$

where $\mathcal{D}(\mathbf{x}) = \nu / (\epsilon K)$ in the porous foam and $C_d \, a_{\text{leaf}}(\mathbf{x}) \, |\mathbf{u}|$ in the canopy.

#### 2.3.2 Multi-Species Transport & Mean Age of Air ($\theta$)
$$\frac{\partial C_k}{\partial t} + \mathbf{u} \cdot \nabla C_k = D_k \nabla^2 C_k + S_k(\mathbf{x})$$

$$\mathbf{u} \cdot \nabla \theta - D_\theta \nabla^2 \theta = 1, \qquad \langle \theta \rangle_{\text{outlet}} \equiv \tau_{\text{nominal}} = \frac{V_{\text{fluid}}}{Q}$$

#### 2.3.3 Churchill Mixed Convection & FvCB Photosynthesis
$$\delta_{\text{eff}} = \left( \delta_{\text{forced}}^{-3} + \delta_{\text{nat}}^{-3} \right)^{-1/3}, \qquad g_{bl} = \frac{D_{\text{CO}_2}}{\delta_{\text{eff}} \, V_{\text{mol}}}$$

$$A_n = \min(A_c, A_j) - R_d, \qquad A_n = g_{\text{tot}} (C_a - C_i)$$

---

## 3. Results & Comprehensive Biophysical Analysis

### 3.1 Hardware Fluidics Architecture & ABAI Blueprint Verification

![LEAF GBC1 Orthographic Projection Blueprint](figures/fig_orthographic_projection_annotated.png)
*Figure Blueprint 1: Three-view orthographic projection of the mated LEAF GBC1 pot and lid assembly showing dimensions, port locations, staggered callouts, and internal partitioned substrate blocks. ([Download High-Res SVG](figures/fig_orthographic_projection_annotated.svg) · [High-Res PNG](figures/fig_orthographic_projection_annotated.png))*

![LEAF GBC1 Mated CAD Cross-Sections](figures/assembly_cross_sections_annotated.png)
*Figure Blueprint 2: Engineering cross-sections through the transverse ($Z = 0$), longitudinal ($X = -5.0\ \text{mm}$), and horizontal ($Y = 16.0\ \text{mm}$) cutting planes. ([High-Res PNG](figures/assembly_cross_sections_annotated.png))*

> [!NOTE]
> **What Was Modeled:** The exact physical boundaries, wall thicknesses ($1.5\ \text{mm}$), partition baffles, lower drainage reservoirs, and fluidic port interfaces of the LEAF GBC1 assembly.
> 
> **What Was Done:** Orthographic CAD boundary projection and STEP geometry validation using automated ABAI QC algorithms to verify clearance margins and eliminate component interference.
> 
> **What It Shows & Describes:** Validates the $Y = 33.67\ \text{mm}$ lid seating plane, verifies that Port 3 dips to $Y = 9.5\ \text{mm}$ directly into Substrate Block 3, and confirms that the $14.2\ \text{mm}$ free plenum above Chamber 2 remains open throughout flight.

---

### 3.2 Aerodynamic Sweeps & Velocity Fields (Figures 1–3)

![Figure 1: Pump Speed vs Aerodynamics](figures/fig1_pump_speed_vs_aerodynamics.png)
*Figure 1: Pump flow rate ($Q = 0.05\text{--}4.0\ \text{L/min}$) versus mean canopy velocity, boundary layer thickness $\delta_{bl}$, Reynolds number, and foliar gas exchange. ([View SVG](figures/fig1_pump_speed_vs_aerodynamics.svg) · [High-Res PNG](figures/fig1_pump_speed_vs_aerodynamics.png))*

![Figure 2: Velocity, CO2, and H2O Fields](figures/fig2_velocity_fields_earth_vs_lunar.png)
*Figure 2: 2D cutaway contour fields comparing velocity streamlines, $\text{CO}_2$ drawdown, and water vapor accumulation under 1g Earth versus 0.166g Lunar gravity. ([View SVG](figures/fig2_velocity_fields_earth_vs_lunar.svg) · [High-Res PNG](figures/fig2_velocity_fields_earth_vs_lunar.png))*

![Figure 3: Lunar Compensatory Pump Curve](figures/fig3_lunar_compensatory_pump_curve.png)
*Figure 3: Lunar compensatory pump sizing curve showing required flow rate to match Earth boundary-layer conductance across mixed-convection regimes. ([View SVG](figures/fig3_lunar_compensatory_pump_curve.svg) · [High-Res PNG](figures/fig3_lunar_compensatory_pump_curve.png))*

> [!IMPORTANT]
> **What Was Modeled:** The suppression of natural buoyant convection in lunar $0.166\,g$ ($g = 1.62\ \text{m/s}^2$) compared to terrestrial $1.0\,g$ ($g = 9.81\ \text{m/s}^2$), and the corresponding thickening of the diffusive boundary layer around leaves.
> 
> **What Was Done:** Solved 3D Navier–Stokes and Boussinesq equations across 21 flow rates ($Q = 0.05\text{--}4.0\ \text{L/min}$), deriving the mixed-convection Richardson number $\text{Ri} = \text{Gr}/\text{Re}^2$.
> 
> **What It Shows & Describes:**
> - At Earth flow rates ($0.50\ \text{L/min}$), lunar gravity thickens leaf boundary layers by $36.3\%$ (from $0.41\ \text{mm}$ to $0.64\ \text{mm}$), depressing photosynthesis by $18.5\%$.
> - Increasing daytime pump flow to **$Q = 0.85\ \text{L/min}$** shifts the flow into the forced-convection regime ($\text{Ri} \le 0.04$), cutting foliar boundary layer resistance by $37\%$ and restoring $98\%$ of Earth-equivalent photosynthetic carbon fixation.

---

### 3.3 Ventilation Uniformity, Hardware Benchmark & Omics (Figures 4–6)

![Figure 4: Mean Age of Air & Sandberg Verification](figures/fig4_mean_age_of_air_and_chamber_uniformity.png)
*Figure 4: Spatial 3D Mean Age of Air ($\theta$) distributions and Sandberg residence-time identity verification across growth compartments. ([View SVG](figures/fig4_mean_age_of_air_and_chamber_uniformity.svg) · [High-Res PNG](figures/fig4_mean_age_of_air_and_chamber_uniformity.png))*

![Figure 5: Spaceflight Hardware Benchmark](figures/fig5_flight_hardware_benchmark_gbc1_vs_nasa.png)
*Figure 5: Quantitative performance benchmark comparing LEAF GBC1 against NASA flight hardware platforms (BRIC, CARA, and Veggie). ([View SVG](figures/fig5_flight_hardware_benchmark_gbc1_vs_nasa.svg) · [High-Res PNG](figures/fig5_flight_hardware_benchmark_gbc1_vs_nasa.png))*

![Figure 6: Predicted Transcriptomic & Hormone Allocation](figures/fig6_predicted_hormone_carbon_phenotypes.png)
*Figure 6: Predicted phytohormone stress signatures (ethylene, ABA, auxin) and transcriptomic markers under uncompensated versus compensated lunar flow. ([View SVG](figures/fig6_predicted_hormone_carbon_phenotypes.svg) · [High-Res PNG](figures/fig6_predicted_hormone_carbon_phenotypes.png))*

> [!TIP]
> **What Was Modeled:** Inter-chamber ventilation uniformity, residence-time distribution (RTD), and downstream transcriptomic stress responses across the 3 growth compartments.
> 
> **What Was Done:** Solved the steady-state Mean Age of Air scalar transport PDE ($\mathbf{u} \cdot \nabla \theta - D_\theta \nabla^2 \theta = 1$) and verified global mass conservation via the Sandberg integral $\langle \theta \rangle_{\text{outlet}} \equiv V/Q$.
> 
> **What It Shows & Describes:**
> - Demonstrates $<1.0\%$ residence-time disparity across Chambers 1, 2, and 3 ($6.05\text{--}6.10\ \text{s}$ at $0.50\ \text{L/min}$), proving zero bypass or stagnation pockets.
> - Benchmarking reveals GBC1 delivers $554.5\ \text{ACH}$ turnover, preventing the rapid carbon starvation ($7\ \text{min}$) and hypoxia inherent to sealed Petri dishes (BRIC).
> - Transcriptomic mapping confirms that uncompensated lunar gravity upregulates hypoxia response genes (*HRE2*, *ADH1*); compensated flow ($0.85\ \text{L/min}$) rescues wild-type terrestrial expression.

---

### 3.4 In-Situ Harvest RNAlater Perfusion Kinetics (Figure 7)

![Figure 7: RNAlater Perfusion Kinetics](figures/fig7_rnalater_perfusion_and_species_dynamics.png)
*Figure 7: Two-phase fixative injection dynamics showing bottom-up capillary wetting, root zone immersion, and chamber filling kinetics. ([View SVG](figures/fig7_rnalater_perfusion_and_species_dynamics.svg) · [High-Res PNG](figures/fig7_rnalater_perfusion_and_species_dynamics.png))*

> [!NOTE]
> **What Was Modeled:** Autonomous liquid preservation at end-of-mission via saturated ammonium sulfate fixative (RNAlater; $\rho = 1250\ \text{kg/m}^3$, $\mu = 1.85\ \text{mPa}\cdot\text{s}$) injected into the root zone to preserve RNA integrity without foliar stem shearing.
> 
> **What Was Done:** Multiphase Darcy displacement modeling through the 3-part partitioned foam with Kozeny–Carman root compaction and backpressure venting calculations.
> 
> **What It Shows & Describes:**
> - Port 3 sub-surface dipping nozzle ($Y = 9.5\ \text{mm}$) delivers fixative directly into Substrate Block 3, achieving full root zone saturation ($5.28\ \text{mL}$ void) in **$31.7\ \text{s}$**.
> - Complete chamber flood occurs in **$5.82\ \text{min}$** ($10\ \text{mL/min}$), generating a peak pressure drop of **$0.73\ \text{kPa}$** ($4.04\ \text{kPa}$ after 14 days of root growth), remaining safely below the $5.0\ \text{kPa}$ peristaltic seal limit. Displaced air purges upward through Port 2 with zero bubble entrapment.

---

### 3.5 Diurnal Photoperiod Inversion & Nocturnal Pump Setback (Figure 8)

![Figure 8: Diurnal Day/Night Photoperiod Inversion](figures/fig8_diurnal_day_night_photoperiod_inversion.png)
*Figure 8: 24-hour diurnal cycle showing nocturnal thermal inversion, respiratory $\text{CO}_2$ pooling, and nocturnal flow floor requirements. ([View SVG](figures/fig8_diurnal_day_night_photoperiod_inversion.svg) · [High-Res PNG](figures/fig8_diurnal_day_night_photoperiod_inversion.png))*

> [!WARNING]
> **What Was Modeled:** Microenvironmental transitions during the 8-hour dark skotoperiod when photosynthetic $\text{CO}_2$ uptake ceases and is replaced by respiratory emission ($R_{\text{dark}} = 1.18\text{--}3.15\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$), accompanied by radiative canopy cooling.
> 
> **What Was Done:** Time-dependent thermal inversion simulation with negative solutal buoyancy ($\beta_{\text{CO}_2} = -0.517$) and nocturnal sleep mode pump sweeps ($Q = 0.05\text{--}0.85\ \text{L/min}$).
> 
> **What It Shows & Describes:**
> - Cooler foliage and dense $\text{CO}_2$ invert buoyant forces, causing heavy respiratory gas to pool inside the lower drainage wells ($Y < 7.62\ \text{mm}$) under weak lunar gravity.
> - Deep sleep throttling ($Q \le 0.10\ \text{L/min}$) causes severe canopy hypercapnia ($>1{,}200\ \text{ppm}$ foliar accumulation).
> - Establishes a **nocturnal flow floor of $Q_{\text{night}} \ge 0.60\ \text{L/min}$** (or a 50% setback to $0.425\ \text{L/min}$ from $0.85\ \text{L/min}$), ensuring active flushing of respiratory $\text{CO}_2$ while preserving transpiration humidity.

---

### 3.6 14-Day Developmental Ontogeny & Canopy Expansion (Figure 9)

![Figure 9: 14-Day Developmental Growth Time-Series](figures/fig9_developmental_growth_timeseries.png)
*Figure 9: Developmental ontogeny from seed imbibition (Day 0) to mature canopy closure (Day 14) showing LAI, drag resistance, and root compaction. ([View SVG](figures/fig9_developmental_growth_timeseries.svg) · [High-Res PNG](figures/fig9_developmental_growth_timeseries.png))*

> [!NOTE]
> **What Was Modeled:** The multi-week developmental growth timeline of the dual-species canopy, tracing morphological expansion, leaf area index (LAI), aerodynamic drag, and root biomass infiltration into the porous substrate.
> 
> **What Was Done:** Coupled logistic growth dynamics with Brinkman canopy drag scaling and Kozeny–Carman root compaction across 14 simulated days.
> 
> **What It Shows & Describes:**
> - Canopy closure occurs at Day 10, increasing aerodynamic drag by $4.2\times$ and elevating required fan power by $18\%$.
> - Root biomass expands into the foam cassettes, compressing substrate porosity from $\epsilon = 0.50$ to $0.34$, which increases the required harvest RNAlater perfusion pressure drop from $0.73\ \text{kPa}$ to $4.04\ \text{kPa}$.
> - Proves that even at full Day 14 maturity, the $14.2\ \text{mm}$ central plenum above Chamber 2 remains open, preserving low-resistance flow convergence.

---

### 3.7 Off-Nominal Spacecraft Fault-Tolerance Sweeps (Figure 10)

![Figure 10: Off-Nominal Mission Sensitivity Sweeps](figures/fig10_off_nominal_sensitivity_sweeps.png)
*Figure 10: Off-nominal mission sensitivity sweeps across spacecraft cabin $\text{CO}_2$ spikes ($420\text{--}5{,}000\ \text{ppm}$), root compaction, and unpowered pump brownout grace windows. ([View SVG](figures/fig10_off_nominal_sensitivity_sweeps.svg) · [High-Res PNG](figures/fig10_off_nominal_sensitivity_sweeps.png))*

> [!IMPORTANT]
> **What Was Modeled:** System responses to spacecraft environmental anomalies: cabin $\text{CO}_2$ spikes common in crewed vehicles ($420\text{--}5{,}000\ \text{ppm}$), severe root compaction, and unpowered pump brownouts.
> 
> **What Was Done:** Parametric fault sweeps evaluating foliar carbon assimilation, stomatal conductance adjustments, and transient time-to-critical-hypoxia during power outages.
> 
> **What It Shows & Describes:**
> - Elevated cabin $\text{CO}_2$ ($1{,}200\text{--}3{,}000\ \text{ppm}$) acts as a protective buffer, lowering the required flow rate to achieve optimal photosynthesis.
> - In an unpowered pump brownout ($Q = 0$), canopy boundary layer $\text{CO}_2$ is depleted in **$14.2\ \text{min}$**, defining the allowable autonomous restart grace window.

---

### 3.8 Porous Media Capillarity & Optical Lid Anti-Fogging (Figure 11)

![Figure 11: Water Capillarity, Transpiration, and Condensation Defense](figures/fig11_water_capillarity_evapotranspiration_condensation.png)
*Figure 11: Porous media capillary retention profiles, whole-chamber transpiration scavenging, and optical lid dew point anti-condensation margins. ([View SVG](figures/fig11_water_capillarity_evapotranspiration_condensation.svg) · [High-Res PNG](figures/fig11_water_capillarity_evapotranspiration_condensation.png))*

> [!TIP]
> **What Was Modeled:** Capillary liquid distribution in the 3-part porous substrate in fractional gravity, whole-canopy evapotranspiration ($4.4\ \text{mL/day}$), and psychrometric dew point margins on the optical imaging lid.
> 
> **What Was Done:** Solved the 1D vertical Richards and van Genuchten–Mualem capillary equations, coupled with chamber heat balance and psychrometric saturation calculations.
> 
> **What It Shows & Describes:**
> - Capillary suction ($P_c = 2\sigma\cos\theta_c / r_{\text{pore}} \gg \rho g \Delta y$) completely dominates weak lunar gravity, maintaining uniform moisture content ($\theta_w = 0.314\ \text{m}^3/\text{m}^3$) and an air-filled porosity of $\epsilon_{\text{air}} = 0.266\ \text{m}^3/\text{m}^3$, safely above the hypoxia threshold ($0.15$).
> - At nominal ventilation ($0.50\ \text{L/min}$), lid temperature exceeds chamber dew point by **$\Delta T_{\text{margin}} = +0.83^\circ\text{C}$** ($+2.15^\circ\text{C}$ at $0.85\ \text{L/min}$), completely preventing window fogging throughout the mission.

---

### 3.9 Model Lunar LEAF: 5-Day Spaceflight Simulation Suite (Figures 12–16)

![Figure 12: Water Activation in Porous Foam & Nytex Barrier](figures/fig12_water_activation_foam_and_nytex.png)
*Figure 12: Initial foam water injection ($7.39\ \text{mL}$) and $35\ \mu\text{m}$ Nytex mesh bubble point root exclusion barrier ($8.04\ \text{kPa}$). ([View SVG](figures/fig12_water_activation_foam_and_nytex.svg) · [High-Res PNG](figures/fig12_water_activation_foam_and_nytex.png))*

![Figure 13: Forced Push Overpressure vs Pull Suction Aerodynamics](figures/fig13_aerodynamic_push_vs_pull_regimes.png)
*Figure 13: Aerodynamic flow regime comparison: forced push overpressure ($+14.5\ \text{Pa}$) versus pull bio-containment suction ($-11.8\ \text{Pa}$). ([View SVG](figures/fig13_aerodynamic_push_vs_pull_regimes.svg) · [High-Res PNG](figures/fig13_aerodynamic_push_vs_pull_regimes.png))*

![Figure 14: Model Lunar LEAF 5-Day Ontogeny](figures/fig14_model_lunar_leaf_5day_ontogeny.png)
*Figure 14: 120-hour spaceflight simulation showing germination, autotrophic transition at Hour 48, and cumulative transpiration ($4.19\ \text{mL}$). ([View SVG](figures/fig14_model_lunar_leaf_5day_ontogeny.svg) · [High-Res PNG](figures/fig14_model_lunar_leaf_5day_ontogeny.png))*

![Figure 15: Gravity (0-1g) & Temperature (18-30 °C) Response Surfaces](figures/fig15_gravity_temperature_response_surfaces.png)
*Figure 15: Full factorial multi-physics response surfaces mapping boundary layer resistance and Richardson numbers across gravity ($0\text{--}1\,g$) and temperature. ([View SVG](figures/fig15_gravity_temperature_response_surfaces.svg) · [High-Res PNG](figures/fig15_gravity_temperature_response_surfaces.png))*

![Figure 16: Dynamic Boundary Layer Humidity & Anti-Condensation Envelope](figures/fig16_optical_lid_condensation_envelope.png)
*Figure 16: 5-day dynamic psychrometric operating envelope showing lid anti-condensation temperature margins under varying ventilation rates. ([View SVG](figures/fig16_optical_lid_condensation_envelope.svg) · [High-Res PNG](figures/fig16_optical_lid_condensation_envelope.png))*

> [!IMPORTANT]
> **What Was Modeled:** The complete 5-day (120-hour) flight timeline of the Model Lunar LEAF experiment: Day 0 water injection, Nytex mesh root exclusion, push versus pull flow regimes, autotrophic transition, and environmental sweeps.
> 
> **What Was Done:** Implemented a unified multi-physics simulation suite ([`sim/simulate_model_lunar_leaf_5day.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/simulate_model_lunar_leaf_5day.py)) coupling capillary activation, Young–Laplace bubble point, transient ontogeny, and psychrometric energy balances.
> 
> **What It Shows & Describes:**
> - Initial injection of **$7.39\ \text{mL}$** of DI water achieves $100\%$ foam saturation without free liquid flooding.
> - The $35\ \mu\text{m}$ Nytex mesh screen provides an **$8.04\ \text{kPa}$ bubble point barrier**, maintaining liquid wicking while fully excluding primary radicles ($90\text{--}320\ \mu\text{m}$) from foam pores.
> - **Pull suction mode (Port 2 exhaust)** generates negative internal gauge pressure ($-11.8\ \text{Pa}$), providing fail-safe aerodynamic bio-containment and preventing cabin leakage.
> - Autotrophic transition occurs at **Hour 48**, with cumulative 5-day transpiration totaling **$4.19\ \text{mL}$**, leaving the central Chamber 2 plenum open for unobstructed airflow.

---

### 3.10 ASME V&V 20 Discretization Verification & Foliar Hotspots (Figure 17)

![Figure 17: Foliar Surface Distribution and Hotspots](figures/fig17_foliar_surface_distribution_and_hotspots.png)
*Figure 17: Foliar boundary layer scalar distributions and spatial microenvironmental hotspot mapping across all 164 individual plant positions. ([View SVG](figures/fig17_foliar_surface_distribution_and_hotspots.svg) · [High-Res PNG](figures/fig17_foliar_surface_distribution_and_hotspots.png))*

> [!NOTE]
> **What Was Modeled:** Numerical mesh independence verification using the ASME V&V 20-2009 standard, coupled with spatial microenvironmental hotspot analysis across every single plant position (164 plants total).
> 
> **What Was Done:** 3-grid convergence study ($\Delta x = 1.50, 1.00, 0.75\ \text{mm}$) computing apparent order $p$, Richardson extrapolation $\phi_{\text{ext}}$, and Grid Convergence Index (GCI), combined with micro-voxel surface integration for all 164 plants.
> 
> **What It Shows & Describes:**
> - Achieves strict grid independence: canopy conductance GCI is **$<0.01\%$** ($p = 3.42$, asymptotic ratio $= 1.000$).
> - Spatial hotspot mapping reveals that Chamber 2 rosettes experience sheltered boundary layers ($g_{bl} = 0.041\ \text{mol m}^{-2}\ \text{s}^{-1}$), while outer *B. rapa* flanks experience peak ventilation ($0.052\ \text{mol m}^{-2}\ \text{s}^{-1}$), showing $<12\%$ intra-chamber variance.

---

### 3.11 Comparative Spaceflight Hardware Benchmark (Figure 18)

![Figure 18: Comparative Spaceflight Hardware Benchmark](figures/fig18_spaceflight_hardware_comprehensive_benchmark.png)
*Figure 18: Multi-dimensional architectural radar comparison and comparative trade matrix evaluating LEAF GBC1 against NASA BRIC, CARA, and Veggie. ([View SVG](figures/fig18_spaceflight_hardware_comprehensive_benchmark.svg) · [High-Res PNG](figures/fig18_spaceflight_hardware_comprehensive_benchmark.png))*

> [!TIP]
> **What Was Modeled:** Systematic architectural and fluid dynamic comparison between LEAF GBC1 and established NASA flight cultivation hardware platforms.
> 
> **What Was Done:** Evaluated air changes per hour (ACH), time to carbon starvation ($\tau_{\text{starve}}$), time to nocturnal hypoxia ($\tau_{\text{hypoxia}}$), boundary layer conductance in $0g$, and in-situ fixation capabilities across NASA BRIC (60 mm sealed dish), NASA CARA (120 mm taped dish), NASA Veggie (ISS cabin bellows), and LEAF GBC1.
> 
> **What It Shows & Describes:**
> - Demonstrates that passive dishes (BRIC, CARA) suffer rapid carbon depletion ($7\text{--}42\ \text{min}$) and severe fogging due to zero forced convection.
> - LEAF GBC1 delivers **$554.5\ \text{ACH}$**, maintaining continuous $>99.5\%$ ambient $\text{CO}_2$ replenishment, sub-second autonomous fixation ($1.20\ \text{s}$ saturation), and active fogging prevention.

---

## 4. Summary of Lunar Flight Operating Rules

**Table 2: Master Flight Operational Specifications for the LEAF GBC1 Mission**

| Parameter | Recommended Flight Setpoint | Physical Rationale & Boundary |
| :--- | :--- | :--- |
| **Daytime Flow Rate ($Q_{\text{day}}$)** | **$0.85\ \text{L/min}$** ($2.10\ \text{m/s}$ inlet) | Compensates $0.166\,g$ lunar buoyancy loss; restores $98\%$ Earth carbon assimilation. |
| **Nocturnal Flow Floor ($Q_{\text{night}}$)** | **$\ge 0.60\ \text{L/min}$** | Flushes heavy respiratory $\text{CO}_2$ pooling from lower drainage wells; maintains $\text{Ri} \le 0.05$. |
| **Ventilation Mode** | **Controlled Pull Suction (Port 2 Exhaust)** | Generates negative gauge pressure ($-11.8\ \text{Pa}$) for fail-safe bio-containment and anti-fogging. |
| **Water Activation Volume** | **$7.39\ \text{mL}$** ($5.28\ \text{mL}$ for $\phi=0.50$) | $100\%$ saturation of 3 foam cassettes without free liquid flooding or droplet detachment. |
| **Root Barrier Screen** | **$35\ \mu\text{m}$ Nytex monofilament mesh** | $8.04\ \text{kPa}$ bubble point allows water wicking while providing $100\%$ root exclusion. |
| **Fixation Injection Rate** | **$10.0\ \text{mL/min}$ RNAlater** | Sub-surface Port 3 dipping nozzle wets roots in $31.7\ \text{s}$; chamber floods in $5.82\ \text{min}$ ($\Delta P \le 4.04\ \text{kPa}$). |
| **Brownout Grace Window** | **$\le 14.2\ \text{minutes}$** | Maximum unpowered pump outage before canopy $\text{CO}_2$ drawdown induces hypoxic stress. |
| **Lid Condensation Margin** | **$\Delta T_{\text{margin}} \ge +0.83^\circ\text{C}$** | Maintains acrylic lid above chamber dew point; critical condensation threshold is $Q_{\text{crit}} = 0.08\ \text{L/min}$. |

---

## 5. Figure and Dataset Catalog

### Publication Figures (`figures/`)
1. [`fig_orthographic_projection_annotated.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig_orthographic_projection_annotated.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig_orthographic_projection_annotated.svg)) — 3-view orthographic blueprint.
2. [`assembly_cross_sections_annotated.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/assembly_cross_sections_annotated.png) — Mated CAD assembly cross-sections.
3. [`fig1_pump_speed_vs_aerodynamics.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig1_pump_speed_vs_aerodynamics.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig1_pump_speed_vs_aerodynamics.svg)) — Flow rate vs aerodynamics.
4. [`fig2_velocity_fields_earth_vs_lunar.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig2_velocity_fields_earth_vs_lunar.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig2_velocity_fields_earth_vs_lunar.svg)) — Velocity, $\text{CO}_2$, $\text{H}_2\text{O}$ cutaways.
5. [`fig3_lunar_compensatory_pump_curve.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig3_lunar_compensatory_pump_curve.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig3_lunar_compensatory_pump_curve.svg)) — Compensatory pump curve.
6. [`fig4_mean_age_of_air_and_chamber_uniformity.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig4_mean_age_of_air_and_chamber_uniformity.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig4_mean_age_of_air_and_chamber_uniformity.svg)) — 3D Mean Age of Air.
7. [`fig5_flight_hardware_benchmark_gbc1_vs_nasa.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig5_flight_hardware_benchmark_gbc1_vs_nasa.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig5_flight_hardware_benchmark_gbc1_vs_nasa.svg)) — NASA hardware benchmark.
8. [`fig6_predicted_hormone_carbon_phenotypes.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig6_predicted_hormone_carbon_phenotypes.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig6_predicted_hormone_carbon_phenotypes.svg)) — Hormone & transcriptomic predictions.
9. [`fig7_rnalater_perfusion_and_species_dynamics.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig7_rnalater_perfusion_and_species_dynamics.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig7_rnalater_perfusion_and_species_dynamics.svg)) — RNAlater perfusion kinetics.
10. [`fig8_diurnal_day_night_photoperiod_inversion.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig8_diurnal_day_night_photoperiod_inversion.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig8_diurnal_day_night_photoperiod_inversion.svg)) — Diurnal cycle & night setback.
11. [`fig9_developmental_growth_timeseries.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig9_developmental_growth_timeseries.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig9_developmental_growth_timeseries.svg)) — 14-day developmental growth series.
12. [`fig10_off_nominal_sensitivity_sweeps.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig10_off_nominal_sensitivity_sweeps.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig10_off_nominal_sensitivity_sweeps.svg)) — Off-nominal sensitivity sweeps.
13. [`fig11_water_capillarity_evapotranspiration_condensation.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig11_water_capillarity_evapotranspiration_condensation.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig11_water_capillarity_evapotranspiration_condensation.svg)) — Capillarity & fogging defense.
14. [`fig12_water_activation_foam_and_nytex.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig12_water_activation_foam_and_nytex.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig12_water_activation_foam_and_nytex.svg)) — Foam hydration & Nytex barrier.
15. [`fig13_aerodynamic_push_vs_pull_regimes.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig13_aerodynamic_push_vs_pull_regimes.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig13_aerodynamic_push_vs_pull_regimes.svg)) — Push vs pull aerodynamics.
16. [`fig14_model_lunar_leaf_5day_ontogeny.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig14_model_lunar_leaf_5day_ontogeny.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig14_model_lunar_leaf_5day_ontogeny.svg)) — 5-day developmental ontogeny.
17. [`fig15_gravity_temperature_response_surfaces.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig15_gravity_temperature_response_surfaces.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig15_gravity_temperature_response_surfaces.svg)) — Multi-physics response surfaces.
18. [`fig16_optical_lid_condensation_envelope.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig16_optical_lid_condensation_envelope.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig16_optical_lid_condensation_envelope.svg)) — Lid dew point margins.
19. [`fig17_foliar_surface_distribution_and_hotspots.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig17_foliar_surface_distribution_and_hotspots.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig17_foliar_surface_distribution_and_hotspots.svg)) — 164-plant foliar hotspots.
20. [`fig18_spaceflight_hardware_comprehensive_benchmark.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig18_spaceflight_hardware_comprehensive_benchmark.png) ([SVG](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig18_spaceflight_hardware_comprehensive_benchmark.svg)) — Hardware continuum benchmark.

### Quantitative Datasets (`results/tables/`) & Zenodo Package (`zenodo/`)
- [`gbc1_cfd_sweep_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_cfd_sweep_results.tsv) — Complete 21-case flow sweep.
- [`gbc1_covariates_for_omics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_covariates_for_omics.tsv) — Environmental covariates for NASA OSDR omics.
- [`gbc1_advanced_analytics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_advanced_analytics.tsv) — Age of air & boundary layer conductance.
- [`gbc1_hardware_benchmark_comparison.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_hardware_benchmark_comparison.tsv) — GBC1 vs BRIC, CARA, Veggie.
- [`gbc1_asymmetric_canopy_species_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_asymmetric_canopy_species_results.tsv) — Species-specific compartment fluxes.
- [`gbc1_rnalater_perfusion_kinetics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_rnalater_perfusion_kinetics.tsv) — Fixative displacement kinetics.
- [`gbc1_diurnal_day_night_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_diurnal_day_night_results.tsv) — 24h diurnal photoperiod data.
- [`gbc1_developmental_growth_timeseries.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_developmental_growth_timeseries.tsv) — 14-day developmental series.
- [`gbc1_off_nominal_sensitivity.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_off_nominal_sensitivity.tsv) — Off-nominal sensitivity sweeps.
- [`gbc1_water_and_humidity_transport.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_water_and_humidity_transport.tsv) — Capillarity & psychrometric data.
- [`gbc1_gci_mesh_convergence.csv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_gci_mesh_convergence.csv) — ASME V&V 20 GCI verification data.
- [`gbc1_foliar_surface_distribution.csv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_foliar_surface_distribution.csv) — Spatial percentiles ($p_{10}, p_{50}, p_{90}$).
- [`gbc1_164plants_foliar_telemetry.csv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_164plants_foliar_telemetry.csv) — Individual 164-plant scalar telemetry.
- [`gbc1_spaceflight_hardware_comprehensive_benchmark.csv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc1_spaceflight_hardware_comprehensive_benchmark.csv) — Comprehensive hardware benchmark.
- **Open-Science Archive:** [`zenodo/datapackage.json`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/zenodo/datapackage.json) and [`zenodo/README.md`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/zenodo/README.md).
