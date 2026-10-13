# Computational Fluid Dynamics, Multi-Species Canopy Gas Exchange, and Harvest Fixation Architecture of the LEAF GBC2 Lunar Plant Growth Chamber

**Authors:** Richard Barker et al.  
**Target Hardware:** LEAF Growth Chamber 2 (`M09024F-MDL-5` Pot Body & `M09025F-MDL-2` Lid Assembly)  
**Interactive Digital Twin & Research Portal:** [https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/](https://dr-richard-barker.github.io/LEAF_GBC2_pot_model/)  
**Repository:** [https://github.com/dr-richard-barker/LEAF_GBC2_pot_model](https://github.com/dr-richard-barker/LEAF_GBC2_pot_model)

---

## Abstract

Plant growth habitats deployed on lunar surface landers and orbital platforms operate in fractional gravity ($0.166\,g$ on the Moon) or microgravity ($\mu g$), where buoyancy-driven natural thermal and solutal convection is suppressed by $83.4\%$ to $100\%$. Without compensatory forced ventilation, stagnant boundary layers envelop foliar surfaces, throttling carbon dioxide ($\text{CO}_2$) assimilation, trapping transpired water vapor ($\text{H}_2\text{O}$), elevating leaf temperature, and inducing nocturnal respiratory hypercapnia. Here, we present a three-dimensional Computational Fluid Dynamics (CFD), porous-media multiphasic transport, and biochemical gas-exchange framework built directly from the flight-like Computer-Aided Design (CAD) geometry of the **Lunar Effects on Agricultural Flora (LEAF) Growth Chamber 2 (GBC2)** assembly (`M09024F-MDL-5.STEP` and `M09025F-MDL-2.STEP`). We model the exact three-compartment interior featuring a **3-part partitioned porous rooting substrate** ($Y = 7.62\text{--}14.62\ \text{mm}$, total void volume $5.28\ \text{mL}$ at porosity $\epsilon = 0.50$) and asymmetric multi-species planting densities: **32 *Brassica rapa* seeds** in each lateral compartment (Chambers 1 and 3; tall erect canopy reaching $Y = 26.0\ \text{mm}$, total foliar area $64.0\ \text{cm}^2/\text{chamber}$, $\text{LAI} = 9.44$) flanking **100 *Arabidopsis thaliana* seeds** in the central compartment beneath the exhaust port (Chamber 2; miniature prostrate rosette lawn reaching $Y = 19.5\ \text{mm}$, total foliar area $15.0\ \text{cm}^2$, $\text{LAI} = 3.09$).

Across 21 parametric flow sweeps ($Q = 0.05\text{--}4.0\ \text{L/min}$), 3D Mean Age-of-Air ($\theta$) partial differential equation (PDE) inversions, 14-day developmental growth time-series, and off-nominal spacecraft fault evaluations, we establish four core engineering findings:
1. **Aerodynamic "Dip-and-Converge" Channelling & Uniformity:** The stature differential between the tall *B. rapa* canopies ($\Delta Y = 11.4\ \text{mm}$) and the low-lying *A. thaliana* rosette lawn ($\Delta Y = 4.9\ \text{mm}$) creates a $14.2\ \text{mm}$ open central plenum under Port 2 ($Z = 0\ \text{mm}$). Forced air entering Port 1 ($Z = -21.59\ \text{mm}$) skims across Chamber 1, dips into the low-resistance central Chamber 2 plenum, and converges with return flow from Chamber 3, achieving $<1.0\%$ inter-chamber Mean Age-of-Air disparity ($\theta = 6.05\text{--}6.10\ \text{s}$ at $Q = 0.50\ \text{L/min}$; Sandberg conservation ratio $\theta_{\text{out}}/\tau = 0.938$).
2. **Lunar Buoyancy Suppression & Daytime Compensatory Flow:** At the nominal Earth pump rate ($Q = 0.50\ \text{L/min}$), Lunar gravity ($1.62\ \text{m/s}^2$) reduces leaf boundary-layer conductance ($g_{bl}$) by $36.3\%$ (from $0.0604$ to $0.0385\ \text{mol m}^{-2}\ \text{s}^{-1}$ in *B. rapa* and $0.0759$ to $0.0484\ \text{mol m}^{-2}\ \text{s}^{-1}$ in *A. thaliana*), depressing Farquhar–von Caemmerer–Berry (FvCB) net photosynthesis ($A_{\text{net}}$) to $81.5\%$ ($10.63\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$) and $83.8\%$ ($8.43\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$) of Earth $1g$ baselines while elevating leaf temperatures by $+2.42\ ^\circ\text{C}$. Increasing daytime flow to **$Q = 0.85\ \text{L/min}$** reduces canopy $\text{CO}_2$ drawdown and humidity accumulation by $37\%$, restoring $98\%$ of terrestrial microclimate scalar uniformity.
3. **Diurnal Photoperiod Inversion & Nocturnal Pump Setback Floor:** During the 8-hour dark skotoperiod, radiative leaf cooling ($\Delta T = -0.3\ \text{K}$) and dense mitochondrial $\text{CO}_2$ emission ($\beta_{\text{CO}_2} = -0.517$) invert buoyancy forces ($\text{Ri} < 0$), causing respiratory $\text{CO}_2$ to pool inside the lower drainage wells ($Y < 7.62\ \text{mm}$) and porous substrate. Under Lunar $0.166\,g$, low-power sleep setbacks ($Q \le 0.10\ \text{L/min}$) trigger severe canopy hypercapnia; maintaining a **nocturnal flow floor of $Q_{\text{night}} \ge 0.60\ \text{L/min}$** (or a $50\%$ setback to $0.425\ \text{L/min}$ from a $0.85\ \text{L/min}$ day baseline) prevents respiratory pooling while conserving transpired moisture ($+0.88\ \text{mmol m}^{-3}$).
4. **Bottom-Up RNAlater Perfusion & Fault Resilience:** Injecting saturated ammonium sulfate fixative (RNAlater; $\rho = 1250\ \text{kg/m}^3$, $\mu = 1.85\ \text{mPa}\cdot\text{s}$) at $Q_{\text{inj}} = 10\ \text{mL/min}$ through the Port 3 dipping nozzle ($Y = 9.5\ \text{mm}$ inside Substrate Block 3) saturates the 3-part root matrix in **$31.7\ \text{s}$** and submerges the entire $64.6\ \text{cm}^3$ chamber in **$6.46\ \text{min}$** at a peak pressure drop of **$0.73\ \text{kPa}$** ($4.04\ \text{kPa}$ after 14 days of root-driven porosity compaction to $\epsilon = 0.34$), remaining safely below the $5.0\ \text{kPa}$ peristaltic seal safety limit.

---

## 1. Introduction

Spaceflight plant biology experiments conducted aboard the International Space Station (ISS)—including NASA's Biological Research in Canisters (BRIC), Petri-dish fixation hardware (CARA), and the Vegetable Production System (Veggie)—have demonstrated that reduced gravity alters plant transcriptomes, cell wall remodeling, phytohormone signaling (auxin, abscisic acid, and ethylene), and photosynthetic carbon allocation. A major physical driver of these spaceflight phenotypes is not direct gravitropism alone, but the secondary **fluid-dynamic microclimate artifact** caused by the loss of buoyancy-driven natural convection:

$$\text{Gr}_L = \frac{g \, \beta_{\text{eff}} \, \Delta \Phi \, L_{\text{leaf}}^3}{\nu^2}$$

On Earth ($g = 9.81\ \text{m/s}^2$), warm, humidified air parcels adjacent to illuminated leaves rise rapidly via thermal ($\beta_T = 3.39 \times 10^{-3}\ \text{K}^{-1}$) and solutal water-vapor ($\beta_{\text{H}_2\text{O}} = +0.379$) buoyancy, continuously stripping away the diffusive boundary layer even when mechanical fan velocities are low. On the lunar surface ($g = 1.62\ \text{m/s}^2$, $0.166\,g$), buoyant acceleration drops by $83.4\%$, thickening the laminar sublayer $\delta_{bl}$ and shifting the mixed-convection Richardson number ($\text{Ri} = \text{Gr}/\text{Re}^2$) toward diffusion-limited regimes unless forced ventilation is specifically engineered to compensate.

The **LEAF Growth Chamber 2 (GBC2)** hardware addresses these challenges through an actively ventilated, hermetically sealed polycarbonate growth pot (`M09024F-MDL-5`) mated with a three-port fluidic distribution lid (`M09025F-MDL-2`). Unlike single-species chambers, the GBC2 mission architecture co-cultivates two model Brassicaceae species with starkly contrasting growth habits inside three baffled compartments sharing a common overhead air plenum:
- **Chambers 1 and 3 (Lateral Compartments):** Each seeded with **32 seeds of *Brassica rapa*** (fast-growing, broad-leaved, upright stems reaching $Y = 26.0\ \text{mm}$).
- **Chamber 2 (Central Compartment):** Seeded with **100 seeds of *Arabidopsis thaliana*** (miniature, high-density prostrate rosette lawn remaining below $Y = 19.5\ \text{mm}$).

Furthermore, end-of-mission transcriptomic preservation requires in-situ liquid fixation via **RNAlater** injected into the root zone to displace headspace gas through the central exhaust port without causing premature foliar splashing, gas pocket trapping, or mechanical stem shearing.

In this work, we construct a high-resolution 3D voxelized CFD and multiphasic transport digital twin directly from the STEP CAD assemblies of the LEAF GBC2 chamber. We quantify the coupled aerodynamics, residence-time distributions, species-specific carbon/water fluxes, diurnal photoperiod inversions, 14-day developmental ontogeny, and off-nominal spacecraft fault envelopes required to certify the GBC2 flight software and pump operating schedule.

---

## 2. Materials and Methods

### 2.1 CAD Geometry Voxelization & Coordinate Reference Frame

Exact 3D solid models of the GBC2 pot body ([`geometry/M09024F-MDL-5.STEP`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/geometry/M09024F-MDL-5.STEP)) and lid assembly ([`geometry/M09025F-MDL-2.STEP`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/geometry/M09025F-MDL-2.STEP)) were tessellated and mated in a right-handed Cartesian world frame ($X$ = transverse width, $Y$ = vertical elevation against gravity, $Z$ = longitudinal length across the 3 chambers) via [`sim/build_geometry.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/build_geometry.py):

- **Pot Body (`M09024F-MDL-5`):** External envelope $47.63\ \text{mm}$ ($X$) $\times 33.67\ \text{mm}$ ($Y$) $\times 67.95\ \text{mm}$ ($Z$). Total internal fluid cavity volume $V_{\text{total}} = 64.60\ \text{cm}^3$ ($54.05\ \text{cm}^3$ free headspace plus $10.55\ \text{cm}^3$ bulk porous substrate containing $5.28\ \text{cm}^3$ pore void volume).
- **Lid Assembly (`M09025F-MDL-2`):** External envelope $45.09\ \text{mm}$ ($X$) $\times 5.71\ \text{mm}$ ($Y$) $\times 65.41\ \text{mm}$ ($Z$), rotated $90^\circ$ about the $Y$-axis so its three fluidic ports align longitudinally along $Z$ at transverse offset $X = -15.49\ \text{mm}$, seated on the upper pot rim at $Y = 33.67\ \text{mm}$ (top of lid at $Y = 39.38\ \text{mm}$).
- **Lower Drainage & Aeration Troughs ($Y = 0.00\text{--}7.62\ \text{mm}$):** Molded bottom channels ($V_{\text{drain}} = 8.74\ \text{cm}^3$) beneath the structural chamber ledges.
- **Internal Baffles & 3-Part Partitioned Porous Substrate ($Y = 7.62\text{--}14.62\ \text{mm}$):** Two vertical structural dividers at $Z \in [-9.5, -7.5]\ \text{mm}$ and $Z \in [+7.5, +9.5]\ \text{mm}$ ($Y = 7.62\text{--}24.0\ \text{mm}$) partition the interior into three chambers. Resting on the $Y = 7.62\ \text{mm}$ internal ledges are **three discrete porous media blocks** ($\Delta Y = 7.00\ \text{mm}$ thickness, initial porosity $\epsilon_0 = 0.50$, Darcy permeability $K_0 = 1.0 \times 10^{-10}\ \text{m}^2$):
  1. **Substrate Block 1 (Chamber 1):** $X \in [-11.0, +21.3]\ \text{mm}$, $Z \in [-30.5, -9.5]\ \text{mm}$ ($V_{\text{bulk}} = 4.02\ \text{cm}^3$, $V_{\text{void}} = 2.01\ \text{mL}$).
  2. **Substrate Block 2 (Chamber 2):** $X \in [-11.0, +21.3]\ \text{mm}$, $Z \in [-7.5, +7.5]\ \text{mm}$ ($V_{\text{bulk}} = 2.51\ \text{cm}^3$, $V_{\text{void}} = 1.26\ \text{mL}$).
  3. **Substrate Block 3 (Chamber 3):** $X \in [-11.0, +21.3]\ \text{mm}$, $Z \in [+9.5, +30.5]\ \text{mm}$ ($V_{\text{bulk}} = 4.02\ \text{cm}^3$, $V_{\text{void}} = 2.01\ \text{mL}$).
- **Fluidic Ports ($\varnothing 3.18\ \text{mm}$ / $1/8\text{"}$ ID along $X = -15.49\ \text{mm}$):**
  - **Port 1 ($Z = -21.59\ \text{mm}$):** Forced air inlet discharging into the upper perimeter distribution channel ($Y = 28\text{--}32\ \text{mm}$).
  - **Port 2 ($Z = 0.00\ \text{mm}$):** Central air exhaust and harvest gas vent positioned directly over Chamber 2.
  - **Port 3 ($Z = +21.59\ \text{mm}$):** Water and RNAlater injection port fitted with an internal vertical dipping nozzle ($\varnothing 2.0\ \text{mm}$, length $24.2\ \text{mm}$) terminating at $Y = 9.5\ \text{mm}$ inside the core of Substrate Block 3.

The domain was discretized onto a Cartesian finite-volume grid ([`geometry/gbc2_grid.npz`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/geometry/gbc2_grid.npz)) with isotropic voxel pitch $\Delta x = 1.0\ \text{mm}$ ($N_x \times N_y \times N_z = 48 \times 34 \times 68 = 110{,}976$ voxels, of which $64{,}600$ are active fluid/porous/canopy cells).

---

### 2.2 Asymmetric Multi-Species Canopy Morphometrics

Table 1 summarizes the physical and physiological parameterization of the three growth chambers at Day 14 harvest maturity.

**Table 1: Multi-Species Compartment Morphometrics and Aerodynamic Parameters in the LEAF GBC2 Assembly**

| Parameter | Symbol | Chamber 1 (*Brassica rapa*) | Chamber 2 (*Arabidopsis thaliana*) | Chamber 3 (*Brassica rapa*) | Units |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Seed Count** | $N_{\text{seeds}}$ | **32** | **100** | **32** | seeds |
| **Longitudinal Span** | $Z$ | $[-30.5, -9.5]$ | $[-7.5, +7.5]$ | $[+9.5, +30.5]$ | $\text{mm}$ |
| **Floor Footprint Area** | $A_{\text{floor}}$ | $6.78$ ($678.3$) | $4.85$ ($484.5$) | $6.78$ ($678.3$) | $\text{cm}^2$ ($\text{mm}^2$) |
| **Planting Density** | $\rho_{\text{seed}}$ | $47.2$ ($472$) | $206.2$ ($2{,}064$) | $47.2$ ($472$) | $\text{seeds/dm}^2$ ($\text{m}^{-2}$) |
| **Mean Inter-Seed Spacing** | $s_{\text{seed}}$ | $4.60$ | $2.20$ | $4.60$ | $\text{mm}$ |
| **Canopy Elevation Span** | $Y_{\text{canopy}}$ | $[14.62, 26.00]$ | $[14.62, 19.50]$ | $[14.62, 26.00]$ | $\text{mm}$ |
| **Canopy Stature Height** | $h_{\text{canopy}}$ | **$11.38$** | **$4.88$** | **$11.38$** | $\text{mm}$ |
| **Free Headspace to Lid** | $\Delta Y_{\text{plenum}}$ | $7.67$ | **$14.17$** | $7.67$ | $\text{mm}$ |
| **Individual Plant Leaf Area** | $A_{\text{ind}}$ | $2.00$ | $0.15$ | $2.00$ | $\text{cm}^2/\text{plant}$ |
| **Total Chamber Leaf Area** | $A_{\text{tot}}$ | **$64.00$** | **$15.00$** | **$64.00$** | $\text{cm}^2$ |
| **Share of Total Pot Leaf Area** | $f_{\text{area}}$ | $44.76\%$ | $10.49\%$ | $44.76\%$ | $\%$ (of $143.0\ \text{cm}^2$) |
| **Leaf Area Index** | $\text{LAI}$ | $9.44$ | $3.09$ | $9.44$ | $\text{m}^2\ \text{leaf}/\text{m}^2\ \text{floor}$ |
| **Leaf Area Density** | $a_{\text{leaf}}$ | $828.0$ | $630.0$ | $828.0$ | $\text{m}^2/\text{m}^3$ |
| **Characteristic Leaf Length** | $L_{\text{leaf}}$ | $20.0$ | $8.0$ | $20.0$ | $\text{mm}$ |
| **Canopy Drag Coefficient** | $C_d$ | $0.20$ | $0.15$ | $0.20$ | dimensionless |
| **Linearized Brinkman Drag** | $D_{\text{canopy}}$ | $8.28$ | $4.72$ | $8.28$ | $\text{s}^{-1}$ |

---

### 2.3 Governing Fluid, Transport, and Biochemical Equations

#### 2.3.1 Incompressible Navier–Stokes with Boussinesq Thermal & Solutal Buoyancy
Airflow within the headspace, porous substrate, and plant canopies ([`sim/gbc2_cfd_engine.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/gbc2_cfd_engine.py)) is governed by continuity and the Brinkman-penalized momentum equations:

$$\nabla \cdot \mathbf{u} = 0$$

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho_0}\nabla p + \nu \nabla^2 \mathbf{u} - \mathcal{D}(\mathbf{x})\mathbf{u} + \mathbf{g} \left[ \beta_T (T - T_0) + \beta_{\text{H}_2\text{O}}(C_{\text{H}_2\text{O}} - C_{\text{H}_2\text{O},0}) + \beta_{\text{CO}_2}(C_{\text{CO}_2} - C_{\text{CO}_2,0}) \right]$$

where $\rho_0 = 1.204\ \text{kg/m}^3$, $\nu = 1.51 \times 10^{-5}\ \text{m}^2/\text{s}$ at $22\ ^\circ\text{C}$ and $1\ \text{atm}$, and the spatially varying momentum sink $\mathcal{D}(\mathbf{x})$ equals $\frac{\nu}{\epsilon K}$ inside the three porous substrate blocks ($Y \in [7.62, 14.62]\ \text{mm}$) and $C_d \, a_{\text{leaf}}(\mathbf{x}) \, U_{\text{eff}}$ inside the differentiated plant canopies. Gravity $\mathbf{g} = [0, -g, 0]^T$ is evaluated across $g = 9.81\ \text{m/s}^2$ (Earth $1.0\,g$), $g = 1.62\ \text{m/s}^2$ (Lunar $0.166\,g$), and $g = 0.0\ \text{m/s}^2$ ($\mu g$).

#### 2.3.2 Multi-Species Advection–Diffusion & Mean Age of Air ($\theta$)
Dissolved/gaseous scalar deviations for $\text{CO}_2$ ($D_{\text{CO}_2} = 1.60 \times 10^{-5}\ \text{m}^2/\text{s}$) and water vapor ($D_{\text{H}_2\text{O}} = 2.42 \times 10^{-5}\ \text{m}^2/\text{s}$) satisfy:

$$\frac{\partial C_k}{\partial t} + \mathbf{u} \cdot \nabla C_k = D_k \nabla^2 C_k + S_k(\mathbf{x})$$

where volumetric sink/source rates $S_{\text{CO}_2}$ and $S_{\text{H}_2\text{O}}$ scale with total chamber foliar area ($-40.0\ \mu\text{mol m}^{-3}\ \text{s}^{-1}$ and $+80.0\ \text{mmol m}^{-3}\ \text{s}^{-1}$ in Chambers 1 & 3 vs. $-15.0\ \mu\text{mol m}^{-3}\ \text{s}^{-1}$ and $+30.0\ \text{mmol m}^{-3}\ \text{s}^{-1}$ in Chamber 2 during daytime illumination).

To rigorously verify chamber ventilation efficiency and detect dead zones ([`sim/advanced_analytics.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/advanced_analytics.py)), we solve the steady-state **Mean Age of Air** transport equation:

$$\mathbf{u} \cdot \nabla \theta - D_\theta \nabla^2 \theta = 1, \qquad \theta\big|_{\text{Port 1}} = 0, \qquad \frac{\partial \theta}{\partial n}\Bigg|_{\text{walls}} = 0$$

and validate global mass conservation against the **Sandberg Residence-Time Identity**:

$$\langle \theta \rangle_{\text{outlet, flux}} = \frac{\int_{\text{Port 2}} (\mathbf{u} \cdot \mathbf{n}) \, \theta \, dA}{\int_{\text{Port 2}} (\mathbf{u} \cdot \mathbf{n}) \, dA} \equiv \tau_{\text{nominal}} = \frac{V_{\text{fluid}}}{Q}$$

#### 2.3.3 Mixed-Convection Boundary-Layer Conductance & FvCB Photosynthesis
Within each compartment $c \in \{1, 2, 3\}$, the effective laminar leaf boundary-layer thickness $\delta_{\text{eff}, c}$ blends forced Blasius flat-plate convection ($\delta_{\text{forced}, c} = 4.92 \sqrt{\nu L_{\text{leaf}, c} / |\mathbf{u}_c|}$) and natural buoyant Grashof convection ($\delta_{\text{nat}, c} = 3.5 L_{\text{leaf}, c} \text{Gr}_c^{-1/4}$) via Churchill's asymptotic cube-root law:

$$\delta_{\text{eff}, c} = \left( \delta_{\text{forced}, c}^{-3} + \delta_{\text{nat}, c}^{-3} \right)^{-1/3}, \qquad g_{bl, c} = \frac{D_{\text{CO}_2}}{\delta_{\text{eff}, c} \, V_{\text{mol}}}$$

Net carbon assimilation $A_{\text{net}}$ ($\mu\text{mol CO}_2\ \text{m}^{-2}\ \text{s}^{-1}$) is obtained by analytically solving the coupled quadratic **Farquhar–von Caemmerer–Berry (FvCB)** biochemical demand and boundary-plus-stomatal supply equations ($g_{\text{tot}}^{-1} = g_{bl}^{-1} + g_s^{-1}$):

$$A_{\text{net}} = \min(A_c, A_j, A_p) - R_d, \qquad A = g_{\text{tot}}(C_s - C_i)$$

where $A_c$ is the Rubisco carboxylation-limited rate, $A_j$ is the RuBP electron-transport-limited rate, $A_p$ is the triose-phosphate utilization (TPU) ceiling at super-elevated cabin $\text{CO}_2$, and $C_s = C_{\text{amb}} - \Delta[\text{CO}_2]_{\text{canopy}}$.

#### 2.3.4 Two-Stage Harvest RNAlater Perfusion & Kozeny–Carman Root Compaction
During end-of-mission harvest ([`sim/simulate_rnalater_perfusion.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/simulate_rnalater_perfusion.py)), aqueous RNAlater ($\rho_{\text{RNA}} = 1250\ \text{kg/m}^3$, $\mu_{\text{RNA}} = 1.85 \times 10^{-3}\ \text{Pa}\cdot\text{s}$, surface tension $\sigma = 0.078\ \text{N/m}$) is pumped into Port 3 and discharged at $Y = 9.5\ \text{mm}$ inside Substrate Block 3. Stage 1 models radial and inter-chamber Darcy imbibition across the three substrate blocks ($V_{\text{void}} = 5.28\ \text{mL}$), while Stage 2 tracks upward Volume-of-Fluid (VOF) free-surface displacement through the $54.05\ \text{mL}$ headspace toward the Port 2 exhaust vent. Root volume expansion over the 14-day mission compresses substrate porosity from $\epsilon_0 = 0.50$ to $\epsilon(t)$, reducing permeability according to the **Kozeny–Carman relation**:

$$K(\epsilon) = K_0 \frac{\epsilon^3 (1 - \epsilon_0)^2}{\epsilon_0^3 (1 - \epsilon)^2}$$

---

## 3. Results

### 3.1 Hardware Fluidics Architecture & ABAI-Verified Orthographic Blueprints

Figure 0 ([`figures/fig_orthographic_projection_annotated.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig_orthographic_projection_annotated.png)) and the exact 3D world-coordinate CAD cross-sections ([`figures/assembly_cross_sections_annotated.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/assembly_cross_sections_annotated.png) and [`figures/assembly_cross_sections_exact.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/assembly_cross_sections_exact.png)) resolve the mated mechanical interfaces of the GBC2 pot and lid. Key architectural features verified by three-plane slicing include:
- **Anti-Splash Internal Dipping Nozzle (Port 3, $Z = +21.59\ \text{mm}$):** Terminates at $Y = 9.5\ \text{mm}$ inside Substrate Block 3, ensuring that both nutrient water and harvest-phase RNAlater wet the porous matrix from the bottom up rather than dripping onto leaves.
- **3-Part Partitioned Substrate Bed ($Y = 7.62\text{--}14.62\ \text{mm}$):** Seated cleanly on the molded horizontal ledges above the $7.62\ \text{mm}$ lower drainage troughs, isolated into three compartments by the solid vertical baffles at $Z = \pm 8.5\ \text{mm}$.
- **Central Plenum Dip Geometry:** Because the 100 *A. thaliana* rosettes in Chamber 2 remain below $Y = 19.5\ \text{mm}$ while the 32 *B. rapa* plants in Chambers 1 and 3 reach $Y = 26.0\ \text{mm}$, an unobstructed $14.17\ \text{mm}$ vertical plenum sits directly beneath the Port 2 exhaust stack ($Z = 0.0\ \text{mm}$).

---

### 3.2 Parametric Pump Speed Sweeps & Lunar Buoyancy Compensation (Figures 1–3)

Simulations across 16 flow rates ($Q = 0.05\text{--}4.0\ \text{L/min}$) under Earth ($1.0\,g$), Lunar ($0.166\,g$), and microgravity ($\mu g$) conditions ([`results/tables/gbc2_cfd_sweep_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_cfd_sweep_results.tsv), [`results/tables/gbc2_asymmetric_canopy_species_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_asymmetric_canopy_species_results.tsv)) reveal a pronounced gravity-dependent divergence at low-to-moderate pump speeds ([`figures/fig1_pump_speed_vs_aerodynamics.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig1_pump_speed_vs_aerodynamics.png), [`figures/fig2_velocity_fields_earth_vs_lunar.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig2_velocity_fields_earth_vs_lunar.png)):

- **Earth $1.0\,g$ Baseline ($Q = 0.50\ \text{L/min}$):** Buoyant thermal and water-vapor plumes ($\text{Gr} \approx 1.38 \times 10^3$) augment forced airflow, maintaining thin boundary layers ($\delta_{bl} = 10.81\ \text{mm}$ in *B. rapa*, $8.60\ \text{mm}$ in *A. thaliana*) and high boundary-layer conductance ($g_{bl} = 0.0604$ and $0.0759\ \text{mol m}^{-2}\ \text{s}^{-1}$). Net photosynthesis reaches **$13.04\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$** in Chambers 1 & 3 (*B. rapa*) and **$10.06\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$** in Chamber 2 (*A. thaliana*). Owing to its $4.27\times$ larger foliar area per chamber, the 64 *B. rapa* plants across Chambers 1 and 3 account for **$89.5\%$ of total chamber carbon fixation**, whereas the 100 *A. thaliana* seedlings account for **$10.5\%$**.
- **Lunar $0.166\,g$ Uncompensated ($Q = 0.50\ \text{L/min}$):** The $83.4\%$ reduction in gravity suppresses buoyant plume acceleration, expanding $\delta_{bl}$ to $16.97\ \text{mm}$ (*B. rapa*) and $13.50\ \text{mm}$ (*A. thaliana*), cutting $g_{bl}$ by $36.3\%$ to $0.0385$ and $0.0484\ \text{mol m}^{-2}\ \text{s}^{-1}$, and reducing net photosynthesis to **$10.63\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$ ($81.5\%$ of Earth)** in *B. rapa* and **$8.43\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$ ($83.8\%$ of Earth)** in *A. thaliana*.
- **Lunar $0.166\,g$ Compensated Setpoint ($Q = 0.85\ \text{L/min}$):** Increasing the air pump speed by $1.7\times$ to $Q = 0.85\ \text{L/min}$ ([`figures/fig3_lunar_compensatory_pump_curve.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig3_lunar_compensatory_pump_curve.png)) shifts the chamber into the forced-convection dominant regime ($\text{Ri} = 0.048 < 0.10$, mass Péclet number $\text{Pe} = 28.4$), slashing canopy $\text{CO}_2$ depletion from $-3.77$ to $-2.35\ \text{ppm}$ ($-37\%$) and transpired vapor accumulation from $+7.12$ to $+4.50\ \text{mmol m}^{-3}$ ($-37\%$) in *B. rapa*.

---

### 3.3 3D Mean Age of Air, Hardware Benchmarking & Omics Projections (Figures 4–6)

Solving the 3D Mean Age-of-Air PDE across the voxelized domain ([`figures/fig4_mean_age_of_air_and_chamber_uniformity.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig4_mean_age_of_air_and_chamber_uniformity.png), [`results/tables/gbc2_advanced_analytics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_advanced_analytics.tsv)) confirms that the upper perimeter manifold and central plenum dip deliver exceptional longitudinal symmetry:
- At $Q = 0.50\ \text{L/min}$, the theoretical flushing time is $\tau = V_{\text{fluid}}/Q = 7.75\ \text{s}$, and the numerical flux-weighted outlet age at Port 2 is $\theta_{\text{out}} = 7.27\ \text{s}$ (**Sandberg ratio $= 0.938$**).
- Mean canopy residence times across the three compartments are **$\theta_{\text{Ch1}} = 6.05\ \text{s}$**, **$\theta_{\text{Ch2}} = 6.10\ \text{s}$**, and **$\theta_{\text{Ch3}} = 6.06\ \text{s}$**, representing an inter-chamber ventilation disparity of **$< 1.0\%$**.
- Benchmarking GBC2 against passive ISS habitats ([`figures/fig5_flight_hardware_benchmark_gbc2_vs_nasa.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig5_flight_hardware_benchmark_gbc2_vs_nasa.png), [`results/tables/gbc2_hardware_benchmark_comparison.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_hardware_benchmark_comparison.tsv)) shows that passive Petri canisters (BRIC-19, CARA; $\text{Pe} < 0.1$, $g_{bl} \approx 0.004\ \text{mol m}^{-2}\ \text{s}^{-1}$) suffer severe hypoxia and ethylene/ABA stress, whereas actively pumped GBC2 ($Q = 0.50\text{--}0.85\ \text{L/min}$, $\text{Pe} = 16.7\text{--}28.4$) suppresses predicted transcriptomic stress z-scores for ABA, ethylene, and oxidative chaperones by $>65\%$ ([`figures/fig6_predicted_hormone_carbon_phenotypes.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig6_predicted_hormone_carbon_phenotypes.png)).

---

### 3.4 Harvest-Phase RNAlater Perfusion Kinetics (Figure 7)

Multiphasic Darcy–VOF perfusion simulations ([`figures/fig7_rnalater_perfusion_and_species_dynamics.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig7_rnalater_perfusion_and_species_dynamics.png), [`results/tables/gbc2_rnalater_perfusion_kinetics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_rnalater_perfusion_kinetics.tsv)) demonstrate that partitioning the rooting matrix into three discrete blocks on the $Y = 7.62\ \text{mm}$ ledges cuts root saturation time by more than half compared to a continuous bottom-filled pot:
- **Nominal Injection ($Q_{\text{inj}} = 10\ \text{mL/min}$):** The $5.28\ \text{mL}$ substrate pore volume saturates in **$t_{\text{root}} = 31.7\ \text{s}$ ($0.53\ \text{min}$)**, followed by a coherent, capillary-stabilized upward liquid front ($\text{Bo} = 0.106$ on the Moon vs. $0.640$ on Earth) that submerges the Chamber 2 *Arabidopsis* lawn at $t = 2.1\ \text{min}$, the Chamber 1 & 3 *B. rapa* canopies at $t = 4.8\ \text{min}$, and reaches Port 2 exhaust ($64.6\ \text{cm}^3$ total volume) in **$t_{\text{flood}} = 387.6\ \text{s}$ ($6.46\ \text{min}$)** with a total pressure drop of only **$\Delta P = 0.73\ \text{kPa}$**.
- **Rapid Protocol ($Q_{\text{inj}} = 20\ \text{mL/min}$):** Completes root saturation in **$15.8\ \text{s}$** and full chamber fixation in **$3.23\ \text{min}$** at **$\Delta P = 1.47\ \text{kPa}$**.

---

### 3.5 Diurnal Photoperiod Inversion & Nocturnal Pump Setback Rules (Figure 8)

Simulating the 16-hour day photoperiod ($\Delta T_{\text{leaf}} = +1.5\ \text{K}$, active photosynthesis and transpiration) versus the 8-hour night skotoperiod ($\Delta T_{\text{leaf}} = -0.3\ \text{K}$, $85\%$ stomatal closure, dark mitochondrial respiration $S_{\text{CO}_2} > 0$) reveals a fundamental **buoyancy sign inversion** ([`figures/fig8_diurnal_day_night_photoperiod_inversion.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig8_diurnal_day_night_photoperiod_inversion.png), [`results/tables/gbc2_diurnal_day_night_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_diurnal_day_night_results.tsv)):
- At night, cool leaf surfaces and dense respired $\text{CO}_2$ ($M_{\text{CO}_2} = 44.01\ \text{g/mol} > M_{\text{air}} = 28.97\ \text{g/mol}$, $\beta_{\text{CO}_2} = -0.517$) sink downward into the substrate blocks and lower drainage troughs ($Y < 7.62\ \text{mm}$).
- If the air pump is throttled down to a severe low-power sleep mode ($Q = 0.10\ \text{L/min}$) under Lunar $0.166\,g$, respiratory $\text{CO}_2$ pools in the canopy and root zone ($\Delta[\text{CO}_2] > +9.5\ \text{ppm}$ locally).
- **Flight Rule:** Implementing a **nocturnal flow floor of $Q_{\text{night}} \ge 0.60\ \text{L/min}$** (or a controlled $50\%$ nocturnal setback from $0.85$ to $0.425\text{--}0.60\ \text{L/min}$) keeps $|\text{Ri}| < 0.05$, completely flushing respired $\text{CO}_2$ while retaining $+0.88\ \text{mmol m}^{-3}$ of beneficial nocturnal humidity.

---

### 3.6 14-Day Developmental Growth Time-Series: Imbibition to Harvest Closure (Figure 9)

To guide pump scheduling across the entire mission timeline, we simulated the daily ontogenetic trajectory from Day 0 (seed imbibition) through Day 14 (canopy closure) ([`sim/simulate_developmental_growth_timeseries.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/simulate_developmental_growth_timeseries.py), [`figures/fig9_developmental_growth_timeseries.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig9_developmental_growth_timeseries.png), [`results/tables/gbc2_developmental_growth_timeseries.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_developmental_growth_timeseries.tsv)):
1. **Days 0–2 (Imbibition & Radicle Emergence):** Seeds exhibit pure dark respiration ($A_{\text{net}} < 0$, whole-chamber flux $-0.45$ to $-0.72\ \text{nmol CO}_2/\text{s}$). Operating the pump at a gentle **germination setback ($Q = 0.25\ \text{L/min}$)** flushes metabolic $\text{CO}_2$ while preventing premature desiccation of the exposed seed coats on the top surface of Substrate Blocks 1, 2, and 3 ($Y = 14.62\ \text{mm}$).
2. **Day 2.8–3.2 (Autotrophic Transition):** Following cotyledon unfolding and chlorophyll synthesis, daily integrated carbon balance crosses zero at **Day 2.8** in *B. rapa* and **Day 3.2** in *A. thaliana*.
3. **Days 7–14 (Canopy Blockage Channelling & Cumulative Carbon Gain):** As the 32 *B. rapa* plants elongate to $h = 11.4\ \text{mm}$ ($\text{LAI} = 9.44$), aerodynamic drag in Chambers 1 and 3 increases velocity channelling into the central Chamber 2 plenum from $|\mathbf{u}|_{\text{Ch2}}/|\mathbf{u}|_{\text{Ch1}} = 0.443$ (Day 0) to $0.794$ (Day 14). Over the 14-day mission, cumulative net carbon fixed per *B. rapa* chamber reaches **$383.7\ \text{mg C}$ on Earth $1g$** ($65.7\ \text{mg C}$ in Chamber 2 *A. thaliana*), drops to **$311.7\ \text{mg C}$ ($81.2\%$)** under uncompensated Lunar $0.166\,g$ ($Q = 0.50\ \text{L/min}$), and recovers to **$312.6\ \text{mg C}$ with $41\%$ lower boundary-layer $\text{CO}_2$/humidity drawdown** under the **Adaptive Lunar Pump Schedule** ($Q(t) = 0.25 \rightarrow 0.50 \rightarrow 0.65 \rightarrow 0.85\ \text{L/min}$).
4. **Root-Driven Porosity Compaction:** Root volume expansion reduces substrate porosity from $\epsilon_0 = 0.50$ at Day 0 to **$\epsilon(14) = 0.34$** in Substrate Blocks 1 & 3 (*B. rapa*) and **$\epsilon(14) = 0.41$** in Substrate Block 2 (*A. thaliana*). By Kozeny–Carman scaling, Darcy permeability in Chambers 1 & 3 decreases from $1.00 \times 10^{-10}\ \text{m}^2$ to $1.81 \times 10^{-11}\ \text{m}^2$, elevating harvest RNAlater injection pressure at $Q_{\text{inj}} = 10\ \text{mL/min}$ from **$0.73\ \text{kPa}$ (Day 0) to $4.04\ \text{kPa}$ (Day 14)**—remaining safely below the $5.0\ \text{kPa}$ hardware seal ceiling while shortening root-zone saturation time from $31.6\ \text{s}$ to **$22.6\ \text{s}$**.

---

### 3.7 Off-Nominal Mission Sensitivity & Fault-Tolerance Sweeps (Figure 10)

Finally, we evaluated three off-nominal spacecraft environmental and hardware fault modes ([`sim/simulate_off_nominal_sensitivity.py`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/sim/simulate_off_nominal_sensitivity.py), [`figures/fig10_off_nominal_sensitivity_sweeps.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig10_off_nominal_sensitivity_sweeps.png), [`results/tables/gbc2_off_nominal_sensitivity.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_off_nominal_sensitivity.tsv)):
1. **Elevated Spacecraft Cabin $\text{CO}_2$ ($420\text{--}5{,}000\ \text{ppm}$):** In crewed spacecraft and lunar lander cabins ($2{,}000\text{--}4{,}500\ \text{ppm}\ \text{CO}_2$), elevated ambient $\text{CO}_2$ saturates Rubisco carboxylation and reaches the triose-phosphate utilization (TPU) plateau ($A_{\text{net}} \approx 29.9\text{--}31.6\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$ in *B. rapa* and $23.4\text{--}24.4\ \mu\text{mol m}^{-2}\ \text{s}^{-1}$ in *A. thaliana* above $2{,}500\ \text{ppm}$), masking the photosynthetic carbon penalty of thinner lunar boundary-layer conductance. **However**, super-ambient $\text{CO}_2$ induces Ball–Berry/Medlyn stomatal closure ($g_s$ drops by $50\%$ from $0.165$ to $0.082\ \text{mol m}^{-2}\ \text{s}^{-1}$ at $4{,}000\ \text{ppm}$), curtailing latent evaporative cooling and driving Lunar leaf temperature elevation up to **$\Delta T_{\text{leaf}} = +3.05\ ^\circ\text{C}$** at $Q = 0.50\ \text{L/min}$ vs. **$+2.19\ ^\circ\text{C}$** at $Q = 0.85\ \text{L/min}$. Thus, compensatory high-flow ventilation ($0.85\ \text{L/min}$) remains essential at elevated cabin $\text{CO}_2$ for **foliar thermal regulation**.
2. **Severe Root Compaction ($\epsilon = 0.50 \rightarrow 0.20$):** Below $\epsilon = 0.32$, Millington–Quirk gas tortuosity ($\tau_{\text{gas}} \propto \epsilon^{4/3}$) drives lunar rhizosphere $\text{O}_2$ depletion above $14\%$ and pushes $10\ \text{mL/min}$ RNAlater injection pressure past the $5.0\ \text{kPa}$ seal limit at $\epsilon \le 0.32$. For over-compacted root masses ($\epsilon \in [0.25, 0.32]$), stepping down harvest RNAlater pump flow from $10$ to **$5\ \text{mL/min}$** halves $\Delta P$ back below $3.5\ \text{kPa}$.
3. **Pump Brownout & Power Outage Grace Window ($Q = 0.85 \rightarrow 0.02\ \text{L/min}$):** Under Lunar $0.166\,g$, reducing pump flow below **$Q_{\text{crit}} = 0.15\ \text{L/min}$** causes canopy $\text{CO}_2$ drawdown ($|\Delta[\text{CO}_2]| > 16.8\ \text{ppm}$) and vapor accumulation ($\Delta[\text{H}_2\text{O}] > 30.8\ \text{mmol m}^{-3}$) to spike sharply. In the event of a total pump power outage ($Q \le 0.04\ \text{L/min}$), the **Maximum Allowable Outage Grace Window is $t_{\text{grace}} = 10.4\text{--}18.5\ \text{min}$** before transpired water vapor reaches $92\%$ relative humidity condensation onset on the polycarbonate lid optics.

---

## 4. Summary of Lunar Flight Operating Rules

Table 2 synthesizes the recommended LEAF GBC2 pump and fluidic actuator setpoints across all mission phases on the Moon ($0.166\,g$).

**Table 2: Recommended LEAF GBC2 Lunar Surface ($0.166\,g$) Flight Operating Schedule**

| Mission Phase | Timeline | Photoperiod | Recommended Flow Setpoint | Primary Physical / Biological Rationale |
| :--- | :--- | :--- | :---: | :--- |
| **Phase 1: Seed Imbibition & Germination** | Days 0–2 | 16h Day / 8h Night | $Q = 0.25\ \text{L/min}$ | Flushes dark seed respiration ($-0.45\ \text{nmol/s}$) while preventing seed-coat desiccation on exposed substrate top ($Y = 14.62\ \text{mm}$). |
| **Phase 2: Cotyledon & Rosette Establishment** | Days 3–6 | 16h Day / 8h Night | $Q_{\text{day}} = 0.50\ \text{L/min}$<br>$Q_{\text{night}} = 0.40\ \text{L/min}$ | Supports autotrophic transition (Day 2.8–3.2) across small cotyledons ($\text{LAI} < 2.0$) without excessive mechanical wind shear. |
| **Phase 3: Rapid Vegetative Expansion** | Days 7–9 | 16h Day / 8h Night | $Q_{\text{day}} = 0.65\ \text{L/min}$<br>$Q_{\text{night}} = 0.50\ \text{L/min}$ | Counteracts rising Brinkman drag as 32 *B. rapa* plants elongate past $Y = 20\ \text{mm}$ ($\text{LAI} \approx 5.0$). |
| **Phase 4: Full Canopy Maturity & Closure** | Days 10–14 | 16h Day (LED On)<br>8h Night (LED Off) | **$Q_{\text{day}} = 0.85\ \text{L/min}$**<br>**$Q_{\text{night}} = 0.60\ \text{L/min}$** | Compensates $83.4\%$ lunar buoyancy loss ($\text{Ri} = 0.048$), cuts $\Delta[\text{CO}_2]$ & $\Delta[\text{H}_2\text{O}]$ by $37\%$, and prevents nocturnal $\text{CO}_2$ pooling. |
| **Phase 5: Off-Nominal Power Brownout** | Contingency | Any | $Q_{\text{min}} \ge 0.15\ \text{L/min}$<br>($t_{\text{grace}} = 10.4\text{--}18.5\ \text{min}$ at $Q=0$) | Prevents foliar thermal runaway ($\Delta T > +3.5\ ^\circ\text{C}$) and lid optical condensation ($\text{RH} > 92\%$). |
| **Phase 6: End-of-Mission RNAlater Fixation** | Day 14 Harvest | Lights Off / Vent Open | **$Q_{\text{inj}} = 10.0\ \text{mL/min}$**<br>via Port 3 Nozzle ($Y = 9.5\ \text{mm}$) | Saturates 3-part root matrix in $22.6\text{--}31.7\ \text{s}$ and floods full $64.6\ \text{cm}^3$ pot in $6.46\ \text{min}$ at $\Delta P = 0.73\text{--}4.04\ \text{kPa}$ ($< 5\ \text{kPa}$ limit). |

---

## 5. Figure and Dataset Catalog

### Publication Figures (`figures/`)
- **Figure 0A ([`fig_orthographic_projection_annotated.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig_orthographic_projection_annotated.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig_orthographic_projection_annotated.svg)):** Three-view engineering orthographic projection (Top Plan Z–X, Transverse Cutaway X–Y at $Z=0$, Longitudinal Section Z–Y at $X=-5\ \text{mm}$) showing Port 1 Air Inlet, Port 2 Central Exhaust, Port 3 Dipping Nozzle ($Y = 9.5\ \text{mm}$), 3-part porous substrate ($Y = 7.62\text{--}14.62\ \text{mm}$), and 32 *B. rapa* vs. 100 *A. thaliana* canopies.
- **Figure 0B ([`assembly_cross_sections_annotated.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/assembly_cross_sections_annotated.png) & [`assembly_cross_sections_exact.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/assembly_cross_sections_exact.png)):** Exact 3D world-coordinate CAD cross-sections verifying mating of `M09024F-MDL-5` and `M09025F-MDL-2`.
- **Figure 1 ([`fig1_pump_speed_vs_aerodynamics.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig1_pump_speed_vs_aerodynamics.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig1_pump_speed_vs_aerodynamics.svg)):** Parametric pump speed sweep ($Q = 0.05\text{--}4.0\ \text{L/min}$) across Earth $1g$ and Lunar $0.166g$.
- **Figure 2 ([`fig2_velocity_fields_earth_vs_lunar.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig2_velocity_fields_earth_vs_lunar.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig2_velocity_fields_earth_vs_lunar.svg)):** Spatial velocity, $\Delta[\text{CO}_2]$, and $\Delta[\text{H}_2\text{O}]$ midplane slices comparing Earth $1g$ ($0.50\ \text{L/min}$) vs. Lunar $0.166g$ ($0.50$ and $0.85\ \text{L/min}$).
- **Figure 3 ([`fig3_lunar_compensatory_pump_curve.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig3_lunar_compensatory_pump_curve.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig3_lunar_compensatory_pump_curve.svg)):** Lunar compensatory pump sizing curve and Richardson/Péclet convection regime boundaries.
- **Figure 4 ([`fig4_mean_age_of_air_and_chamber_uniformity.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig4_mean_age_of_air_and_chamber_uniformity.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig4_mean_age_of_air_and_chamber_uniformity.svg)):** 3D Mean Age-of-Air ($\theta$) spatial distribution, inter-chamber uniformity ($<1\%$ disparity), and Sandberg conservation validation.
- **Figure 5 ([`fig5_flight_hardware_benchmark_gbc2_vs_nasa.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig5_flight_hardware_benchmark_gbc2_vs_nasa.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig5_flight_hardware_benchmark_gbc2_vs_nasa.svg)):** Comparative spaceflight hardware benchmark (LEAF GBC2 vs. NASA BRIC-19, CARA, and Veggie).
- **Figure 6 ([`fig6_predicted_hormone_carbon_phenotypes.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig6_predicted_hormone_carbon_phenotypes.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig6_predicted_hormone_carbon_phenotypes.svg)):** Predicted transcriptomic and phytohormone (Auxin, ABA, Ethylene) perturbation trajectories across ventilation regimes.
- **Figure 7 ([`fig7_rnalater_perfusion_and_species_dynamics.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig7_rnalater_perfusion_and_species_dynamics.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig7_rnalater_perfusion_and_species_dynamics.svg)):** Two-stage RNAlater perfusion front progression and asymmetric 32 *B. rapa* vs. 100 *A. thaliana* gas exchange telemetry.
- **Figure 8 ([`fig8_diurnal_day_night_photoperiod_inversion.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig8_diurnal_day_night_photoperiod_inversion.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig8_diurnal_day_night_photoperiod_inversion.svg)):** Diurnal (16h Day / 8h Night) photoperiod buoyancy inversion, nocturnal respiratory $\text{CO}_2$ pooling, and lunar night setback envelope.
- **Figure 9 ([`fig9_developmental_growth_timeseries.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig9_developmental_growth_timeseries.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig9_developmental_growth_timeseries.svg)):** 14-day developmental growth time-series (Day 0 seed imbibition to Day 14 harvest closure), tracking $\text{LAI}(t)$, cumulative carbon gain, aerodynamic channelling, adaptive pump ramp $Q(t)$, and root-compaction porosity $\epsilon(t)$.
- **Figure 10 ([`fig10_off_nominal_sensitivity_sweeps.png`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig10_off_nominal_sensitivity_sweeps.png) / [`.svg`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/figures/fig10_off_nominal_sensitivity_sweeps.svg)):** Off-nominal mission sensitivity sweeps across elevated spacecraft cabin $\text{CO}_2$ ($420\text{--}5{,}000\ \text{ppm}$), root porosity compaction ($\epsilon = 0.50 \rightarrow 0.20$), and pump brownout/outage grace windows.

### Quantitative Tabular Datasets (`results/tables/`)
1. [`gbc2_asymmetric_canopy_species_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_asymmetric_canopy_species_results.tsv) — Chamber-resolved gas exchange for 32 *B. rapa* vs. 100 *A. thaliana*.
2. [`gbc2_cfd_sweep_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_cfd_sweep_results.tsv) — 16-condition parametric pump speed sweep at Earth $1g$ and Lunar $0.166g$.
3. [`gbc2_advanced_analytics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_advanced_analytics.tsv) — Mean Age of Air ($\theta$), Sandberg ratios, and hormone perturbation indices.
4. [`gbc2_diurnal_day_night_results.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_diurnal_day_night_results.tsv) — 16h Day vs. 8h Night photoperiod inversion metrics.
5. [`gbc2_rnalater_perfusion_kinetics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_rnalater_perfusion_kinetics.tsv) — Harvest RNAlater Darcy + VOF filling times and pressure drops ($2\text{--}25\ \text{mL/min}$).
6. [`gbc2_hardware_benchmark_comparison.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_hardware_benchmark_comparison.tsv) — LEAF GBC2 vs. NASA BRIC-19, CARA, and Veggie metrics.
7. [`gbc2_covariates_for_omics.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_covariates_for_omics.tsv) — Dimensionless physical microclimate covariates ($\text{Re}, \text{Gr}, \text{Ri}, \text{Pe}, \text{Sh}$) for multi-omics regression.
8. [`gbc2_developmental_growth_timeseries.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_developmental_growth_timeseries.tsv) — 14-day developmental growth time-series across Earth $1g$, Lunar uncompensated, and Lunar adaptive pump schedules.
9. [`gbc2_off_nominal_sensitivity.tsv`](file:///Users/drb_laptop/Documents/LEAF_GBC2_pot_model/results/tables/gbc2_off_nominal_sensitivity.tsv) — Off-nominal mission sensitivity sweeps (cabin $\text{CO}_2$, root compaction, and pump brownout grace windows).
