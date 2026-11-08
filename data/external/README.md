# External solver outputs used by this repository

These tables are copied unchanged from two sibling repositories so that the hardware
comparison (`sim/generate_spaceflight_hardware_comparison.py`, Fig. 18) and the advanced
pipeline's benchmark rows (`sim/run_advanced_pipeline.py`, Fig. 5) read their numbers from
solver output instead of typing them in. Do not edit them here; re-copy from the source.

| File | Source repository and path | Source commit | SHA-256 |
|---|---|---|---|
| `lunarleaf_cfd/T8_enclosure_timescales.csv` | [LunarLeaf-CFD](https://github.com/dr-richard-barker/LunarLeaf-CFD) `results/tables/T8_enclosure_timescales.csv` | `b5c18bd` (2026-10-30) | `c26107f8…bdd49b` |
| `lunarleaf_cfd/T11_photosynthesis_feedback.csv` | LunarLeaf-CFD `results/tables/T11_photosynthesis_feedback.csv` | `b5c18bd` | `7f2c0908…674c07` |
| `lunarleaf_cfd/T19_dish_hardware.csv` | LunarLeaf-CFD `results/tables/T19_dish_hardware.csv` | `b5c18bd` | `604419df…05c52f` |
| `lunarleaf_cfd/T20_dish_feedback_summary.csv` | LunarLeaf-CFD `results/tables/T20_dish_feedback_summary.csv` | `b5c18bd` | `e8e351eb…143b0e` |
| `deepspace_seed_twin/seed_air_speed.csv` | deepspace-seed-twin `results/tables/seed_air_speed.csv` (produced by `scripts/lbm_seed_runs.ts`) | `168e52c` (2026-11-06; local commit, repository not yet public) | `11a65cfd…71c66d` |

Full hashes: `shasum -a 256 data/external/*/*.csv`.

## What each table is, and its limits (as stated by the source repository)

- **T19:** a leaf 1 mm above the agar inside 2-D cross-sections of the BRIC 60 mm round dish
  and the CARA 120 mm square dish, at 1 g and 0 g only (no lunar-gravity run). `g_bl_vs_ambient`
  is conductance against cabin air; `g_bl_vs_enclosure_bulk` is against the dish's own air.
  LunarLeaf compares a **sealed** BRIC dish against its own air, because the whole dish
  atmosphere drifts; its 1 g sealed value (0.330) is marked "indicative only" — not steady, and
  the depleted dish air exceeds the Boussinesq limit (LunarLeaf `analysis/plot_dish_hardware.py`).
- **T8:** analytic sealed-dish (BRIC) atmosphere timescales; CARA and VEGGIE vent to the cabin
  and stay near ambient ("bounded").
- **T20:** sealed BRIC and taped CARA with CO2-limited photosynthesis feedback (LBM).
- **T11:** 0-D closed-loop model, 12 h carbon fixed as % of Earth per enclosure.
- **seed_air_speed:** 2-D lattice-Boltzmann run of the GBC1 Z–Y cross-section in the
  germination phase (no canopy), vendored LunarLeaf solver. Air speed in the first fluid row
  above the substrate (0.144 mm) over the central compartment (|Z| < 7.5 mm), for
  0.25 / 0.5 / 0.85 L/min. No buoyancy source is present, so the same steady flow is written for
  0, 0.165 and 1 g. The slice carries the pump flow per unit of pot width, taken as the 45.09 mm
  lid width matching the CAD geometry. Mass balance is reported as `outflow_over_inflow`.

LunarLeaf's forced-flow (fan/VEGGIE) scenes use a velocity inlet whose delivered flow is under
review, so no VEGGIE conductance from LunarLeaf is used here.
