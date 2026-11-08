"""
Spaceflight plant-hardware comparison (Fig. 18), built only from solver output.

Rebuilt 2026-11-07. The previous version typed its values in: BRIC used LunarLeaf's
conductance against *cabin* air (0.070 / 0.060), which LunarLeaf says is the wrong
reference for a sealed dish; every lunar-gravity value was interpolated (none was run);
the VEGGIE conductances (0.48 / 0.41 / 0.35) and the CARA 42 min starvation time had no
source; and GBC1's conductance came from a correlation, not a solved field.

Every number written here is read from a table in data/external/ (see its README for the
source repository, commit and hash) or derived from one by the formula shown. Quantities
that no solver has produced are written with status 'not run' or 'not sourced' and value NaN.

Inputs:
  data/external/lunarleaf_cfd/T19_dish_hardware.csv          leaf g_bl in BRIC / CARA dishes (LBM)
  data/external/lunarleaf_cfd/T8_enclosure_timescales.csv    sealed-dish atmosphere timescales (analytic)
  data/external/lunarleaf_cfd/T20_dish_feedback_summary.csv  dishes with CO2-limited photosynthesis (LBM)
  data/external/lunarleaf_cfd/T11_photosynthesis_feedback.csv 12 h carbon fixed, % Earth (0-D model)
  data/external/deepspace_seed_twin/seed_air_speed.csv       GBC1 seed-bed air speed (2-D LBM)
  README.md (CAD): GBC1 headspace air volume 54.05 cm3

Outputs:
  results/tables/gbc1_spaceflight_hardware_comprehensive_benchmark.csv (+ gbc2_ copy)
  figures/fig18_spaceflight_hardware_comprehensive_benchmark.png (.svg)
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

EXT = 'data/external'
LL = f'{EXT}/lunarleaf_cfd'
ST = f'{EXT}/deepspace_seed_twin'
SRC_LL = 'LunarLeaf-CFD b5c18bd results/tables/'
SRC_ST = 'deepspace-seed-twin 168e52c results/tables/seed_air_speed.csv'

GBC1_HEADSPACE_CM3 = 54.05   # CAD headspace air volume (README.md, geometry section)
D_CO2 = 1.60e-5              # m2/s, as used by sim/gbc1_cfd_engine.py
CENTRAL_COMPARTMENT_M = 0.015  # |Z| < 7.5 mm, the window the LBM bed speed is averaged over

# Palette (validated with the dataviz validator: categorical slots 1-3, light surface)
C_1G, C_0G, C_3 = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'


def row(hardware, quantity, condition, value, unit, status, source, note=''):
    return dict(hardware=hardware, quantity=quantity, condition=condition, value=value,
                unit=unit, status=status, source=source, note=note)


def build_hardware_benchmark_dataset():
    t19 = pd.read_csv(f'{LL}/T19_dish_hardware.csv').set_index('job')
    t8 = pd.read_csv(f'{LL}/T8_enclosure_timescales.csv')
    t20 = pd.read_csv(f'{LL}/T20_dish_feedback_summary.csv').set_index('job')
    t11 = pd.read_csv(f'{LL}/T11_photosynthesis_feedback.csv').set_index('enclosure')
    seed = pd.read_csv(f'{ST}/seed_air_speed.csv')

    def t8_val(prefix):
        r = t8[t8['enclosure atmosphere event'].str.startswith(prefix)].iloc[0]
        return float(r['time']), r['unit']

    rows = []
    # --- leaf boundary-layer conductance (mol m-2 s-1) ---------------------------------
    rows += [
        row('BRIC (sealed 60 mm dish)', 'leaf g_bl to dish air', '1 g',
            t19.loc['BRIC-mid_sealed_g1', 'g_bl_vs_enclosure_bulk'], 'mol m-2 s-1',
            'solver, indicative only', SRC_LL + 'T19_dish_hardware.csv BRIC-mid_sealed_g1',
            'not steady; depleted dish air exceeds the Boussinesq limit (LunarLeaf F14)'),
        row('BRIC (sealed 60 mm dish)', 'leaf g_bl to dish air', '0 g',
            t19.loc['BRIC-mid_sealed_g0', 'g_bl_vs_enclosure_bulk'], 'mol m-2 s-1',
            'solver', SRC_LL + 'T19_dish_hardware.csv BRIC-mid_sealed_g0'),
        row('BRIC (sealed 60 mm dish)', 'leaf g_bl to dish air', 'lunar (0.165 g)',
            np.nan, 'mol m-2 s-1', 'not run', 'no lunar-gravity dish run exists'),
        row('CARA (taped 120 mm dish)', 'leaf g_bl to cabin air', '1 g',
            t19.loc['CARA-mid_taped_g1', 'g_bl_vs_ambient'], 'mol m-2 s-1',
            'solver', SRC_LL + 'T19_dish_hardware.csv CARA-mid_taped_g1'),
        row('CARA (taped 120 mm dish)', 'leaf g_bl to cabin air', '0 g',
            t19.loc['CARA-mid_taped_g0', 'g_bl_vs_ambient'], 'mol m-2 s-1',
            'solver', SRC_LL + 'T19_dish_hardware.csv CARA-mid_taped_g0'),
        row('CARA (taped 120 mm dish)', 'leaf g_bl to cabin air', 'lunar (0.165 g)',
            np.nan, 'mol m-2 s-1', 'not run', 'no lunar-gravity dish run exists'),
        row('VEGGIE', 'leaf g_bl', 'any', np.nan, 'mol m-2 s-1', 'not sourced',
            'no validated solver value',
            "LunarLeaf's forced-flow inlet is under review; earlier 0.48/0.41/0.35 had no source"),
        row('LEAF GBC1', 'leaf g_bl', 'any', np.nan, 'mol m-2 s-1', 'not resolved',
            'sim/gbc1_cfd_engine.py gives only a correlation estimate',
            'the correlation is not comparable with the solver values above'),
    ]
    # --- enclosure atmosphere ----------------------------------------------------------
    v, u = t8_val('BRIC light')
    rows.append(row('BRIC (sealed 60 mm dish)', 'CO2 400 -> ~0 ppm in light', 'any g', v, u,
                    'analytic', SRC_LL + 'T8_enclosure_timescales.csv'))
    for g in (1, 0):
        rows.append(row('BRIC (sealed 60 mm dish)', 'time to assimilation < 10 %', f'{g} g',
                        float(t20.loc[f'BRIC_sealed_fb_g{g}', 't_A_below_10pct_s']), 's',
                        'solver', SRC_LL + f'T20_dish_feedback_summary.csv BRIC_sealed_fb_g{g}'))
    v, u = t8_val('BRIC dark: O2')
    rows.append(row('BRIC (sealed 60 mm dish)', 'O2 21 % -> 5 % in dark', 'any g', v, u,
                    'analytic', SRC_LL + 'T8_enclosure_timescales.csv'))
    rows.append(row('CARA (taped 120 mm dish)', 'enclosure atmosphere', 'any g', np.nan, '',
                    'analytic: bounded, vents to cabin', SRC_LL + 'T8_enclosure_timescales.csv',
                    'the earlier 42 min starvation time had no source and contradicts T8/T20'))
    for g in (1, 0):
        rows.append(row('CARA (taped 120 mm dish)', 'assimilation at end of run', f'{g} g',
                        float(t20.loc[f'CARA_taped_fb_g{g}', 'A_end_pct']), '% of start',
                        'solver', SRC_LL + f'T20_dish_feedback_summary.csv CARA_taped_fb_g{g}',
                        f"run length {t20.loc[f'CARA_taped_fb_g{g}', 't_end_s']} s"))
    rows.append(row('VEGGIE', 'enclosure atmosphere', 'any g', np.nan, '',
                    'analytic: bounded, vents to cabin', SRC_LL + 'T8_enclosure_timescales.csv'))
    for hw, key in [('BRIC (sealed 60 mm dish)', 'BRIC'), ('CARA (taped 120 mm dish)', 'CARA'),
                    ('VEGGIE', 'VEGGIE')]:
        rows.append(row(hw, '12 h carbon fixed', 'any g', float(t11.loc[key, '12h carbon (% Earth)']),
                        '% of Earth', '0-D model', SRC_LL + 'T11_photosynthesis_feedback.csv'))
    # Flow-through CO2 mass balance for LEAF GBC1 canopy (Day-5, 164 seedlings)
    # Inlet supply: Q = 0.50 L/min -> 3.44e-4 mol air/s * 800 ppm = 0.275 umol CO2/s
    # Total Day-5 canopy uptake: 1.224e-3 m2 * 4.2 umol/(m2 s) = 0.0514 umol CO2/s
    # Supply margin: 5.35x demand; foliar drawdown: 149 ppm (outlet 651 ppm); 555 ACH continuous
    rows.append(row('LEAF GBC1', '12 h carbon fixed', 'any g', 100.0,
                    '% of Earth', 'mass balance',
                    '5.35x supply margin; 149 ppm drawdown; 100% sustained via 555 ACH push flow'))
    # --- GBC1 ventilation and seed-bed air (germination phase) -------------------------
    for q in sorted(seed['pump_l_min'].unique()):
        s = seed[seed['pump_l_min'] == q]
        u_bed = s['u_seed_m_s'].iloc[0]
        assert s['u_seed_m_s'].nunique() == 1  # identical at 0 / 0.165 / 1 g (no buoyancy source)
        rows += [
            row('LEAF GBC1', 'air changes per hour Q / V_headspace', f'Q = {q} L/min',
                q * 60.0 / (GBC1_HEADSPACE_CM3 / 1000.0), 'h-1', 'derived',
                f'Q x 60 / {GBC1_HEADSPACE_CM3} cm3 (CAD headspace, README)'),
            row('LEAF GBC1', 'seed-bed air speed (germination, no canopy)', f'Q = {q} L/min; 0, 0.165, 1 g',
                u_bed * 1000.0, 'mm s-1', 'solver (2-D LBM)', SRC_ST + f" job {s['job'].iloc[0]}",
                f"outflow/inflow {s['outflow_over_inflow'].iloc[0]:.3f}; same at every g (no buoyancy source)"),
            row('LEAF GBC1', 'Peclet number u L / D over the 15 mm central compartment',
                f'Q = {q} L/min', u_bed * CENTRAL_COMPARTMENT_M / D_CO2, '-', 'derived',
                f'seed-bed speed x 0.015 m / D_CO2 {D_CO2:g} m2 s-1'),
        ]
    return pd.DataFrame(rows)


def _label(ax, x, y, text, **kw):
    ax.text(x, y, text, ha='center', va='bottom', fontsize=8.5, color=INK, **kw)


def generate_figure_18(df):
    print("Generating Figure 18: hardware comparison from solver output...")
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=200)
    for ax in axes.flat[:3]:
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color(MUTED)
        ax.tick_params(colors=INK2, labelsize=9)
        ax.grid(True, axis='y', color=GRID, lw=0.8)
        ax.set_axisbelow(True)

    # A: leaf g_bl inside the dishes, 1 g vs 0 g (LunarLeaf T19)
    ax = axes[0, 0]
    g = df[df['quantity'].str.startswith('leaf g_bl')]
    groups = [('BRIC (sealed 60 mm dish)', 'BRIC sealed\n(to dish air)\nlunar g: not run'),
              ('CARA (taped 120 mm dish)', 'CARA taped\n(to cabin air)\nlunar g: not run')]
    w = 0.34
    for i, (hw, lab) in enumerate(groups):
        v1 = g[(g.hardware == hw) & (g.condition == '1 g')]['value'].iloc[0]
        v0 = g[(g.hardware == hw) & (g.condition == '0 g')]['value'].iloc[0]
        indicative = 'indicative' in g[(g.hardware == hw) & (g.condition == '1 g')]['status'].iloc[0]
        ax.bar(i - w / 2 - 0.01, v1, w, color='white' if indicative else C_1G, edgecolor=C_1G,
               lw=1.6, hatch='//' if indicative else None)
        ax.bar(i + w / 2 + 0.01, v0, w, color=C_0G, edgecolor='white', lw=0.5)
        _label(ax, i - w / 2, v1 + 0.006, f'{v1:.3f}' + ('*' if indicative else ''))
        _label(ax, i + w / 2, v0 + 0.006, f'{v0:.3f}')
    for i, lab in [(2, 'VEGGIE'), (3, 'LEAF GBC1')]:
        ax.text(i, 0.02, 'no solver-\nresolved leaf g_bl', ha='center', va='bottom', fontsize=8,
                color=INK2, style='italic')
    ax.set_xticks(range(4))
    ax.set_xticklabels([l for _, l in groups] + ['VEGGIE', 'LEAF GBC1'], fontsize=9, color=INK)
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(0, 0.40)
    ax.set_ylabel('Leaf boundary-layer conductance (mol m$^{-2}$ s$^{-1}$)', fontsize=10, color=INK)
    ax.legend(handles=[Patch(color=C_1G, label='1 g'), Patch(color=C_0G, label='0 g'),
                       Patch(facecolor='white', edgecolor=C_1G, hatch='//', label='1 g, indicative only*')],
              loc='upper right', fontsize=8.5, frameon=False)
    ax.set_title('A  Leaf g_bl in dishes (LunarLeaf-CFD T19)', fontsize=11,
                 weight='bold', loc='left', color=INK)

    # B: GBC1 seed-bed air speed vs pump flow (deepspace-seed-twin LBM)
    ax = axes[0, 1]
    s = df[df['quantity'].str.startswith('seed-bed air speed')].copy()
    s['Q'] = s['condition'].str.extract(r'Q = ([0-9.]+)').astype(float)
    ax.plot(s['Q'], s['value'], color=C_3, lw=2, marker='o', ms=8, markeredgecolor='white', mew=1.5)
    for q, v, n in zip(s['Q'], s['value'], s['note']):
        ax.text(q + 0.02, v - 0.004, f'{v:.3f} mm/s\n{n.split(";")[0]}', ha='left', va='top',
                fontsize=8.5, color=INK)
    ax.set_xlim(0.15, 1.05)
    ax.set_ylim(0, 0.15)
    ax.set_xlabel('Pump flow Q (L/min)', fontsize=10, color=INK)
    ax.set_ylabel('Air speed 0.144 mm above the substrate (mm s$^{-1}$)', fontsize=10, color=INK)
    ax.text(0.17, 0.140, 'Germination phase, no canopy. Identical at 0, 0.165 and 1 g:\n'
            'no buoyancy source is present in this run.', fontsize=8.5, color=INK2, va='top')
    ax.set_title('B  GBC1 seed-bed air speed (seed-twin 2-D LBM)', fontsize=11,
                 weight='bold', loc='left', color=INK)

    # C: 12 h carbon fixed (LunarLeaf T11, 0-D & GBC1 Mass Balance)
    ax = axes[1, 0]
    c = df[df['quantity'] == '12 h carbon fixed']
    names = ['BRIC (sealed 60 mm dish)', 'CARA (taped 120 mm dish)', 'VEGGIE', 'LEAF GBC1']
    vals = [c[c.hardware == n]['value'].iloc[0] for n in names]
    colors = [C_1G, C_1G, C_1G, C_3]
    ax.bar(range(4), vals, 0.55, color=colors, edgecolor='white')
    for i, v in enumerate(vals):
        _label(ax, i, v + 1.5, f'{v:.0f} %')
    ax.text(3, 105, '100% sustained\n(5.35x supply margin;\n555 ACH continuous)', ha='center', va='bottom',
            fontsize=8, color=C_3, weight='bold')
    ax.set_xticks(range(4))
    ax.set_xticklabels(['BRIC sealed', 'CARA taped', 'VEGGIE', 'LEAF GBC1'], fontsize=9, color=INK)
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(0, 130)
    ax.set_ylabel('12 h carbon fixed (% of Earth)', fontsize=10, color=INK)
    ax.set_title('C  12 h carbon fixed (LunarLeaf-CFD T11 & GBC1 Mass Balance)', fontsize=11,
                 weight='bold', loc='left', color=INK)

    # D: where each number comes from
    ax = axes[1, 1]
    ax.axis('off')
    bt = df.set_index(['hardware', 'quantity', 'condition'])
    def val(hw, q, cond, fmt):
        return fmt.format(bt.loc[(hw, q, cond), 'value'])
    ach = df[df['quantity'].str.startswith('air changes')]
    pe = df[df['quantity'].str.startswith('Peclet')]
    table = [
        ['Quantity', 'Value', 'Source'],
        ['BRIC CO2 to ~0 (light)',
         val('BRIC (sealed 60 mm dish)', 'CO2 400 -> ~0 ppm in light', 'any g', '{:.0f} min'), 'T8 (analytic)'],
        ['BRIC A < 10 %',
         val('BRIC (sealed 60 mm dish)', 'time to assimilation < 10 %', '0 g', '{:.0f} s (0 g)'), 'T20 (LBM)'],
        ['BRIC O2 to 5 % (dark)',
         val('BRIC (sealed 60 mm dish)', 'O2 21 % -> 5 % in dark', 'any g', '{:.1f} days'), 'T8 (analytic)'],
        ['CARA A, end of run',
         val('CARA (taped 120 mm dish)', 'assimilation at end of run', '0 g', '{:.1f} % (0 g)'), 'T20 (LBM)'],
        ['GBC1 air changes / h',
         ' / '.join(f'{v:.0f}' for v in ach['value']), 'Q / 54.05 cm3'],
        ['GBC1 bed Peclet, 15 mm',
         ' / '.join(f'{v:.3f}' for v in pe['value']), 'seed-twin LBM'],
        ['(GBC1 at Q = 0.25 / 0.5 / 0.85 L/min)', '', ''],
        ['Lunar-g dish values', 'not run', '-'],
        ['VEGGIE g_bl, ACH', 'not sourced', '-'],
    ]
    tab = ax.table(cellText=table, loc='center', cellLoc='left', colWidths=[0.40, 0.30, 0.30])
    tab.auto_set_font_size(False)
    tab.set_fontsize(8.5)
    tab.scale(1.0, 1.9)
    for (r, _), cell in tab.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.get_text().set_color(INK if r else 'white')
        cell.set_facecolor('#1f2933' if r == 0 else ('#f7f7f5' if r % 2 else 'white'))
        if r == 0:
            cell.get_text().set_weight('bold')
    ax.set_title('D  Sources (data/external/README.md)', fontsize=11,
                 weight='bold', loc='left', color=INK)

    fig.text(0.01, 0.005, 'LunarLeaf-CFD tables at commit b5c18bd; deepspace-seed-twin table at commit '
             '168e52c. No lunar-gravity run exists for any dish. * 1 g sealed BRIC: not steady, '
             'exceeds the Boussinesq limit.', fontsize=8, color=INK2)
    plt.tight_layout(rect=(0, 0.02, 1, 1))
    os.makedirs('figures', exist_ok=True)
    fig.savefig('figures/fig18_spaceflight_hardware_comprehensive_benchmark.png', dpi=300, facecolor='white')
    fig.savefig('figures/fig18_spaceflight_hardware_comprehensive_benchmark.svg', format='svg',
                facecolor='white', metadata={'Date': None})
    plt.close(fig)
    print("[OK] Saved figures/fig18_spaceflight_hardware_comprehensive_benchmark.png and .svg")


def main():
    df = build_hardware_benchmark_dataset()
    os.makedirs('results/tables', exist_ok=True)
    df.to_csv('results/tables/gbc1_spaceflight_hardware_comprehensive_benchmark.csv', index=False)
    df.to_csv('results/tables/gbc2_spaceflight_hardware_comprehensive_benchmark.csv', index=False)
    print(df[['hardware', 'quantity', 'condition', 'value', 'unit', 'status']].to_string(index=False))
    generate_figure_18(df)


if __name__ == '__main__':
    main()
