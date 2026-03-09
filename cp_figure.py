#!/usr/bin/env python3
"""
cp_figure.py — figures for circular permutation analysis of S6 (1RIS).

main figure: three-column layout
  left:   WT vs CP54 chain/arc diagrams (reorganization)
  middle: L2→L3 energetic mismatch (bar chart)
  right:  insulation comparison (four bars)

SI figure: four-row arc diagrams (WT, CP13, CP54, CP68)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from pathlib import Path

RESULTS_DIR = Path("/storage/kiran-stuff/codon_project/cp_results")
N_RES = 97

# ── segment definitions ──
WT_SEGMENTS = [
    (0, 9, 'F0', 'B'), (10, 43, 'F1', 'A'), (44, 71, 'F2', 'A'),
    (72, 85, 'F3', 'B'), (86, 96, 'F4', 'B'),
]

CP_SEGMENTS_WT = {
    54: [
        ([range(54, 80)], 'F0\''), ([range(80, 92)], 'F1\''),
        ([range(0, 9), range(92, 97)], 'F2\''),
        ([range(9, 33)], 'F3\''), ([range(33, 43)], 'F4\''),
        ([range(43, 54)], 'F5\''),
    ],
    13: [
        ([range(13, 37)], 'F0\''), ([range(37, 47)], 'F1\''),
        ([range(47, 57)], 'F2\''), ([range(57, 81)], 'F3\''),
        ([range(81, 97)], 'F4\''), ([range(0, 13)], 'F5\''),
    ],
    68: [
        ([range(68, 80)], 'F0\''), ([range(80, 96)], 'F1\''),
        ([range(0, 14), range(96, 97)], 'F2\''),
        ([range(17, 49)], 'F3\''), ([range(49, 68)], 'F4\''),
    ],
}

SEG_COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377']
L2_COLOR = '#2166AC'
L2_TO_L3_COLOR = '#E66100'
L3_TO_L2_COLOR = '#5D9B3A'
CUT_COLOR = '#CC0000'


def draw_arc(ax, i, j, height_scale=1.0, above=True, **kwargs):
    """draw bezier arc between residues i and j."""
    mid = (i + j) / 2.0
    span = abs(j - i)
    h = np.sqrt(span) * 0.7 * height_scale
    if not above:
        h = -h
    verts = [(i, 0), (i, h*0.5), (mid, h), (j, h*0.5), (j, 0)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CURVE4, MplPath.LINETO]
    path = MplPath(verts, codes)
    patch = mpatches.PathPatch(path, fill=False, **kwargs)
    ax.add_patch(patch)


def draw_chain_wt(ax, y_center, fontsize=6):
    """draw WT chain with segment colors."""
    bar_h = 1.0
    for idx, (start, end, label, stype) in enumerate(WT_SEGMENTS):
        w = end - start + 1
        rect = plt.Rectangle((start, y_center - bar_h/2), w, bar_h,
                              facecolor=SEG_COLORS[idx], edgecolor='black',
                              linewidth=0.5, alpha=0.75, zorder=5)
        ax.add_patch(rect)
        mid = (start + end) / 2
        ax.text(mid, y_center, f'{label}\n({stype})',
                ha='center', va='center', fontsize=fontsize,
                fontweight='bold', zorder=6)


def draw_chain_cp(ax, y_center, cut_site, fontsize=6):
    """draw CP chain segments in WT coordinates."""
    bar_h = 1.0
    segs = CP_SEGMENTS_WT[cut_site]
    for idx, (range_list, label) in enumerate(segs):
        all_res = []
        for r in range_list:
            all_res.extend(list(r))
        if not all_res:
            continue
        all_res.sort()
        blocks = []
        bs = all_res[0]
        prev = bs
        for res in all_res[1:]:
            if res == prev + 1:
                prev = res
            else:
                blocks.append((bs, prev))
                bs = res
                prev = res
        blocks.append((bs, prev))
        for b_start, b_end in blocks:
            w = b_end - b_start + 1
            rect = plt.Rectangle((b_start, y_center - bar_h/2), w, bar_h,
                                  facecolor=SEG_COLORS[idx % len(SEG_COLORS)],
                                  edgecolor='black', linewidth=0.5,
                                  alpha=0.75, zorder=5)
            ax.add_patch(rect)
        largest = max(blocks, key=lambda b: b[1] - b[0])
        mid = (largest[0] + largest[1]) / 2
        ax.text(mid, y_center, label,
                ha='center', va='center', fontsize=fontsize,
                fontweight='bold', zorder=6)


def setup_chain_axis(ax, xlabel=True):
    """format axis for chain diagram."""
    ax.set_xlim(-3, N_RES + 3)
    ax.set_ylim(-7.5, 7.5)
    if xlabel:
        ax.set_xlabel('Residue (WT numbering)', fontsize=8)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks(range(0, N_RES + 1, 20))
    ax.tick_params(axis='x', labelsize=7)


def make_main_figure():
    """
    three-column main text figure:
      left:   WT vs CP54 arc diagrams
      middle: energetic mismatch bar chart
      right:  insulation comparison
    """
    comp = pd.read_csv(RESULTS_DIR / "cp54_comparison.csv")
    wt_contacts = pd.read_csv(RESULTS_DIR / "wt_contacts.csv")
    cp_contacts = pd.read_csv(RESULTS_DIR / "cp54_contacts.csv")

    l2_to_l3 = comp[(comp['wt_layer'] == 'L2') & (comp['cp_layer'] == 'L3')]
    l3_to_l2 = comp[(comp['wt_layer'] == 'L3') & (comp['cp_layer'] == 'L2')]

    fig = plt.figure(figsize=(16, 7.5))

    # ── LEFT COLUMN: arc diagrams (takes ~55% width) ──
    ax_wt = fig.add_axes([0.03, 0.54, 0.52, 0.40])
    ax_cp = fig.add_axes([0.03, 0.08, 0.52, 0.40])

    # WT panel
    ax_wt.set_title('Wild-type S6 — 5 foldons, 28 L2 contacts',
                     fontsize=9, fontweight='bold', loc='left', pad=6)
    draw_chain_wt(ax_wt, 0, fontsize=5.5)
    wt_l2 = wt_contacts[wt_contacts['layer'] == 'L2']
    for _, row in wt_l2.iterrows():
        draw_arc(ax_wt, row['i'], row['j'], height_scale=0.9, above=True,
                 color=L2_COLOR, linewidth=0.9, alpha=0.6)
    for _, row in l2_to_l3.iterrows():
        draw_arc(ax_wt, row['wt_i'], row['wt_j'], height_scale=0.6,
                 above=False, color=L2_TO_L3_COLOR, linewidth=1.0, alpha=0.65)
    ax_wt.axvline(x=54, color=CUT_COLOR, linestyle='--', linewidth=1.2,
                  alpha=0.4, zorder=3)
    ax_wt.text(54, 6.5, 'cut', ha='center', va='bottom',
               fontsize=6.5, color=CUT_COLOR, fontweight='bold')
    # F2 bracket
    ax_wt.annotate('', xy=(44, -4.5), xytext=(71, -4.5),
                   arrowprops=dict(arrowstyle='<->', color=L2_TO_L3_COLOR, lw=1.2))
    ax_wt.text(57.5, -5.3, 'F2: 17 L2, E$_d$=−0.186',
               ha='center', va='top', fontsize=5.5, color=L2_TO_L3_COLOR,
               fontstyle='italic')
    setup_chain_axis(ax_wt, xlabel=False)

    # CP54 panel
    ax_cp.set_title('CP54 — 6 foldons, 11 L2 contacts',
                     fontsize=9, fontweight='bold', loc='left', pad=6)
    draw_chain_cp(ax_cp, 0, 54, fontsize=5.5)
    for _, row in cp_contacts[cp_contacts['layer'] == 'L2'].iterrows():
        wt_i = (int(row['i']) + 54) % N_RES
        wt_j = (int(row['j']) + 54) % N_RES
        if wt_i > wt_j:
            wt_i, wt_j = wt_j, wt_i
        draw_arc(ax_cp, wt_i, wt_j, height_scale=0.9, above=True,
                 color=L2_COLOR, linewidth=0.9, alpha=0.6)
    for _, row in l2_to_l3.iterrows():
        draw_arc(ax_cp, row['wt_i'], row['wt_j'], height_scale=0.6,
                 above=False, color=L2_TO_L3_COLOR, linewidth=1.0, alpha=0.65)
    for _, row in l3_to_l2.iterrows():
        draw_arc(ax_cp, row['wt_i'], row['wt_j'], height_scale=0.6,
                 above=False, color=L3_TO_L2_COLOR, linewidth=1.0, alpha=0.65)
    ax_cp.axvline(x=54, color=CUT_COLOR, linestyle='--', linewidth=1.2,
                  alpha=0.4, zorder=3)
    # wrapping module annotation
    ax_cp.annotate('wrapping\nmodule',
                   xy=(4, -1.0), xytext=(15, -5.0),
                   fontsize=5.5, color=SEG_COLORS[2], fontweight='bold',
                   ha='center', va='top',
                   arrowprops=dict(arrowstyle='->', color=SEG_COLORS[2], lw=1.0))
    setup_chain_axis(ax_cp, xlabel=True)

    # legend for left column
    leg_elements = [
        mpatches.Patch(facecolor=L2_COLOR, alpha=0.6, label='L2 (intra-foldon)'),
        mpatches.Patch(facecolor=L2_TO_L3_COLOR, alpha=0.65,
                       label=f'L2→L3 (n={len(l2_to_l3)})'),
        mpatches.Patch(facecolor=L3_TO_L2_COLOR, alpha=0.65,
                       label=f'L3→L2 (n={len(l3_to_l2)})'),
    ]
    ax_cp.legend(handles=leg_elements, loc='upper right', fontsize=6,
                 frameon=True, edgecolor='gray', fancybox=False)

    # ── MIDDLE COLUMN: energetic mismatch ──
    ax_energy = fig.add_axes([0.60, 0.15, 0.16, 0.75])

    # bar data: mean E_direct for different contact classes
    categories = ['WT L2\n(intra-module)', 'L2→L3\nswitched', 'WT L3\n(interface)']
    values = [
        wt_contacts[wt_contacts['layer'] == 'L2']['E_direct'].mean(),
        l2_to_l3['E_direct'].mean(),
        wt_contacts[wt_contacts['layer'] == 'L3']['E_direct'].mean(),
    ]
    colors_bar = [L2_COLOR, L2_TO_L3_COLOR, '#BBBBBB']

    bars = ax_energy.bar(range(3), values, color=colors_bar, edgecolor='black',
                          linewidth=0.5, width=0.7, alpha=0.8)
    ax_energy.set_xticks(range(3))
    ax_energy.set_xticklabels(categories, fontsize=7)
    ax_energy.set_ylabel('Mean E$_{direct}$ (kJ/mol)', fontsize=8)
    ax_energy.set_title('Energetic mismatch', fontsize=9, fontweight='bold', pad=8)
    ax_energy.axhline(y=0, color='black', linewidth=0.5, linestyle='-')
    ax_energy.spines['top'].set_visible(False)
    ax_energy.spines['right'].set_visible(False)
    ax_energy.tick_params(axis='y', labelsize=7)

    # value labels on bars
    for bar, val in zip(bars, values):
        y_pos = val - 0.015 if val < 0 else val + 0.005
        ax_energy.text(bar.get_x() + bar.get_width()/2, y_pos,
                       f'{val:.3f}', ha='center', va='top' if val < 0 else 'bottom',
                       fontsize=7, fontweight='bold')

    # annotation: arrow between L2→L3 and WT L3
    ax_energy.annotate('', xy=(2, values[2] + 0.01), xytext=(1, values[1] - 0.01),
                       arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
    ax_energy.text(1.5, (values[1] + values[2])/2 + 0.02, 'mismatch',
                   ha='center', va='bottom', fontsize=6.5, color='red',
                   fontweight='bold')

    # ── RIGHT COLUMN: insulation comparison ──
    ax_insul = fig.add_axes([0.82, 0.15, 0.16, 0.75])

    labels = ['WT', 'CP13', 'CP54', 'CP68']
    insulation_values = [0.995, 0.969, 0.929, 0.978]
    bar_colors = ['#555555', '#888888', CUT_COLOR, '#888888']
    bar_alphas = [0.7, 0.5, 0.8, 0.5]

    bars2 = ax_insul.bar(range(4), insulation_values,
                          color=bar_colors, edgecolor='black',
                          linewidth=0.5, width=0.65)
    for bar, alpha in zip(bars2, bar_alphas):
        bar.set_alpha(alpha)

    ax_insul.set_xticks(range(4))
    ax_insul.set_xticklabels(labels, fontsize=7.5)
    ax_insul.set_ylabel('Mean independence (1 − |ρ|)', fontsize=8)
    ax_insul.set_title('Modular insulation', fontsize=9, fontweight='bold', pad=8)
    ax_insul.set_ylim(0.90, 1.005)
    ax_insul.spines['top'].set_visible(False)
    ax_insul.spines['right'].set_visible(False)
    ax_insul.tick_params(axis='y', labelsize=7)

    # value labels
    for bar, val in zip(bars2, insulation_values):
        ax_insul.text(bar.get_x() + bar.get_width()/2, val + 0.002,
                      f'{val:.3f}', ha='center', va='bottom',
                      fontsize=7, fontweight='bold')

    # highlight CP54 as outlier
    ax_insul.annotate('factorizability\nbroken',
                       xy=(2, 0.929), xytext=(2, 0.912),
                       fontsize=6, color=CUT_COLOR, fontweight='bold',
                       ha='center', va='top',
                       arrowprops=dict(arrowstyle='->', color=CUT_COLOR, lw=1.0))

    # panel labels
    fig.text(0.02, 0.96, 'a', fontsize=14, fontweight='bold')
    fig.text(0.58, 0.96, 'b', fontsize=14, fontweight='bold')
    fig.text(0.80, 0.96, 'c', fontsize=14, fontweight='bold')

    for ext in ['png', 'pdf']:
        fig.savefig(RESULTS_DIR / f"cp_main_figure.{ext}",
                    dpi=300, bbox_inches='tight')
    print(f"saved: cp_main_figure.png/pdf")
    plt.close()


def make_si_figure():
    """SI figure: all three CP variants + WT, four-row arc diagrams."""
    wt_contacts = pd.read_csv(RESULTS_DIR / "wt_contacts.csv")
    wt_l2 = wt_contacts[wt_contacts['layer'] == 'L2']

    fig, axes = plt.subplots(4, 1, figsize=(14, 13),
                              gridspec_kw={'hspace': 0.28})

    variants = [
        (None, 'Wild-type S6 (1RIS) — 28 L2 contacts, insulation=0.995'),
        (13, 'CP13 (cut at residue 13) — 6 L2, insulation=0.969'),
        (54, 'CP54 (cut at residue 54) — 11 L2, insulation=0.929'),
        (68, 'CP68 (cut at residue 68) — 9 L2, insulation=0.978'),
    ]

    for idx, (cut_site, title) in enumerate(variants):
        ax = axes[idx]
        label = chr(ord('a') + idx)
        ax.set_title(f'{label}) {title}', fontsize=9, fontweight='bold',
                     loc='left', pad=6)

        if cut_site is None:
            draw_chain_wt(ax, 0, fontsize=5.5)
            for _, row in wt_l2.iterrows():
                draw_arc(ax, row['i'], row['j'], height_scale=0.8, above=True,
                         color=L2_COLOR, linewidth=0.8, alpha=0.55)
        else:
            draw_chain_cp(ax, 0, cut_site, fontsize=5.5)
            cp_contacts = pd.read_csv(RESULTS_DIR / f"cp{cut_site}_contacts.csv")
            comp = pd.read_csv(RESULTS_DIR / f"cp{cut_site}_comparison.csv")

            for _, row in cp_contacts[cp_contacts['layer'] == 'L2'].iterrows():
                wt_i = (int(row['i']) + cut_site) % N_RES
                wt_j = (int(row['j']) + cut_site) % N_RES
                if wt_i > wt_j:
                    wt_i, wt_j = wt_j, wt_i
                draw_arc(ax, wt_i, wt_j, height_scale=0.8, above=True,
                         color=L2_COLOR, linewidth=0.8, alpha=0.55)

            l2l3 = comp[(comp['wt_layer'] == 'L2') & (comp['cp_layer'] == 'L3')]
            for _, row in l2l3.iterrows():
                draw_arc(ax, row['wt_i'], row['wt_j'], height_scale=0.55,
                         above=False, color=L2_TO_L3_COLOR, linewidth=0.9,
                         alpha=0.55)

            l3l2 = comp[(comp['wt_layer'] == 'L3') & (comp['cp_layer'] == 'L2')]
            for _, row in l3l2.iterrows():
                draw_arc(ax, row['wt_i'], row['wt_j'], height_scale=0.55,
                         above=False, color=L3_TO_L2_COLOR, linewidth=0.9,
                         alpha=0.55)

            ax.axvline(x=cut_site, color=CUT_COLOR, linestyle='--',
                       linewidth=1.0, alpha=0.35, zorder=3)

            # stats box
            ax.text(N_RES + 1, 0,
                    f'L2={len(cp_contacts[cp_contacts["layer"]=="L2"])}\n'
                    f'{len(l2l3)} L2→L3\n{len(l3l2)} L3→L2',
                    fontsize=6, va='center', ha='left',
                    bbox=dict(boxstyle='round', facecolor='white',
                              edgecolor='gray', alpha=0.8))

        setup_chain_axis(ax, xlabel=(idx == 3))

    leg = [
        mpatches.Patch(facecolor=L2_COLOR, alpha=0.55, label='L2 (intra-foldon)'),
        mpatches.Patch(facecolor=L2_TO_L3_COLOR, alpha=0.55, label='L2→L3 switched'),
        mpatches.Patch(facecolor=L3_TO_L2_COLOR, alpha=0.55, label='L3→L2 switched'),
        plt.Line2D([0], [0], color=CUT_COLOR, linestyle='--', lw=1.0, label='cut site'),
    ]
    fig.legend(handles=leg, loc='lower center', ncol=4, fontsize=8,
               frameon=True, edgecolor='gray', bbox_to_anchor=(0.5, -0.01))

    for ext in ['png', 'pdf']:
        fig.savefig(RESULTS_DIR / f"cp_all_variants_si.{ext}",
                    dpi=300, bbox_inches='tight')
    print(f"saved: cp_all_variants_si.png/pdf")
    plt.close()


if __name__ == "__main__":
    make_main_figure()
    make_si_figure()
