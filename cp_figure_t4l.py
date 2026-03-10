#!/usr/bin/env python3
"""
cp_figure_t4l.py — figures for T4 lysozyme circular permutation analysis.

main figure: three-column layout
  left:   WT vs CP13 chain/arc diagrams (subdomain decoupling)
  middle: insulation comparison (4 bars)
  right:  insulation vs experimental ΔΔG (quantitative prediction)

SI figure: four-row arc diagrams (WT, CP37, CP13, CP75)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from pathlib import Path

RESULTS_DIR = Path("/storage/kiran-stuff/codon_project/cp_results_t4l")
N_RES = 164

# ── WT segment definitions from LayerScan ──
WT_SEGMENTS = [
    (0, 37, 'F0', 'A'), (38, 47, 'F1', 'B'), (48, 75, 'F2', 'A'),
    (76, 85, 'F3', 'B'), (86, 95, 'F4', 'B'), (96, 105, 'F5', 'B'),
    (106, 131, 'F6', 'A'), (132, 141, 'F7', 'B'), (142, 163, 'F8', 'A'),
]

# CP segments in WT coordinates (from analysis output)
CP_SEGMENTS_WT = {
    13: [
        ([range(13, 37)], "F0'"),      # CP[0-23] = WT[13-36]
        ([range(37, 49)], "F1'"),      # CP[24-35] = WT[37-48]
        ([range(49, 85)], "F2'"),      # CP[36-71] = WT[49-84]
        ([range(85, 103)], "F3'"),     # CP[72-89] = WT[85-102]
        ([range(103, 115)], "F4'"),    # CP[90-101] = WT[103-114]
        ([range(115, 131)], "F5'"),    # CP[102-117] = WT[115-130]
        ([range(131, 147)], "F6'"),    # CP[118-133] = WT[131-146]
        ([range(147, 159)], "F7'"),    # CP[134-145] = WT[147-158]
        ([range(0, 13), range(159, 164)], "F8'"),  # CP[146-163] = WT[0-12 + 159-163] WRAPS
    ],
    37: [
        ([range(37, 47)], "F0'"),      # CP[0-9] = WT[37-46]
        ([range(47, 77)], "F1'"),      # CP[10-39] = WT[47-76]
        ([range(77, 93)], "F2'"),      # CP[40-55] = WT[77-92]
        ([range(93, 109)], "F3'"),     # CP[56-71] = WT[93-108]
        ([range(109, 121)], "F4'"),    # CP[72-83] = WT[109-120]
        ([range(121, 133)], "F5'"),    # CP[84-95] = WT[121-132]
        ([range(133, 147)], "F6'"),    # CP[96-109] = WT[133-146]
        ([range(0, 25), range(147, 164)], "F7'"),  # CP[110-151] = WT[0-24 + 147-163] WRAPS
        ([range(25, 37)], "F8'"),      # CP[152-163] = WT[25-36]
    ],
    75: [
        ([range(75, 85)], "F0'"),      # CP[0-9] = WT[75-84]
        ([range(85, 99)], "F1'"),      # CP[10-23] = WT[85-98]
        ([range(99, 113)], "F2'"),     # CP[24-37] = WT[99-112]
        ([range(113, 129)], "F3'"),    # CP[38-53] = WT[113-128]
        ([range(129, 141)], "F4'"),    # CP[54-65] = WT[129-140]
        ([range(0, 7), range(141, 164)], "F5'"),   # CP[66-95] = WT[0-6 + 141-163] WRAPS
        ([range(7, 35)], "F6'"),       # CP[96-123] = WT[7-34]
        ([range(35, 49)], "F7'"),      # CP[124-137] = WT[35-48]
        ([range(49, 75)], "F8'"),      # CP[138-163] = WT[49-74]
    ],
}

SEG_COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE',
              '#AA3377', '#999933', '#882255', '#44AA99']
L2_COLOR = '#2166AC'
L2_TO_L3_COLOR = '#E66100'
L3_TO_L2_COLOR = '#5D9B3A'
CUT_COLOR = '#CC0000'

# experimental data (from literature)
# CP37: Zhang et al. 1993 — ~0.8 kcal/mol
# CP13: Llinas & Marqusee 1998 — ~3 kcal/mol
# CP75: Llinas & Marqusee 1998 — ~9 kcal/mol
EXP_DDG = {37: 0.8, 13: 3.0, 75: 9.0}
INSULATION = {'WT': 0.9906, 37: 0.9872, 13: 0.9585, 75: 0.9268}


def draw_arc(ax, i, j, height_scale=1.0, above=True, **kwargs):
    """draw bezier arc between residues i and j."""
    mid = (i + j) / 2.0
    span = abs(j - i)
    h = np.sqrt(span) * 0.5 * height_scale
    if not above:
        h = -h
    verts = [(i, 0), (i, h*0.5), (mid, h), (j, h*0.5), (j, 0)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CURVE4, MplPath.LINETO]
    path = MplPath(verts, codes)
    patch = mpatches.PathPatch(path, fill=False, **kwargs)
    ax.add_patch(patch)


def draw_chain_wt(ax, y_center, fontsize=5):
    """draw WT chain with segment colors."""
    bar_h = 1.0
    for idx, (start, end, label, stype) in enumerate(WT_SEGMENTS):
        w = end - start + 1
        rect = plt.Rectangle((start, y_center - bar_h/2), w, bar_h,
                              facecolor=SEG_COLORS[idx], edgecolor='black',
                              linewidth=0.4, alpha=0.75, zorder=5)
        ax.add_patch(rect)
        mid = (start + end) / 2
        # only label segments wide enough
        if w >= 14:
            ax.text(mid, y_center, f'{label}\n({stype})',
                    ha='center', va='center', fontsize=fontsize,
                    fontweight='bold', zorder=6)
        elif w >= 10:
            ax.text(mid, y_center, f'{label}',
                    ha='center', va='center', fontsize=fontsize - 0.5,
                    fontweight='bold', zorder=6)


def draw_chain_cp(ax, y_center, cut_site, fontsize=5):
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
                                  edgecolor='black', linewidth=0.4,
                                  alpha=0.75, zorder=5)
            ax.add_patch(rect)
        largest = max(blocks, key=lambda b: b[1] - b[0])
        w_largest = largest[1] - largest[0] + 1
        mid = (largest[0] + largest[1]) / 2
        if w_largest >= 10:
            ax.text(mid, y_center, label,
                    ha='center', va='center', fontsize=fontsize,
                    fontweight='bold', zorder=6)


def setup_chain_axis(ax, xlabel=True):
    """format axis for chain diagram."""
    ax.set_xlim(-5, N_RES + 5)
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
    three-column figure for T4L:
      left:   WT vs CP13 arc diagrams
      middle: insulation bars (WT + 3 CP variants)
      right:  insulation vs experimental ΔΔG (quantitative)
    """
    wt_contacts = pd.read_csv(RESULTS_DIR / "wt_contacts.csv")
    cp13_contacts = pd.read_csv(RESULTS_DIR / "cp13_contacts.csv")
    comp13 = pd.read_csv(RESULTS_DIR / "cp13_comparison.csv")

    l2_to_l3 = comp13[(comp13['wt_layer'] == 'L2') & (comp13['cp_layer'] == 'L3')]
    l3_to_l2 = comp13[(comp13['wt_layer'] == 'L3') & (comp13['cp_layer'] == 'L2')]

    fig = plt.figure(figsize=(17, 7.5))

    # ── LEFT COLUMN: arc diagrams ──
    ax_wt = fig.add_axes([0.02, 0.54, 0.52, 0.40])
    ax_cp = fig.add_axes([0.02, 0.08, 0.52, 0.40])

    # WT panel
    ax_wt.set_title('Wild-type T4 lysozyme (2LZM) — 9 foldons, 113 L2 contacts',
                     fontsize=9, fontweight='bold', loc='left', pad=6)
    draw_chain_wt(ax_wt, 0)
    wt_l2 = wt_contacts[wt_contacts['layer'] == 'L2']
    for _, row in wt_l2.iterrows():
        draw_arc(ax_wt, row['i'], row['j'], height_scale=0.7, above=True,
                 color=L2_COLOR, linewidth=0.5, alpha=0.4)

    # mark key modules
    # F2 [48-75]: N-subdomain core
    ax_wt.annotate('', xy=(48, -4.5), xytext=(75, -4.5),
                   arrowprops=dict(arrowstyle='<->', color='#228833', lw=1.0))
    ax_wt.text(61.5, -5.3, 'N-sub core (F2)', ha='center', va='top',
               fontsize=5.5, color='#228833', fontstyle='italic')
    # F6 [106-131]: C-subdomain core
    ax_wt.annotate('', xy=(106, -4.5), xytext=(131, -4.5),
                   arrowprops=dict(arrowstyle='<->', color='#AA3377', lw=1.0))
    ax_wt.text(118.5, -5.3, 'C-sub core (F6)', ha='center', va='top',
               fontsize=5.5, color='#AA3377', fontstyle='italic')
    # helix A annotation
    ax_wt.annotate('helix A', xy=(6, 1.0), xytext=(6, 5.5),
                   fontsize=5.5, color='#4477AA', fontweight='bold',
                   ha='center', va='bottom',
                   arrowprops=dict(arrowstyle='->', color='#4477AA', lw=0.8))

    # cut line
    ax_wt.axvline(x=13, color=CUT_COLOR, linestyle='--', linewidth=1.2,
                  alpha=0.4, zorder=3)
    ax_wt.text(15, 6.5, 'cut (CP13)', ha='left', va='bottom',
               fontsize=6.5, color=CUT_COLOR, fontweight='bold')
    setup_chain_axis(ax_wt, xlabel=False)

    # CP13 panel
    ax_cp.set_title('CP13 — 9 foldons, 84 L2 contacts',
                     fontsize=9, fontweight='bold', loc='left', pad=6)
    draw_chain_cp(ax_cp, 0, 13)

    # draw CP L2 arcs
    for _, row in cp13_contacts[cp13_contacts['layer'] == 'L2'].iterrows():
        wt_i = (int(row['i']) + 13) % N_RES
        wt_j = (int(row['j']) + 13) % N_RES
        if wt_i > wt_j:
            wt_i, wt_j = wt_j, wt_i
        draw_arc(ax_cp, wt_i, wt_j, height_scale=0.7, above=True,
                 color=L2_COLOR, linewidth=0.5, alpha=0.4)

    # draw switched contacts below
    for _, row in l2_to_l3.iterrows():
        draw_arc(ax_cp, row['wt_i'], row['wt_j'], height_scale=0.5,
                 above=False, color=L2_TO_L3_COLOR, linewidth=0.6, alpha=0.5)
    for _, row in l3_to_l2.iterrows():
        draw_arc(ax_cp, row['wt_i'], row['wt_j'], height_scale=0.5,
                 above=False, color=L3_TO_L2_COLOR, linewidth=0.6, alpha=0.5)

    ax_cp.axvline(x=13, color=CUT_COLOR, linestyle='--', linewidth=1.2,
                  alpha=0.4, zorder=3)

    # annotate wrapping module
    ax_cp.annotate("F8': helix A +\nC-terminal tail",
                   xy=(6, -0.8), xytext=(25, -5.5),
                   fontsize=5.5, color=SEG_COLORS[8], fontweight='bold',
                   ha='center', va='top',
                   arrowprops=dict(arrowstyle='->', color=SEG_COLORS[8], lw=0.8))

    setup_chain_axis(ax_cp, xlabel=True)

    # legend
    leg_elements = [
        mpatches.Patch(facecolor=L2_COLOR, alpha=0.5, label='L2 (intra-foldon)'),
        mpatches.Patch(facecolor=L2_TO_L3_COLOR, alpha=0.5,
                       label=f'L2→L3 (n={len(l2_to_l3)})'),
        mpatches.Patch(facecolor=L3_TO_L2_COLOR, alpha=0.5,
                       label=f'L3→L2 (n={len(l3_to_l2)})'),
    ]
    ax_cp.legend(handles=leg_elements, loc='upper right', fontsize=6,
                 frameon=True, edgecolor='gray', fancybox=False)

    # ── MIDDLE COLUMN: insulation comparison ──
    ax_insul = fig.add_axes([0.59, 0.15, 0.17, 0.75])

    labels = ['WT', 'CP37', 'CP13', 'CP75']
    insul_vals = [INSULATION['WT'], INSULATION[37], INSULATION[13], INSULATION[75]]
    bar_colors = ['#555555', '#888888', '#CC6600', CUT_COLOR]
    bar_alphas = [0.7, 0.5, 0.7, 0.8]

    bars = ax_insul.bar(range(4), insul_vals,
                         color=bar_colors, edgecolor='black',
                         linewidth=0.5, width=0.65)
    for bar, alpha in zip(bars, bar_alphas):
        bar.set_alpha(alpha)

    ax_insul.set_xticks(range(4))
    ax_insul.set_xticklabels(labels, fontsize=7.5)
    ax_insul.set_ylabel('Mean independence (1 − |ρ|)', fontsize=8)
    ax_insul.set_title('Modular insulation', fontsize=9, fontweight='bold', pad=8)
    ax_insul.set_ylim(0.90, 1.005)
    ax_insul.spines['top'].set_visible(False)
    ax_insul.spines['right'].set_visible(False)
    ax_insul.tick_params(axis='y', labelsize=7)

    for bar, val in zip(bars, insul_vals):
        ax_insul.text(bar.get_x() + bar.get_width()/2, val + 0.002,
                      f'{val:.3f}', ha='center', va='bottom',
                      fontsize=7, fontweight='bold')

    # experimental ΔΔG annotations
    for i, (label, ddg) in enumerate([(None, None), ('0.8', 0.8),
                                       ('3.0', 3.0), ('9.0', 9.0)]):
        if ddg is not None:
            ax_insul.text(i, 0.908,
                          f'ΔΔG={ddg}',
                          ha='center', va='top', fontsize=5.5,
                          color='red', fontstyle='italic')

    # ── RIGHT COLUMN: insulation vs ΔΔG scatter ──
    ax_scatter = fig.add_axes([0.82, 0.15, 0.16, 0.75])

    ddg_vals = [0.8, 3.0, 9.0]
    ins_vals = [INSULATION[37], INSULATION[13], INSULATION[75]]
    labels_sc = ['CP37', 'CP13', 'CP75']
    sc_colors = ['#888888', '#CC6600', CUT_COLOR]

    for x, y, lab, col in zip(ddg_vals, ins_vals, labels_sc, sc_colors):
        ax_scatter.scatter(x, y, c=col, s=80, zorder=5, edgecolors='black',
                           linewidth=0.8)
        ax_scatter.annotate(lab, (x, y), textcoords="offset points",
                            xytext=(8, -3), fontsize=7, fontweight='bold',
                            color=col)

    # add WT point at ΔΔG=0
    ax_scatter.scatter(0, INSULATION['WT'], c='#555555', s=80, zorder=5,
                       edgecolors='black', linewidth=0.8, marker='D')
    ax_scatter.annotate('WT', (0, INSULATION['WT']),
                        textcoords="offset points", xytext=(8, -3),
                        fontsize=7, fontweight='bold', color='#555555')

    # fit line through the 4 points
    all_ddg = [0] + ddg_vals
    all_ins = [INSULATION['WT']] + ins_vals
    z = np.polyfit(all_ddg, all_ins, 1)
    p = np.poly1d(z)
    x_line = np.linspace(-0.5, 10, 50)
    ax_scatter.plot(x_line, p(x_line), 'k--', linewidth=0.8, alpha=0.4)

    # correlation
    from scipy.stats import pearsonr
    r, pval = pearsonr(all_ddg, all_ins)
    ax_scatter.text(0.95, 0.95, f'r = {r:.3f}\np = {pval:.3f}',
                    transform=ax_scatter.transAxes,
                    fontsize=7, ha='right', va='top',
                    bbox=dict(boxstyle='round', facecolor='lightyellow',
                              edgecolor='gray', alpha=0.9))

    ax_scatter.set_xlabel('Experimental ΔΔG (kcal/mol)', fontsize=8)
    ax_scatter.set_ylabel('Insulation (1 − |ρ|)', fontsize=8)
    ax_scatter.set_title('Insulation predicts\nstability loss', fontsize=9,
                          fontweight='bold', pad=8)
    ax_scatter.set_xlim(-1, 10.5)
    ax_scatter.set_ylim(0.90, 1.005)
    ax_scatter.spines['top'].set_visible(False)
    ax_scatter.spines['right'].set_visible(False)
    ax_scatter.tick_params(labelsize=7)

    # panel labels
    fig.text(0.01, 0.96, 'a', fontsize=14, fontweight='bold')
    fig.text(0.57, 0.96, 'b', fontsize=14, fontweight='bold')
    fig.text(0.80, 0.96, 'c', fontsize=14, fontweight='bold')

    for ext in ['png', 'pdf']:
        fig.savefig(RESULTS_DIR / f"t4l_cp_main_figure.{ext}",
                    dpi=300, bbox_inches='tight')
    print(f"saved: t4l_cp_main_figure.png/pdf")
    plt.close()


def make_si_figure():
    """SI figure: all three CP variants + WT, four-row arc diagrams."""
    wt_contacts = pd.read_csv(RESULTS_DIR / "wt_contacts.csv")
    wt_l2 = wt_contacts[wt_contacts['layer'] == 'L2']

    fig, axes = plt.subplots(4, 1, figsize=(16, 14),
                              gridspec_kw={'hspace': 0.28})

    variants = [
        (None, f'Wild-type T4 lysozyme (2LZM) — 113 L2 contacts, '
               f'insulation={INSULATION["WT"]:.3f}'),
        (37, f'CP37 (loop cut, ΔΔG=0.8 kcal/mol) — '
             f'insulation={INSULATION[37]:.3f}'),
        (13, f'CP13 (helix A boundary, ΔΔG=3.0 kcal/mol) — '
             f'insulation={INSULATION[13]:.3f}'),
        (75, f'CP75 (subdomain boundary, ΔΔG=9.0 kcal/mol) — '
             f'insulation={INSULATION[75]:.3f}'),
    ]

    for idx, (cut_site, title) in enumerate(variants):
        ax = axes[idx]
        label = chr(ord('a') + idx)
        ax.set_title(f'{label}) {title}', fontsize=9, fontweight='bold',
                     loc='left', pad=6)

        if cut_site is None:
            draw_chain_wt(ax, 0)
            for _, row in wt_l2.iterrows():
                draw_arc(ax, row['i'], row['j'], height_scale=0.6, above=True,
                         color=L2_COLOR, linewidth=0.5, alpha=0.35)
        else:
            draw_chain_cp(ax, 0, cut_site)
            cp_contacts = pd.read_csv(RESULTS_DIR / f"cp{cut_site}_contacts.csv")
            comp = pd.read_csv(RESULTS_DIR / f"cp{cut_site}_comparison.csv")

            for _, row in cp_contacts[cp_contacts['layer'] == 'L2'].iterrows():
                wt_i = (int(row['i']) + cut_site) % N_RES
                wt_j = (int(row['j']) + cut_site) % N_RES
                if wt_i > wt_j:
                    wt_i, wt_j = wt_j, wt_i
                draw_arc(ax, wt_i, wt_j, height_scale=0.6, above=True,
                         color=L2_COLOR, linewidth=0.5, alpha=0.35)

            l2l3 = comp[(comp['wt_layer'] == 'L2') & (comp['cp_layer'] == 'L3')]
            for _, row in l2l3.iterrows():
                draw_arc(ax, row['wt_i'], row['wt_j'], height_scale=0.45,
                         above=False, color=L2_TO_L3_COLOR, linewidth=0.6,
                         alpha=0.45)

            l3l2 = comp[(comp['wt_layer'] == 'L3') & (comp['cp_layer'] == 'L2')]
            for _, row in l3l2.iterrows():
                draw_arc(ax, row['wt_i'], row['wt_j'], height_scale=0.45,
                         above=False, color=L3_TO_L2_COLOR, linewidth=0.6,
                         alpha=0.45)

            ax.axvline(x=cut_site, color=CUT_COLOR, linestyle='--',
                       linewidth=1.0, alpha=0.35, zorder=3)

            # stats box
            n_l2_cp = len(cp_contacts[cp_contacts['layer'] == 'L2'])
            ax.text(N_RES + 2, 0,
                    f'L2={n_l2_cp}\n'
                    f'{len(l2l3)} L2→L3\n{len(l3l2)} L3→L2',
                    fontsize=6, va='center', ha='left',
                    bbox=dict(boxstyle='round', facecolor='white',
                              edgecolor='gray', alpha=0.8))

        setup_chain_axis(ax, xlabel=(idx == 3))

    leg = [
        mpatches.Patch(facecolor=L2_COLOR, alpha=0.4, label='L2 (intra-foldon)'),
        mpatches.Patch(facecolor=L2_TO_L3_COLOR, alpha=0.5, label='L2→L3 switched'),
        mpatches.Patch(facecolor=L3_TO_L2_COLOR, alpha=0.5, label='L3→L2 switched'),
        plt.Line2D([0], [0], color=CUT_COLOR, linestyle='--', lw=1.0, label='cut site'),
    ]
    fig.legend(handles=leg, loc='lower center', ncol=4, fontsize=8,
               frameon=True, edgecolor='gray', bbox_to_anchor=(0.5, -0.01))

    for ext in ['png', 'pdf']:
        fig.savefig(RESULTS_DIR / f"t4l_cp_all_variants_si.{ext}",
                    dpi=300, bbox_inches='tight')
    print(f"saved: t4l_cp_all_variants_si.png/pdf")
    plt.close()


if __name__ == "__main__":
    make_main_figure()
    make_si_figure()
