"""
FIG3: Robustness Checks Errorbar Plot
Author: Hui Liu
Date: 2026-04-15
Description: This script generates a multi-column forest plot comparing baseline
             regression results with various robustness checks. It visualizes
             Scenario II and III impacts using distinct marker styles and
             shaded CI bands to verify result stability.
"""

import os
import re
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from config import __date, __version

# ══════════════════════════════════════════════════════════════════════════
# PATH CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════
BASELINE_DIR  = f"../data/reg_results/{__date}-{__version}/01_main_hetero"
FILE_DATA_DIR = f"../data/reg_results/{__date}-{__version}/04_robustness_checks"
BASELINE_FILE = "pop_All Pool Samples.txt"
OUT_DIR       = f"../output/{__date}-{__version}"
OUT_FILENAME  = f"fig3.png"
OUT_PDF       = f"fig3.pdf"

# ══════════════════════════════════════════════════════════════════════════
# CANVAS DIMENSIONS & OUTPUT
# ══════════════════════════════════════════════════════════════════════════
OUT_DPI         = 500
OUT_TRANSPARENT = False
FIG_WIDTH       = 7
FIG_HEIGHT      = 8

# ══════════════════════════════════════════════════════════════════════════
# MARGINS & SPACING
# ══════════════════════════════════════════════════════════════════════════
MARGIN_LEFT   = 0.19
MARGIN_RIGHT  = 0.99
MARGIN_TOP    = 0.97
MARGIN_BOTTOM = 0.12
COL_WSPACE    = 0.08

# ══════════════════════════════════════════════════════════════════════════
# FONT SIZES
# ══════════════════════════════════════════════════════════════════════════
FONT_FAMILY     = 'Arial'
FONT_YTICK      = 8.2    # Robustness test names
FONT_XTICK      = 8.0    # X-axis ticks
FONT_XLABEL     = 9.0    # X-axis titles
FONT_LETTER     = 10.0   # Subplot labels (a, b, c)
FONT_GROUPLABEL = 7.5    # Left-side group labels
FONT_LEGEND     = 8.0    # Legend text

# ══════════════════════════════════════════════════════════════════════════
# LABEL POSITIONING (Axes Fraction)
# ══════════════════════════════════════════════════════════════════════════
LETTER_X = -0.05   # 0 = left edge; negative values extend further left
LETTER_Y = 1.03    # 1 = top edge; values > 1 extend further up

# ══════════════════════════════════════════════════════════════════════════
# Y-AXIS LABELS
# ══════════════════════════════════════════════════════════════════════════
SHOW_YTICK_LABELS = True   # True=Display the names of robustness tests
YTICK_PAD         = 0      # Distance between tick labels and the axis line (pt)

# ══════════════════════════════════════════════════════════════════════════
# GROUP LABELS (Left-side Categories)
# ══════════════════════════════════════════════════════════════════════════
SHOW_GROUP_LABELS  = False   # True = Show; False = Hide
GROUP_LABEL_XSHIFT = 0.12    # Offset from x-axis left boundary (data units; larger = further left)
GROUP_LABEL_YSHIFT = 0.0     # Extra vertical offset (positive = downward; negative = upward)
GROUP_LABEL_VA     = 'center'
GROUP_LABEL_HA     = 'right'

# ══════════════════════════════════════════════════════════════════════════
# MARKER & ERROR BAR DIMENSIONS
# ══════════════════════════════════════════════════════════════════════════
# Baseline (Diamond)
BASELINE_MARKERSIZE  = 6.0
BASELINE_LINEWIDTH   = 1.8
BASELINE_CAPSIZE     = 3.5
# Robustness Checks (Circle)
ROB_MARKERSIZE       = 5.0
ROB_LINEWIDTH        = 1.5
ROB_CAPSIZE          = 3.5
# Hollow Scenario III: Overlay circle size
HOLLOW_COVER_SIZE    = 3.8   # Closer to ROB_MARKERSIZE makes the hollow area larger

# Vertical jitter to separate Scenario II and III
JITTER = 0.18

# ══════════════════════════════════════════════════════════════════════════
# COLOR
# ══════════════════════════════════════════════════════════════════════════
SCEN2_BASELINE_COLOR  = '#922B21'   # Scenario II baseline point
SCEN3_BASELINE_COLOR  = '#E74C3C'   # Scenario III baseline point

BASELINE_SHADE2_COLOR = '#FADBD8'   # Scenario II CI Shading (Pink)
BASELINE_SHADE3_COLOR = '#D5E8D4'   # Scenario III CI Shading (Light Green)
BASELINE_SHADE_ALPHA  = 0.35        # Shading transparency; 0=Close/Off

REF_LINE_COLOR  = '#555555'   # Reference dashed line
SEPARATOR_COLOR = '#888888'   # Divider between baseline and robustness tests
GROUP_SEP_COLOR = '#cccccc'   # Divider between groups

# Categorical colors for robustness groups
GROUP_DOT_COLORS = {
    "Weighting":     "#E67E22",
    "Controls":      "#2980B9",
    "Sample":        "#27AE60",
    "Threshold":     "#C0392B",
    "Time window":   "#8E44AD",
    "Fixed effects": "#16A085",
    "SE variants":   "#D4AC0D",
    "Alt outcome":   "#2C3E50",
}
GROUP_BG_COLORS = {k: "#FFFFFF" for k in GROUP_DOT_COLORS}
GROUP_BG_ALPHA  = 0.5

# ══════════════════════════════════════════════════════════════════════════
# DEPENDENT VARIABLE (DV) SETTINGS
# ══════════════════════════════════════════════════════════════════════════
DVS = ['drink', 'acc', 'acc_ratio__']
DV_XLABELS = {
    'acc':         'Estimated Coefficient (IRR)',
    'drink':       'Estimated Coefficient (IRR)',
    'acc_ratio__': 'Estimated Coefficient (OLS)',
}
COL_LETTERS = ['a', 'b', 'c']

# Column indexing in raw text files
DV_COL_INDEX = {
    'acc':         0,
    'drink':       1,
    'acc_ratio__': 2,
}

# ══════════════════════════════════════════════════════════════════════════
# TEST FILE LIST
# ══════════════════════════════════════════════════════════════════════════
ENTRIES = [
    ("robustness_1_weighted.txt",         "Population weighted",      "Weighting"),
    ("robustness_2_ctrl_social.txt",  "Controls: +Social",        "Controls"),
    ("robustness_2_ctrl_weather.txt", "Controls: +weather",       "Controls"),
    ("robustness_2_ctrl_air.txt",     "Controls: +Air Quality",   "Controls"),
    ("robustness_2_ctrl_covid.txt",   "Controls: +COVID",         "Controls"),
    ("robustness_2_ctrl_full.txt",    "Controls: full set",       "Controls"),
    ("robustness_3_drop_top1pct.txt",     "Drop top 1% cities",       "Sample"),
    ("robustness_3_drop_top5pct.txt",     "Drop top 5% cities",       "Sample"),
    ("robustness_3_drop_capitals.txt",    "Drop provincial capital",  "Sample"),
    ("robustness_4_threshold_th5.txt",    "Threshold = 5",            "Threshold"),
    ("robustness_4_threshold_th10.txt",   "Threshold = 10",           "Threshold"),
    ("robustness_4_threshold_th20.txt",   "Threshold = 20",           "Threshold"),
    ("robustness_5_window.txt",           "Window: Shorter baseline", "Time window"),
    ("robustness_6_fe_simple.txt",        "FE: city + year + month",  "Fixed effects"),
    ("robustness_6_fe_prov.txt",          "FE: province \u00d7 year", "Fixed effects"),
    ("robustness_8_cluster_prov.txt",    "Cluster: province",        "SE variants"),
]

# ══════════════════════════════════════════════════════════════════════════
# PLOT
# ══════════════════════════════════════════════════════════════════════════

def parse_txt(filepath):
    with open(filepath, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    result = {}
    current_var = None
    for line in lines:
        cells = line.rstrip('\n').split('\t')
        varname = cells[0].strip()
        if varname in ('saps_short', 'saps_medium'):
            current_var = varname
            try:
                coefs = []
                for c in cells[2:5]:
                    c = c.strip().replace(',', '')
                    coefs.append(float(c) if c not in ('', '.', '-') else None)
            except:
                coefs = [None, None, None]
            while len(coefs) < 3:
                coefs.append(None)
            result[current_var] = {'coef': coefs, 'ci': [None, None, None]}
        elif (current_var
              and len(cells) >= 3
              and cells[0].strip() == ''
              and cells[1].strip() == ''):
            cis = []
            for c in cells[2:5]:
                m = re.search(r'\(?\s*([-\d.]+)\s*[-\u2013]\s*([-\d.]+)\s*\)?', c)
                cis.append((float(m.group(1)), float(m.group(2))) if m else None)
            while len(cis) < 3:
                cis.append(None)
            result[current_var]['ci'] = cis
            current_var = None
    return result


def load_all(entries, data_dir):
    parsed = []
    for fname, label, group in entries:
        fpath = os.path.join(data_dir, fname)
        if os.path.exists(fpath):
            try:
                parsed.append((label, group, parse_txt(fpath)))
            except Exception as e:
                print(f"  [error] {fname}: {e}")
                parsed.append((label, group, None))
        else:
            print(f"  [skip]  {fname}")
            parsed.append((label, group, None))
    return parsed


def plot_solid(ax, coef, lo, hi, y, color, marker, ms, lw, caps, zorder):
    """Render solid markers for Scenario II."""
    if coef is None:
        ax.plot(0, y, 'x', color='#b0b0b0', markersize=6, markeredgewidth=1.2, zorder=zorder)
        return
    ax.errorbar(coef, y, xerr=[[lo], [hi]],
                fmt=marker, color=color, ecolor=color,
                elinewidth=lw, capsize=caps, capthick=lw,
                markersize=ms, zorder=zorder,
                markerfacecolor=color,
                markeredgecolor='white', markeredgewidth=0.5)


def plot_hollow(ax, coef, lo, hi, y, color, marker, ms, lw, caps, zorder):
    """
    Render hollow markers for Scenario III using a white-overlay technique.
    """
    if coef is None:
        ax.plot(0, y, 'x', color='#b0b0b0', markersize=6, markeredgewidth=1.2, zorder=zorder)
        return
    # 1. Draw the error bars (without marker)
    ax.errorbar(coef, y, xerr=[[lo], [hi]],
                fmt='none', color=color, ecolor=color,
                elinewidth=lw, capsize=caps, capthick=lw,
                zorder=zorder)
    # 2. Cover the center with white to create hollow effect
    ax.plot(coef, y, marker=marker,
            color='white', markersize=HOLLOW_COVER_SIZE,
            markeredgecolor='white', zorder=zorder + 1)
    # 3. Draw the colored outline
    ax.plot(coef, y, marker=marker,
            color='none', markersize=ms,
            markeredgecolor=color, markeredgewidth=1.4,
            zorder=zorder + 2)


def draw_column(ax, parsed, baseline_data, dv_idx, dv_name,
                show_yticks, col_letter):
    """
    Generate a single column of the errorbar plot for a specific dependent variable.
    """
    is_ols  = (dv_name == 'acc_ratio__')
    ref     = 0.0 if is_ols else 1.0
    n_rob   = len(parsed)
    n_total = n_rob + 1

    # Extract Baseline Data
    def get_bl(scenario):
        if baseline_data and scenario in baseline_data:
            bc  = baseline_data[scenario]['coef'][dv_idx]
            bci = baseline_data[scenario]['ci'][dv_idx]
            if bc is not None and bci is not None:
                return bc, bci[0], bci[1]
        return None, None, None

    bl2_coef, bl2_lo, bl2_hi = get_bl('saps_short')
    bl3_coef, bl3_lo, bl3_hi = get_bl('saps_medium')

    # Extract Robustness Data
    def get_rob(scenario):
        coefs, lo_errs, hi_errs, groups, labels = [], [], [], [], []
        for label, group, data in parsed:
            labels.append(label)
            groups.append(group)
            if data and scenario in data:
                c  = data[scenario]['coef'][dv_idx]
                ci = data[scenario]['ci'][dv_idx]
                if c is not None and ci is not None:
                    coefs.append(c)
                    lo_errs.append(max(0, c - ci[0]))
                    hi_errs.append(max(0, ci[1] - c))
                else:
                    coefs.append(None); lo_errs.append(0); hi_errs.append(0)
            else:
                coefs.append(None); lo_errs.append(0); hi_errs.append(0)
        return coefs, lo_errs, hi_errs, groups, labels

    c2, lo2, hi2, groups, labels = get_rob('saps_short')
    c3, lo3, hi3, _,      _      = get_rob('saps_medium')

    # Calculate group ranges for labeling
    group_ranges = {}
    for i, grp in enumerate(groups):
        yi = i + 1
        if grp not in group_ranges:
            group_ranges[grp] = [yi, yi]
        else:
            group_ranges[grp][1] = yi

    # Plotting Elements
    # Baseline Shaded Bands
    if bl2_lo is not None and BASELINE_SHADE_ALPHA > 0:
        ax.axvspan(bl2_lo, bl2_hi, facecolor=BASELINE_SHADE2_COLOR,
                   alpha=BASELINE_SHADE_ALPHA, zorder=1, lw=0)
    if bl3_lo is not None and BASELINE_SHADE_ALPHA > 0:
        ax.axvspan(bl3_lo, bl3_hi, facecolor=BASELINE_SHADE3_COLOR,
                   alpha=BASELINE_SHADE_ALPHA, zorder=1, lw=0)

    prev_grp, g_start = None, 0
    for i, grp in enumerate(groups):
        if grp != prev_grp:
            if prev_grp is not None:
                ax.axhspan(g_start + 0.5, i + 0.5,
                           facecolor=GROUP_BG_COLORS.get(prev_grp, '#fff'),
                           alpha=GROUP_BG_ALPHA, zorder=0, lw=0)
            prev_grp, g_start = grp, i
    ax.axhspan(g_start + 0.5, n_rob + 0.5,
               facecolor=GROUP_BG_COLORS.get(prev_grp, '#fff'),
               alpha=GROUP_BG_ALPHA, zorder=0, lw=0)

    # Reference Lines
    ax.axvline(ref, color=REF_LINE_COLOR, lw=1.1, ls='--', zorder=2, alpha=0.65)

    # Divider between Baseline and Robustness checks
    ax.axhline(0.5, color=SEPARATOR_COLOR, lw=1.3, zorder=3)

    # Divider between Groups
    prev_grp = None
    for i, grp in enumerate(groups):
        if grp != prev_grp and i > 0:
            ax.axhline(i + 0.5, color=GROUP_SEP_COLOR, lw=0.7, zorder=2)
        prev_grp = grp

    # Baseline two Scenarios
    if bl2_coef is not None:
        plot_solid(ax, bl2_coef,
                   max(0, bl2_coef - bl2_lo), max(0, bl2_hi - bl2_coef),
                   0 - JITTER, SCEN2_BASELINE_COLOR, 'D',
                   BASELINE_MARKERSIZE, BASELINE_LINEWIDTH, BASELINE_CAPSIZE, zorder=6)
    if bl3_coef is not None:
        plot_hollow(ax, bl3_coef,
                    max(0, bl3_coef - bl3_lo), max(0, bl3_hi - bl3_coef),
                    0 + JITTER, SCEN3_BASELINE_COLOR, 'D',
                    BASELINE_MARKERSIZE, BASELINE_LINEWIDTH, BASELINE_CAPSIZE, zorder=6)

    # Robustness checks for both scenarios
    for i, grp in enumerate(groups):
        color = GROUP_DOT_COLORS.get(grp, '#555555')
        yi = i + 1
        plot_solid(ax, c2[i], lo2[i], hi2[i],
                   yi - JITTER, color, 'o',
                   ROB_MARKERSIZE, ROB_LINEWIDTH, ROB_CAPSIZE, zorder=4)
        plot_hollow(ax, c3[i], lo3[i], hi3[i],
                    yi + JITTER, color, 'o',
                    ROB_MARKERSIZE, ROB_LINEWIDTH, ROB_CAPSIZE, zorder=4)

    # Axis Styling
    ax.set_yticks(range(n_total))
    if show_yticks and SHOW_YTICK_LABELS:
        ax.set_yticklabels(['Baseline'] + labels, fontsize=FONT_YTICK)
        ax.tick_params(axis='y', length=0, pad=YTICK_PAD)
    else:
        ax.set_yticklabels([''] * n_total)
        ax.tick_params(axis='y', length=0)

    ax.invert_yaxis()
    ax.set_ylim(n_total - 0.5, -0.5)
    ax.set_xlabel(DV_XLABELS[dv_name], fontsize=FONT_XLABEL)
    ax.tick_params(axis='x', labelsize=FONT_XTICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#aaa')
    ax.grid(axis='x', lw=0.4, color='#e0e0e0', zorder=0)
    ax.patch.set_facecolor('white')

    # Group Labels
    if show_yticks and SHOW_GROUP_LABELS:
        ax.figure.canvas.draw()
        xmin = ax.get_xlim()[0]
        x_label = xmin - GROUP_LABEL_XSHIFT
        for grp, (y0, y1) in group_ranges.items():
            y_center = (y0 + y1) / 2.0 + GROUP_LABEL_YSHIFT
            ax.text(x_label, y_center, grp,
                    fontsize=FONT_GROUPLABEL,
                    color=GROUP_DOT_COLORS.get(grp, '#555'),
                    fontweight='bold',
                    va=GROUP_LABEL_VA, ha=GROUP_LABEL_HA,
                    clip_on=False, zorder=10)

    # Subplot Letter Label
    ax.text(LETTER_X, LETTER_Y, col_letter,
            transform=ax.transAxes,
            fontsize=FONT_LETTER, fontweight='bold',
            va='top', ha='left', color='#111')


def main():
    plt.rcParams['font.family'] = FONT_FAMILY
    plt.rcParams['axes.unicode_minus'] = False

    bl_path = os.path.join(BASELINE_DIR, BASELINE_FILE)
    baseline_data = None
    if os.path.exists(bl_path):
        baseline_data = parse_txt(bl_path)
        print(f"Baseline loaded: {bl_path}")
    else:
        print(f"[warn] Baseline not found: {bl_path}")

    parsed = load_all(ENTRIES, FILE_DATA_DIR)
    print(f"Robustness files: {sum(1 for _,_,d in parsed if d is not None)}/{len(parsed)}")

    os.makedirs(OUT_DIR, exist_ok=True)

    fig, axes = plt.subplots(
        1, 3, figsize=(FIG_WIDTH, FIG_HEIGHT),
        gridspec_kw={
            'wspace': COL_WSPACE,
            'left':   MARGIN_LEFT,
            'right':  MARGIN_RIGHT,
            'top':    MARGIN_TOP,
            'bottom': MARGIN_BOTTOM,
        }
    )
    fig.patch.set_facecolor('white')

    for j, dv in enumerate(DVS):
        draw_column(axes[j], parsed, baseline_data,
                    DV_COL_INDEX[dv], dv,
                    show_yticks=(j == 0),
                    col_letter=COL_LETTERS[j])

    legend_handles = [
        plt.Line2D([0], [0], marker='D', color='w',
                   markerfacecolor='grey',
                   markeredgecolor='white', markersize=6,
                   label='Scenario II'),
        plt.Line2D([0], [0], marker='D', color='w',
                   markerfacecolor='white',
                   markeredgecolor='grey',
                   markeredgewidth=1.4, markersize=6,
                   label='Scenario III'),
        mpatches.Patch(facecolor=BASELINE_SHADE2_COLOR,
                       alpha=0.6, label='Baseline 95% CI (Scen II)'),
        mpatches.Patch(facecolor=BASELINE_SHADE3_COLOR,
                       alpha=0.6, label='Baseline 95% CI (Scen III)'),
    ]
    fig.legend(handles=legend_handles,
               loc='lower center', ncol=2,
               fontsize=FONT_LEGEND,
               framealpha=0, edgecolor='none',
               bbox_to_anchor=(0.6, 0.0))

    outpath = os.path.join(OUT_DIR, OUT_FILENAME)
    plt.savefig(outpath, dpi=OUT_DPI, bbox_inches='tight',
                transparent=OUT_TRANSPARENT)

    outpdf = os.path.join(OUT_DIR, OUT_PDF)
    plt.savefig(outpdf, dpi=OUT_DPI, bbox_inches='tight',
                transparent=OUT_TRANSPARENT)

    plt.show()
    plt.close()
    print(f"Saved: {outpath}")
    print("Done.")


if __name__ == '__main__':
    main()
