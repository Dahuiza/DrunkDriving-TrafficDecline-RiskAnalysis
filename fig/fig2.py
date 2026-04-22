"""
FIG2: Event Study and Scenario Analysis Visualization.
Author: Hui Liu
Date: 2026-04-15
Description: This script visualizes event study regression results, comparing the
             short-term and medium-term impacts of the lockdown on drunk driving
             metrics. It includes point estimates, confidence intervals, and
             summary stage coefficients.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import matplotlib
import re
from config import __version, __date

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


# ═══════════════════════════════════════════════════════════════════════
# PATH CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

stage_file   = f'../data/reg_results/{__date}-{__version}/01_main_hetero/pop_All Pool Samples.txt'
pt_file      = f'../data/reg_results/{__date}-{__version}/02_analysis_event_study/pop_All Pool Samples.txt'
save_path    = f'../output/{__date}-{__version}/fig2.png'
save_pdf     = f'../output/{__date}-{__version}/fig2.pdf'

# ═══════════════════════════════════════════════════════════════════════
# FIGURE DIMENSIONS & RESOLUTION
# ═══════════════════════════════════════════════════════════════════════

FIG_WIDTH    = 14 * 8 / 14   # Total width (inches)
FIG_HEIGHT   = 4.5 * 8 / 14  # Total height (inches)
FIG_DPI      = 500

# ═══════════════════════════════════════════════════════════════════════
# SUBPLOT SPACING
# ═══════════════════════════════════════════════════════════════════════

TIGHT_RECT   = [-0.01, 0.05, 1.03, 1.03]   # [left, bottom, right, top]

# ═══════════════════════════════════════════════════════════════════════
# FONT SIZES
# ═══════════════════════════════════════════════════════════════════════

FONT_TITLE   = 10    # Subplot labels a/b/c
FONT_XLABEL  = 10    # X-axis label
FONT_YLABEL  = 10    # Y-axis label
FONT_TICK    = 7.5   # Axis tick numbers
FONT_ANN     = 6.5   # Coefficient labels near error bars
FONT_LEGEND  = 10    # Legend text

# ═══════════════════════════════════════════════════════════════════════
# COLORS
# ═══════════════════════════════════════════════════════════════════════

COLOR_SHORT  = "#E87461"   # Scenario II (Short-term)
COLOR_MEDIUM = "#4F9E8F"   # Scenario III (Medium-term)
COLOR_PRE    = "gray"      # Pre-policy line & CI band

# ═══════════════════════════════════════════════════════════════════════
# LINE STYLES
# ═══════════════════════════════════════════════════════════════════════

LINE_WIDTH   = 1.5    # Main line width
DOT_SIZE     = 8      # Scatter dot size
CI_ALPHA     = 0.25   # Confidence Interval band transparency
REF_LW       = 0.8    # Reference horizontal line width
VLINE_LW     = 0.8    # Vertical line (x=0) width

# ═══════════════════════════════════════════════════════════════════════
# ERRORBAR (STAGE ESTIMATES) STYLES
# ═══════════════════════════════════════════════════════════════════════

EB_MARKERSIZE      = 4     # X-marker size
EB_MARKEREDGEWIDTH = 1     # X-marker line width
EB_LINEWIDTH       = 1.5   # Error bar line width

# ═══════════════════════════════════════════════════════════════════════
# ERRORBAR X-AXIS POSITIONING
# ═══════════════════════════════════════════════════════════════════════

EB_X_SHORT   = 9.8    # X-coordinate for Scenario II errorbar
EB_X_MEDIUM  = 11.8   # X-coordinate for Scenario III errorbar

# ═══════════════════════════════════════════════════════════════════════
# TEXT ANNOTATION OFFSETS
# ═══════════════════════════════════════════════════════════════════════

ANN_OFFSET_SHORT_AB  = (-0.9, -0.5)   # Offset for panels a/b Scenario II
ANN_OFFSET_MEDIUM_AB = (-0.3, -0.5)   # Offset for panels a/b Scenario III
ANN_OFFSET_SHORT_C   = (-0.9, -0.07)  # Offset for panel c Scenario II
ANN_OFFSET_MEDIUM_C  = (-0.3, -0.07)  # Offset for panel c Scenario III

# ═══════════════════════════════════════════════════════════════════════
# Y-AXIS RANGES
# ═══════════════════════════════════════════════════════════════════════

YLIM = {
    "drink": (-0.4, 2.0),
    "acc":   (-0.4, 2.0),
    "ratio": (-0.1, 0.3),
}

# ═══════════════════════════════════════════════════════════════════════
# LEGEND POSITIONING
# ═══════════════════════════════════════════════════════════════════════

LEGEND_ANCHOR = (0.5, -0.03)   # bbox_to_anchor
LEGEND_NCOL   = 6

# ═══════════════════════════════════════════════════════════════════════
# DATA PARSING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def parse_ci(ci_str):
    nums = re.findall(r'-?\d+\.?\d*', ci_str)
    return float(nums[0]), float(nums[1])


def read_result_txt(filepath):
    stop_keywords = {'observations', 'r-squared', 'robust'}
    data = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = [l.rstrip('\r\n') for l in f.readlines()]
    i = 0
    while i < len(lines):
        parts   = lines[i].split('\t')
        varname = parts[0].strip().lower()
        if (not varname or varname in ('variables', '')
                or any(varname.startswith(k) for k in stop_keywords)):
            i += 1
            continue
        if len(parts) >= 5:
            try:
                acc_val   = float(parts[2])
                drink_val = float(parts[3])
                ratio_val = float(parts[4])
            except ValueError:
                i += 1
                continue
            ci_acc = ci_drink = ci_ratio = (np.nan, np.nan)
            j = i + 1
            while j < len(lines):
                ci_parts = lines[j].split('\t')
                if ci_parts[0].strip() == '' and len(ci_parts) >= 5:
                    ci_acc   = parse_ci(ci_parts[2])
                    ci_drink = parse_ci(ci_parts[3])
                    ci_ratio = parse_ci(ci_parts[4])
                    i = j
                    break
                elif ci_parts[0].strip() != '':
                    break
                j += 1
            data[varname] = {
                'acc':   (acc_val,   ci_acc[0],   ci_acc[1]),
                'drink': (drink_val, ci_drink[0], ci_drink[1]),
                'ratio': (ratio_val, ci_ratio[0], ci_ratio[1]),
            }
        i += 1
    return data


def build_series(pt_data, col):
    ref   = 0.0 if col == 'ratio' else 1.0
    order = ['pre_3', 'pre_2', None, 'current',
             'time_1', 'time_2', 'time_3', 'time_4',
             'time_5', 'time_6', 'time_7', 'time_8']
    x_vals, vals, los, his = [], [], [], []
    for xi, key in zip(range(-3, 9), order):
        x_vals.append(xi)
        if key is None:
            vals.append(ref); los.append(ref); his.append(ref)
        else:
            v, lo, hi = pt_data[key][col]
            vals.append(v); los.append(lo); his.append(hi)
    return np.array(x_vals), np.array(vals), np.array(los), np.array(his)


pt_data    = read_result_txt(pt_file)
stage_data = read_result_txt(stage_file)

x_all,  drink_vals, drink_lo, drink_hi = build_series(pt_data, 'drink')
_,      acc_vals,   acc_lo,   acc_hi   = build_series(pt_data, 'acc')
_,      ratio_vals, ratio_lo, ratio_hi = build_series(pt_data, 'ratio')

stage = {
    'drink': {'short':  stage_data['saps_short']['drink'],
              'medium': stage_data['saps_medium']['drink']},
    'acc':   {'short':  stage_data['saps_short']['acc'],
              'medium': stage_data['saps_medium']['acc']},
    'ratio': {'short':  stage_data['saps_short']['ratio'],
              'medium': stage_data['saps_medium']['ratio']},
}

# ═══════════════════════════════════════════════════════════════════════
# PLOT
# ═══════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(FIG_WIDTH, FIG_HEIGHT))

datasets = [
    ("a", "Drunk driving cases",   x_all, drink_vals, drink_lo, drink_hi, 'drink'),
    ("b", "Drunk driving crashes", x_all, acc_vals,   acc_lo,   acc_hi,   'acc'),
    ("c", "The crash incidence",   x_all, ratio_vals, ratio_lo, ratio_hi, 'ratio'),
]

for ax, (panel, title, xs, vals, los, his, key) in zip(axes, datasets):

    ref = 0.0 if key == 'ratio' else 1.0

    # CI band
    ax.fill_between(xs, los, his, color=COLOR_PRE, alpha=CI_ALPHA, zorder=1)

    # Line segments
    seg_defs = [
        ([0, 1, 2, 3],          COLOR_PRE),
        ([3, 4, 5],             COLOR_SHORT),
        ([5, 6, 7, 8, 9, 10, 11], COLOR_MEDIUM),
    ]
    for idx_list, color in seg_defs:
        ax.plot(xs[idx_list], vals[idx_list], color=color,
                linewidth=LINE_WIDTH, zorder=3)

    # Dots
    dot_colors = (
        [COLOR_PRE] * 3
        + [COLOR_SHORT] * 3
        + [COLOR_MEDIUM] * 6
    )
    for xi, vi, dc in zip(xs, vals, dot_colors):
        ax.scatter([xi], [vi], color=dc, s=DOT_SIZE, zorder=4)

    # Reference lines
    ax.hlines(ref, xmin=-4, xmax=10.5, color='black',
              linestyle='--', linewidth=REF_LW, zorder=2)
    ax.axvline(0, color=COLOR_PRE, linestyle='--',
               linewidth=VLINE_LW, zorder=2)

    # Errorbar annotations
    is_c = (panel == 'c')
    ylabel_txt = 'Estimated coefficients' if is_c else 'Estimated coefficients (IRR)'
    ann_off_s  = ANN_OFFSET_SHORT_C  if is_c else ANN_OFFSET_SHORT_AB
    ann_off_m  = ANN_OFFSET_MEDIUM_C if is_c else ANN_OFFSET_MEDIUM_AB

    ax.set_ylabel(ylabel_txt, fontsize=FONT_YLABEL)

    short_v,  short_lo,  short_hi  = stage[key]['short']
    medium_v, medium_lo, medium_hi = stage[key]['medium']

    short_err  = np.array([[short_v  - short_lo],  [short_hi  - short_v]])
    medium_err = np.array([[medium_v - medium_lo], [medium_hi - medium_v]])

    ax.errorbar(EB_X_SHORT, short_v, yerr=short_err,
                fmt='x', color=COLOR_SHORT,
                markersize=EB_MARKERSIZE, markeredgewidth=EB_MARKEREDGEWIDTH,
                elinewidth=EB_LINEWIDTH, capsize=0, zorder=6, clip_on=False)
    ax.annotate(
        f'{short_v:.3f}\n({short_lo:.3f}\n{short_hi:.3f})',
        xy=(EB_X_SHORT, short_v),
        xytext=(EB_X_SHORT + ann_off_s[0], short_v + ann_off_s[1]),
        fontsize=FONT_ANN, ha='left', va='center',
        color='black', annotation_clip=False)

    ax.errorbar(EB_X_MEDIUM, medium_v, yerr=medium_err,
                fmt='x', color=COLOR_MEDIUM,
                markersize=EB_MARKERSIZE, markeredgewidth=EB_MARKEREDGEWIDTH,
                elinewidth=EB_LINEWIDTH, capsize=0, zorder=6, clip_on=False)
    ax.annotate(
        f'{medium_v:.3f}\n({medium_lo:.3f}\n{medium_hi:.3f})',
        xy=(EB_X_MEDIUM, medium_v),
        xytext=(EB_X_MEDIUM + ann_off_m[0], medium_v + ann_off_m[1]),
        fontsize=FONT_ANN, ha='left', va='center',
        color='black', annotation_clip=False)

    # Axes formatting
    ax.set_xlabel('Months since the Lockdown', fontsize=FONT_XLABEL)
    ax.set_title(panel, loc='left', fontsize=FONT_TITLE, fontweight='bold')
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax.tick_params(labelsize=FONT_TICK)
    ax.set_xlim(-4, 8.5)
    ax.set_ylim(YLIM[key])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# ═══════════════════════════════════════════════════════════════════════
# LEGEND
# ═══════════════════════════════════════════════════════════════════════

dot_gray = plt.Line2D([0], [0], marker='o', color=COLOR_PRE,
                      linestyle='-', markersize=5, label='Point estimate (IRR)')
ci_patch = mpatches.Patch(color=COLOR_PRE, alpha=0.4, label='95% CI')
x_red    = plt.Line2D([0], [0], marker='x', color=COLOR_SHORT,
                      linestyle='None', markersize=EB_MARKERSIZE,
                      markeredgewidth=EB_MARKEREDGEWIDTH, label='Scenario II')
x_blue   = plt.Line2D([0], [0], marker='x', color=COLOR_MEDIUM,
                      linestyle='None', markersize=EB_MARKERSIZE,
                      markeredgewidth=EB_MARKEREDGEWIDTH, label='Scenario III')

fig.legend(handles=[dot_gray, ci_patch, x_red, x_blue],
           loc='lower center', ncol=LEGEND_NCOL, fontsize=FONT_LEGEND,
           bbox_to_anchor=LEGEND_ANCHOR, frameon=False,
           fancybox=False, edgecolor='lightgray', facecolor='white')

plt.tight_layout(rect=TIGHT_RECT)
plt.savefig(save_path, dpi=FIG_DPI, bbox_inches='tight')
plt.savefig(save_pdf, dpi=FIG_DPI, bbox_inches='tight')
plt.show()
print("Saved.")
