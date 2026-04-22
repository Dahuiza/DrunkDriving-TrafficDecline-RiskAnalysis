"""
SI-FIG2: Robustness Checks - Placebo Tests and Spurious Effects Analysis
Author: Hui Liu
Date: 2026-04-15
Description: This script visualizes robustness checks for the main regression
             results. It includes two primary components:
             1. Panels (a-c): Probability density of coefficients from 1,000
                randomized placebo simulations (Scenario II vs. III),
                comparing them against the true point estimates.
             2. Panels (d-f): Event study placebo tests using randomized
                treatment timings to validate the identification strategy
                against pre-existing trends or spurious correlations.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import re
import os
import matplotlib
from config import __date, __version

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family']     = "sans-serif"


# ═══════════════════════════════════════════════════════════════════════
# File path configurations
# ═══════════════════════════════════════════════════════════════════════

txt_path     = f'../data/reg_results/{__date}-{__version}/01_main_hetero/pop_All Pool Samples.txt'
placebo_did_path = f'../data/reg_results/{__date}-{__version}/03_placebo_test/placebo_for_python.xlsx'
placebo_es_path = f'../data/reg_results/{__date}-{__version}/03_placebo_test/placebo_event_study.xlsx'
pt_file = f'../data/reg_results/{__date}-{__version}/02_analysis_event_study/pop_All Pool Samples.txt'
save_path = f'../output/{__date}-{__version}/SI_fig2.png'
save_path_pdf = f'../output/{__date}-{__version}/SI-fig2.pdf'


# ═══════════════════════════════════════════════════════════════════════
# Figure Dimensions and Resolution
# ═══════════════════════════════════════════════════════════════════════

FIG_WIDTH    = 7
FIG_HEIGHT   = 5.5
FIG_DPI      = 500
HSPACE       = 0.45


# ═══════════════════════════════════════════════════════════════════════
# Subplot Layout Margins
# ═══════════════════════════════════════════════════════════════════════

SUBPLOT_LEFT   = 0.065
SUBPLOT_RIGHT  = 0.99
SUBPLOT_BOTTOM = 0.10
SUBPLOT_TOP    = 0.97
SUBPLOT_WSPACE = 0.25


# ═══════════════════════════════════════════════════════════════════════
# Font Size Configuration
# ═══════════════════════════════════════════════════════════════════════

FONT_TITLE   = 10
FONT_XLABEL  = 8
FONT_YLABEL  = 8
FONT_TICK    = 6
FONT_ANNOT   = 6
FONT_LEGEND  = 8


# ═══════════════════════════════════════════════════════════════════════
# Legend Positioning (Figure coordinates relative to subplot rows)
# ═══════════════════════════════════════════════════════════════════════

LEGEND_ABC_Y_OFFSET = 0.12   # Offset for the top row (a, b, c)
LEGEND_DEF_Y_OFFSET = 0.12   # Offset for the bottom row (d, e, f)


# ═══════════════════════════════════════════════════════════════════════
# Parameters for Top Row (a-c): Histogram Plots
# ═══════════════════════════════════════════════════════════════════════

ANNOT_Y_HI       = 0.88
ANNOT_Y_LO       = 0.88
ANNOT_X_PAD_L    = 0.02
ANNOT_X_PAD_R    = 0.02
ANNOT_LINE_SPACE = 1.3
Y_TOP_SCALE      = 1.60

HIST_BINS    = 51
HIST_COLOR_S = "#8FA8D4"      # Scenario II color
HIST_COLOR_M = "#F5C884"      # Scenario III color
HIST_ALPHA   = 0.5

VLINE_COLOR_S = "#4472C4"     # Real estimate line (Scenario II)
VLINE_COLOR_M = "#C07020"     # Real estimate line (Scenario III)
VLINE_WIDTH   = 1.2
VLINE_STYLE   = "--"

X_MARGINS = {
    "acc":   (0.7, 0.2),
    "drink": (0.7, 0.2),
    "ratio": (0.6, 0.6),
}

# ═══════════════════════════════════════════════════════════════════════
# Parameters for Bottom Row (d-f): Event Study Plots
# ═══════════════════════════════════════════════════════════════════════

ES_COLOR_PLACEBO  = "#F5C884"
ES_COLOR_CI       = "#F5C884"
ES_CI_ALPHA       = 0.35
ES_COLOR_REAL     = "#4472C4"
ES_COLOR_REAL_CI  = "lightgray"
ES_REAL_CI_ALPHA  = 0.4

ES_LINE_WIDTH     = 1.2
ES_DOT_SIZE       = 8
ES_REF_LW         = 0.7
ES_VLINE_LW       = 0.7

ES_YLIM = {
    "acc":   (-0.4, 2.0),
    "drink": (-0.4, 2.0),
    "ratio": (-0.10, 0.20),
}

ES_YLABEL = {
    "acc":   "Estimated coefficients (IRR)",
    "drink": "Estimated coefficients (IRR)",
    "ratio": "Estimated coefficients",
}


# ═══════════════════════════════════════════════════════════════════════
# Data Parsing Utilities
# ═══════════════════════════════════════════════════════════════════════

def parse_main_results(txt_path):
    with open(txt_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    rows = {}
    for i, line in enumerate(lines):
        if line.strip().startswith("saps_short"):
            rows["short"]  = {"coef_line": line, "ci_line": lines[i + 1]}
        elif line.strip().startswith("saps_medium"):
            rows["medium"] = {"coef_line": line, "ci_line": lines[i + 1]}
    var_order = ["acc", "drink", "ratio"]
    real = {v: {} for v in var_order}
    for term, data in rows.items():
        coefs = re.findall(r"-?\d+\.\d+", data["coef_line"])
        cis   = re.findall(r"\((-?\d+\.\d+)\s*-\s*(-?\d+\.\d+)\)", data["ci_line"])
        for idx, var in enumerate(var_order):
            real[var][term] = {
                "irr":     float(coefs[idx]),
                "ci_low":  float(cis[idx][0]),
                "ci_high": float(cis[idx][1]),
            }
    return real


def parse_ci(ci_str):
    """Helper to parse '(low, high)' confidence interval strings."""
    nums = re.findall(r'-?\d+\.?\d*', ci_str)
    return float(nums[0]), float(nums[1])


def read_event_study_txt(filepath):
    """Parses event study output into a dictionary of coefficients and CIs."""
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
                cp = lines[j].split('\t')
                if cp[0].strip() == '' and len(cp) >= 5:
                    ci_acc   = parse_ci(cp[2])
                    ci_drink = parse_ci(cp[3])
                    ci_ratio = parse_ci(cp[4])
                    i = j
                    break
                elif cp[0].strip() != '':
                    break
                j += 1
            data[varname] = {
                'acc':   (acc_val,   ci_acc[0],   ci_acc[1]),
                'drink': (drink_val, ci_drink[0], ci_drink[1]),
                'ratio': (ratio_val, ci_ratio[0], ci_ratio[1]),
            }
        i += 1
    return data


def build_real_series(pt_data, col):
    """Constructs time-series arrays for event study plotting."""
    ref   = 0.0 if col == 'ratio' else 1.0
    order = ['pre_3', 'pre_2', None,
             'current', 'time_1', 'time_2', 'time_3', 'time_4',
             'time_5',  'time_6', 'time_7', 'time_8']
    xs, vals, los, his = [], [], [], []
    for xi, key in zip(range(-3, 9), order):
        xs.append(xi)
        if key is None:
            vals.append(ref); los.append(ref); his.append(ref)
        else:
            v, lo, hi = pt_data[key][col]
            vals.append(v); los.append(lo); his.append(hi)
    return np.array(xs), np.array(vals), np.array(los), np.array(his)


# ═══════════════════════════════════════════════════════════════════════
# Data Loading Execution
# ═══════════════════════════════════════════════════════════════════════

real   = parse_main_results(txt_path)
df_did = pd.read_excel(placebo_did_path)
df_did.columns = ["iteration",
                  "p_s_acc", "p_m_acc",
                  "p_s_drink", "p_m_drink",
                  "p_s_ratio", "p_m_ratio"]

df_es = pd.read_excel(placebo_es_path)
df_es.columns = [
    "time_point",
    "mean_irr_acc",  "mean_irr_drink",  "mean_ratio",
    "p5_irr_acc",    "p5_irr_drink",    "p5_ratio",
    "p95_irr_acc",   "p95_irr_drink",   "p95_ratio",
]
df_es = df_es.sort_values("time_point").reset_index(drop=True)

pt_data = read_event_study_txt(pt_file)
x_all, drink_vals, drink_lo, drink_hi = build_real_series(pt_data, 'drink')
_,     acc_vals,   acc_lo,   acc_hi   = build_real_series(pt_data, 'acc')
_,     ratio_vals, ratio_lo, ratio_hi = build_real_series(pt_data, 'ratio')


# ═══════════════════════════════════════════════════════════════════════
# Plotting: Row 1 (Panels a-c) - DID Placebo Histograms
# ═══════════════════════════════════════════════════════════════════════

def plot_did_placebo(ax, short_vals, medium_vals, real_s, real_m,
                     xlabel, title, x_margin_left, x_margin_right):
    short_vals  = short_vals.dropna()
    medium_vals = medium_vals.dropna()

    all_vals   = pd.concat([short_vals, medium_vals])
    data_min   = min(all_vals.min(), real_s["irr"], real_m["irr"])
    data_max   = max(all_vals.max(), real_s["irr"], real_m["irr"])
    data_range = data_max - data_min
    xmin   = data_min - data_range * x_margin_left
    xmax   = data_max + data_range * x_margin_right
    xrange = xmax - xmin

    bins = np.linspace(xmin, xmax, HIST_BINS)

    ax.hist(short_vals,  bins=bins,
            weights=np.ones(len(short_vals))  / len(short_vals),
            color=HIST_COLOR_S, alpha=HIST_ALPHA, zorder=2)
    ax.hist(medium_vals, bins=bins,
            weights=np.ones(len(medium_vals)) / len(medium_vals),
            color=HIST_COLOR_M, alpha=HIST_ALPHA, zorder=1)

    ax.axvline(real_s["irr"], color=VLINE_COLOR_S, linewidth=VLINE_WIDTH,
               linestyle=VLINE_STYLE, zorder=4)
    ax.axvline(real_m["irr"], color=VLINE_COLOR_M, linewidth=VLINE_WIDTH,
               linestyle=VLINE_STYLE, zorder=4)

    ax.relim(); ax.autoscale_view()
    y_top = ax.get_ylim()[1] * Y_TOP_SCALE
    ax.set_ylim(0, y_top)

    ci_s = f"({real_s['ci_low']:.3f},{real_s['ci_high']:.3f})"
    ci_m = f"({real_m['ci_low']:.3f},{real_m['ci_high']:.3f})"

    if real_s["irr"] <= real_m["irr"]:
        lv, lci, lc, ll = real_s["irr"], ci_s, VLINE_COLOR_S, "Scenario II"
        rv, rci, rc, rl = real_m["irr"], ci_m, VLINE_COLOR_M, "Scenario III"
    else:
        lv, lci, lc, ll = real_m["irr"], ci_m, VLINE_COLOR_M, "Scenario III"
        rv, rci, rc, rl = real_s["irr"], ci_s, VLINE_COLOR_S, "Scenario II"

    ax.text(lv - ANNOT_X_PAD_L * xrange, y_top * ANNOT_Y_HI,
            f"Real estimate:\n{ll}\n{lv:.3f}\n{lci}",
            ha="right", va="top", fontsize=FONT_ANNOT, color=lc,
            linespacing=ANNOT_LINE_SPACE,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    ax.text(rv + ANNOT_X_PAD_R * xrange, y_top * ANNOT_Y_LO,
            f"Real estimate:\n{rl}\n{rv:.3f}\n{rci}",
            ha="left", va="top", fontsize=FONT_ANNOT, color=rc,
            linespacing=ANNOT_LINE_SPACE,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(xlabel, fontsize=FONT_XLABEL)
    ax.set_ylabel("Share of estimates", fontsize=FONT_YLABEL)
    ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold", loc="left", pad=6)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.2f}"))
    ax.tick_params(labelsize=FONT_TICK)
    ax.spines[["top", "right"]].set_visible(False)


# ═══════════════════════════════════════════════════════════════════════
# Plotting: Row 2 (Panels d-f) - Event Study Placebo Tests
# ═══════════════════════════════════════════════════════════════════════

def plot_es_placebo(ax, var, xs_real, vals_real, lo_real, hi_real, title):
    ref = 0.0 if var == "ratio" else 1.0

    if var == "ratio":
        col_mean, col_p5, col_p95 = "mean_ratio", "p5_ratio", "p95_ratio"
    else:
        col_mean = f"mean_irr_{var}"
        col_p5   = f"p5_irr_{var}"
        col_p95  = f"p95_irr_{var}"

    xs_es   = df_es["time_point"].values
    mean_es = df_es[col_mean].values
    p5_es   = df_es[col_p5].values
    p95_es  = df_es[col_p95].values

    ax.fill_between(xs_es, p5_es, p95_es,
                    color=ES_COLOR_CI, alpha=ES_CI_ALPHA, zorder=1)
    ax.plot(xs_es, mean_es,
            color=ES_COLOR_PLACEBO, linewidth=ES_LINE_WIDTH,
            marker='o', markersize=np.sqrt(ES_DOT_SIZE),
            zorder=3)

    ax.fill_between(xs_real, lo_real, hi_real,
                    color=ES_COLOR_REAL_CI, alpha=ES_REAL_CI_ALPHA, zorder=2)
    ax.plot(xs_real, vals_real,
            color=ES_COLOR_REAL, linewidth=ES_LINE_WIDTH,
            marker='o', markersize=np.sqrt(ES_DOT_SIZE),
            zorder=4)

    ax.axhline(ref, color='black', linestyle='--',
               linewidth=ES_REF_LW, zorder=2)
    ax.axvline(-1, color='gray', linestyle='--',
               linewidth=ES_VLINE_LW, zorder=2)

    if ES_YLIM[var] is not None:
        ax.set_ylim(ES_YLIM[var])

    ax.set_xlim(-4, 9)
    ax.set_xticks(range(-3, 9))
    ax.set_xlabel("Months relative to fake treatment", fontsize=FONT_XLABEL)
    ax.set_ylabel(ES_YLABEL[var], fontsize=FONT_YLABEL)
    ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold", loc="left", pad=6)
    ax.tick_params(labelsize=FONT_TICK)
    ax.spines[["top", "right"]].set_visible(False)


# ═══════════════════════════════════════════════════════════════════════
# Figure Assembly (2x3 Grid)
# ═══════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 3, figsize=(FIG_WIDTH, FIG_HEIGHT))
fig.subplots_adjust(
    left=SUBPLOT_LEFT, right=SUBPLOT_RIGHT,
    bottom=SUBPLOT_BOTTOM, top=SUBPLOT_TOP,
    wspace=SUBPLOT_WSPACE, hspace=HSPACE,
)

# Render Row 1: abc (Histograms)
abc_panels = [
    ("acc",   "p_s_acc",   "p_m_acc",
     "Estimated coefficients (IRR)",      "a"),
    ("drink", "p_s_drink", "p_m_drink",
     "Estimated coefficients (IRR)",      "b"),
    ("ratio", "p_s_ratio", "p_m_ratio",
     "Estimated coefficients", "c"),
]
for ax, (var, col_s, col_m, xlabel, title) in zip(axes[0], abc_panels):
    ml, mr = X_MARGINS[var]
    plot_did_placebo(
        ax,
        short_vals     = df_did[col_s],
        medium_vals    = df_did[col_m],
        real_s         = real[var]["short"],
        real_m         = real[var]["medium"],
        xlabel         = xlabel,
        title          = title,
        x_margin_left  = ml,
        x_margin_right = mr,
    )

# Render Row 2: def (Event Studies)
real_es = {
    "acc":   (x_all, acc_vals,   acc_lo,   acc_hi),
    "drink": (x_all, drink_vals, drink_lo, drink_hi),
    "ratio": (x_all, ratio_vals, ratio_lo, ratio_hi),
}
def_panels = [("acc", "d"), ("drink", "e"), ("ratio", "f")]
for ax, (var, title) in zip(axes[1], def_panels):
    xs, vs, lo, hi = real_es[var]
    plot_es_placebo(ax, var, xs, vs, lo, hi, title)

# ═══════════════════════════════════════════════════════════════════════
# Legend Configuration: Distinct legends for each row
# ═══════════════════════════════════════════════════════════════════════

# rigger layout calculation for accurate coordinate retrieval
fig.canvas.draw()

pos_tl = axes[0][0].get_position()
pos_tr = axes[0][2].get_position()
pos_bl = axes[1][0].get_position()
pos_br = axes[1][2].get_position()

# Legend for panels abc
patch_s = mpatches.Patch(color=HIST_COLOR_S, alpha=HIST_ALPHA,
                          label="Randomized: Scenario II")
patch_m = mpatches.Patch(color=HIST_COLOR_M, alpha=HIST_ALPHA,
                          label="Randomized: Scenario III")
fig.legend(
    handles=[patch_s, patch_m],
    loc="lower center", ncol=2, fontsize=FONT_LEGEND,
    bbox_to_anchor=((pos_tl.x0 + pos_tr.x1) / 2,
                    pos_tl.y0 - LEGEND_ABC_Y_OFFSET),
    frameon=False)

# Legend for panels def
line_p = mlines.Line2D([], [], color=ES_COLOR_PLACEBO, marker='o',
                        markersize=4, linewidth=ES_LINE_WIDTH,
                        label="Placebo mean (±90% CI)")
line_r = mlines.Line2D([], [], color=ES_COLOR_REAL, marker='o',
                        markersize=4, linewidth=ES_LINE_WIDTH,
                        label="Real estimate (95% CI)")
fig.legend(
    handles=[line_p, line_r],
    loc="lower center", ncol=2, fontsize=FONT_LEGEND,
    bbox_to_anchor=((pos_bl.x0 + pos_br.x1) / 2,
                    pos_bl.y0 - LEGEND_DEF_Y_OFFSET),
    frameon=False)

plt.savefig(save_path, dpi=FIG_DPI, bbox_inches="tight")
plt.savefig(save_path_pdf, dpi=FIG_DPI, bbox_inches="tight")
plt.show()
print(f"✓ Figure successfully saved to {save_path}")
