"""
FIG4: Heterogeneity Analysis
Author: Hui Liu
Date: 2026-04-15
Description: This script generates a comprehensive multi-panel figure visualizing the
             heterogeneous dynamic effects of policy interventions
             across multiple dimensions (e.g., Vehicle Type, Road Type, Region).
             The event study framework tests for parallel trends in the pre-treatment
             period and explores variations in the impact across different
             sub-populations and environmental contexts.
"""

import os
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines
from config import __date, __version

# ═══════════════════════════════════════════════════════════════════════
# File Paths
# ═══════════════════════════════════════════════════════════════════════

DATA_DIR  = f"../data/reg_results/{__date}-{__version}/02_analysis_event_study"
SAVE_PATH = f"../output/{__date}-{__version}/SI_fig4.jpg"
SAVE_PATH_PDF = f"../output/{__date}-{__version}/SI_fig4.pdf"


# ═══════════════════════════════════════════════════════════════════════
# Font Settings
# ═══════════════════════════════════════════════════════════════════════

FONT_FAMILY       = "Arial"
FONT_SIZE_TICK    = 9        # 坐标轴刻度
FONT_SIZE_LABEL   = 9        # 坐标轴标签（左侧 & 底部）
FONT_SIZE_PANEL   = 11       # 子图序号（a/b/c...）
FONT_WEIGHT_PANEL = "bold"


# ═══════════════════════════════════════════════════════════════════════
# Legend Configuration (Bottom of Figure)
# ═══════════════════════════════════════════════════════════════════════

FONT_SIZE_LEGEND_CAT  = 7.5    # Font size for category labels
FONT_SIZE_LEGEND_ITEM = 7.0
LEGEND_TOP_MARGIN     = 0.012  # Vertical position of bottom row (figure coords)
GAP_LEGEND = 0.08              # Global horizontal offset for legend alignment


# ═══════════════════════════════════════════════════════════════════════
# Figure Geometry
# ═══════════════════════════════════════════════════════════════════════
FIG_WIDTH   = 7
FIG_HEIGHT  = 11
N_COLS      = 3
N_ROWS      = 6


# ═══════════════════════════════════════════════════════════════════════
# Subplot Spacing
# ═══════════════════════════════════════════════════════════════════════
HSPACE        = 0.25    # Vertical space between rows
WSPACE        = 0.50    # Horizontal space between columns
LEFT_MARGIN   = 0.07
RIGHT_MARGIN  = 0.99
TOP_MARGIN    = 0.98
BOTTOM_MARGIN = 0.10   # Extended margin to accommodate multi-row legends


# ═══════════════════════════════════════════════════════════════════════
# Coordinate Range Limits (Standardized)
# ═══════════════════════════════════════════════════════════════════════

YLIM_ACC    = (-0.2,  1.7)
YLIM_DRINK  = (-0.2,  1.7)
YLIM_RATIO  = (-0.25, 0.55)
XLIM        = (-3.6,  8.6)


# ═══════════════════════════════════════════════════════════════════════
# Axis Labels
# ═══════════════════════════════════════════════════════════════════════

YLABEL_ACC   = "Estimated COEF. (IRR)"
YLABEL_DRINK = "Estimated COEF. (IRR)"
YLABEL_RATIO = "Estimated COEF. (OLS)"
XLABEL       = "Months"


# ═══════════════════════════════════════════════════════════════════════
# Reference Lines
# ═══════════════════════════════════════════════════════════════════════

# Horizontal benchmarks (1.0 for IRR models, 0.0 for OLS models)
HLINE_Y_ACC   = 1.0
HLINE_Y_DRINK = 1.0
HLINE_Y_RATIO = 0.0
HLINE_COLOR   = "black"
HLINE_ALPHA   = 0.4
HLINE_LW      = 0.8
HLINE_STYLE   = "--"

VLINE_X     = 0
VLINE_COLOR = "grey"
VLINE_ALPHA = 0.5
VLINE_LW    = 0.8
VLINE_STYLE = "--"


# ═══════════════════════════════════════════════════════════════════════
# Plot Styles (Lines and Confidence Intervals)
# ═══════════════════════════════════════════════════════════════════════

LINE_WIDTH   = 1.0
MARKER_SIZE  = 3.0
MARKER_STYLE = "o"
CI_CAPSIZE   = 2.0
CI_LW        = 0.8
CI_ALPHA     = 0.25
USE_FILL_CI  = True


# ═══════════════════════════════════════════════════════════════════════
# Panel Label Positioning (Axes coordinates)
# ═══════════════════════════════════════════════════════════════════════

PANEL_LABEL_X = -0.18
PANEL_LABEL_Y = 1.02


# ═══════════════════════════════════════════════════════════════════════
# Color Mapping by Dimension Group
# ═══════════════════════════════════════════════════════════════════════

GROUP_COLORS = {
    "car": "#E07B54",
    "crl": "#5B8DB8",
    "rak": "#6AAF6A",
    "rat": "#9B6BB5",
    "reg": "#C4A83A",
    "tim": "#C45B5B",
}


# ═══════════════════════════════════════════════════════════════════════
# Category Labels (Legend Titles)
# ═══════════════════════════════════════════════════════════════════════

CATEGORY_NAMES = {
    "car": "Vehicle Type",
    "crl": "Urban–Rural Classification",
    "rak": "Road Type",
    "rat": "Prior DUI Crash Incidence",
    "reg": "Region",
    "tim": "Time Period",
}


# ═══════════════════════════════════════════════════════════════════════
# Display Order for Legend Categories
# ═══════════════════════════════════════════════════════════════════════

LEGEND_CAT_ORDER = ["car", "crl", "rak", "rat", "reg", "tim"]


# ═══════════════════════════════════════════════════════════════════════
# Start positions for legend columns (Figure coordinates)
# ═══════════════════════════════════════════════════════════════════════

LEGEND_X_POSITIONS = None


# ═══════════════════════════════════════════════════════════════════════
# File Stems and Descriptive Labels
# ═══════════════════════════════════════════════════════════════════════

FILE_LABELS = {
    "car_C":                  "Type C",
    "car_D":                  "Type D",
    "crl_Rural areas":        "Rural",
    "crl_Suburban areas":     "Suburban",
    "crl_Urbanization areas": "Urban",
    "rak_mainRoads":          "Main roads",
    "rak_motorway":           "Motorway",
    "rak_residential":        "Residential",
    "rak_service":            "Service",
    "rat_higherRate":         "Higher",
    "rat_lowerRate":          "Lower",
    "reg_East":               "East",
    "reg_Middle":             "Middle",
    "reg_West":               "West",
    "tim_afternoon":          "Afternoon",
    "tim_midnight":           "Midnight",
    "tim_morning":            "Morning",
    "tim_night":              "Night",
}


# ═══════════════════════════════════════════════════════════════════════
# Matrix Row Layout
# ═══════════════════════════════════════════════════════════════════════

ROW_ORDER = [
    ["car_C",                 "car_D"],
    ["crl_Rural areas",       "crl_Suburban areas",       "crl_Urbanization areas"],
    ["rak_mainRoads",         "rak_motorway",              "rak_residential",  "rak_service"],
    ["rat_higherRate",        "rat_lowerRate"],
    ["reg_East",              "reg_Middle",                "reg_West"],
    ["tim_afternoon",         "tim_midnight",              "tim_morning",      "tim_night"],
]


# ═══════════════════════════════════════════════════════════════════════
# Export Configuration
# ═══════════════════════════════════════════════════════════════════════

SAVE_DPI    = 500
SAVE_FORMAT = "pdf"


# ═══════════════════════════════════════════════════════════════════════
# Implementation Logic
# ═══════════════════════════════════════════════════════════════════════

matplotlib.rcParams["font.sans-serif"] = FONT_FAMILY
matplotlib.rcParams["font.family"]     = "sans-serif"
matplotlib.rcParams["pdf.fonttype"]    = 42

# Time points mapping to relative months
X_POINTS  = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
TIME_KEYS = ["pre_3", "pre_2", "current",
             "time_1", "time_2", "time_3", "time_4",
             "time_5", "time_6", "time_7", "time_8"]
X_MAP = {k: v for k, v in zip(TIME_KEYS, [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])}
X_MAP["time_8"] = 8


def parse_txt(filepath):
    """Extracts point estimates and confidence intervals from regression text outputs."""
    result = {"acc": {}, "drink": {}, "ratio": {}}
    with open(filepath, encoding="utf-8") as f:
        lines = [l.rstrip() for l in f if l.strip()]
    i = 0
    while i < len(lines):
        parts   = lines[i].split("\t")
        varname = parts[0].strip()
        if varname in TIME_KEYS:
            try:
                acc_pt   = float(parts[2])
                drink_pt = float(parts[3])
                ratio_pt = float(parts[4])
            except Exception:
                i += 1
                continue
            if i + 1 < len(lines):
                ci_parts = lines[i + 1].split("\t")
                def parse_ci(s):
                    # Regex to extract (low, high) confidence bounds
                    m = re.findall(r"-?\d+\.\d+", s.strip().strip("()"))
                    return (float(m[0]), float(m[1])) if len(m) >= 2 else (None, None)
                acc_lo,   acc_hi   = parse_ci(ci_parts[2] if len(ci_parts) > 2 else "")
                drink_lo, drink_hi = parse_ci(ci_parts[3] if len(ci_parts) > 3 else "")
                ratio_lo, ratio_hi = parse_ci(ci_parts[4] if len(ci_parts) > 4 else "")
                i += 2
            else:
                acc_lo = acc_hi = drink_lo = drink_hi = ratio_lo = ratio_hi = None
                i += 1
            x = X_MAP[varname]
            result["acc"][x]   = (acc_pt,   acc_lo,   acc_hi)
            result["drink"][x] = (drink_pt, drink_lo, drink_hi)
            result["ratio"][x] = (ratio_pt, ratio_lo, ratio_hi)
        else:
            i += 1
    return result


def plot_series(ax, data_dict, color):
    """Renders the event study line and shaded confidence interval."""
    xs  = sorted(data_dict.keys())
    pts = [data_dict[x][0] for x in xs]
    los = [data_dict[x][1] for x in xs]
    his = [data_dict[x][2] for x in xs]
    ax.plot(xs, pts, color=color, linewidth=LINE_WIDTH,
            marker=MARKER_STYLE, markersize=MARKER_SIZE, zorder=3)
    if USE_FILL_CI and None not in los:
        ax.fill_between(xs, los, his, color=color, alpha=CI_ALPHA, zorder=1)
    elif None not in los:
        yerr_lo = [p - l for p, l in zip(pts, los)]
        yerr_hi = [h - p for p, h in zip(pts, his)]
        ax.errorbar(xs, pts, yerr=[yerr_lo, yerr_hi],
                    fmt="none", color=color, capsize=CI_CAPSIZE,
                    linewidth=CI_LW, zorder=2)


def style_ax(ax, ylim, hline_y, show_xlabel, show_ylabel, ylabel_text):
    """Applies standardized academic styling to the axes."""
    ax.set_xlim(XLIM)
    ax.set_ylim(ylim)
    ax.set_xticks(X_POINTS)
    ax.set_xticklabels([str(x) for x in X_POINTS], fontsize=FONT_SIZE_TICK)
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK)
    ax.axhline(hline_y, color=HLINE_COLOR, alpha=HLINE_ALPHA,
               linewidth=HLINE_LW, linestyle=HLINE_STYLE, zorder=0)
    ax.axvline(VLINE_X, color=VLINE_COLOR, alpha=VLINE_ALPHA,
               linewidth=VLINE_LW, linestyle=VLINE_STYLE, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_xlabel:
        ax.set_xlabel(XLABEL, fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel(ylabel_text, fontsize=FONT_SIZE_LABEL)


def draw_bottom_legend(fig):
    """Generates a two-row legend at the bottom of the figure."""
    # Geometry for legend line icons
    line_half = 0.025   # Half-width in figure coordinates
    gap       = 0.018   # Gap between icon and text

    # Row Y-positions
    y_row1 = LEGEND_TOP_MARGIN + 0.028   # Top row
    y_row2 = LEGEND_TOP_MARGIN + 0.008   # Bottom row

    # Calculate horizontal stepping for 3 categories per row
    n_per_row  = 3
    total_span = RIGHT_MARGIN - LEFT_MARGIN
    step       = total_span / n_per_row

    rows = [
        (LEGEND_CAT_ORDER[:3], y_row1),
        (LEGEND_CAT_ORDER[3:], y_row2),
    ]

    for cats, y_legend in rows:
        for j, cat in enumerate(cats):
            x_center = LEFT_MARGIN + (j + 0.5) * step
            color    = GROUP_COLORS[cat]
            label    = CATEGORY_NAMES[cat]

            # Draw representative line + marker
            line = matplotlib.lines.Line2D(
                [x_center - line_half - GAP_LEGEND, x_center + line_half - GAP_LEGEND],
                [y_legend, y_legend],
                color=color,
                linewidth=LINE_WIDTH * 1.5,
                marker=MARKER_STYLE,
                markersize=MARKER_SIZE + 1,
                transform=fig.transFigure,
                figure=fig,
            )
            fig.lines.append(line)

            # Draw text label
            fig.text(x_center + line_half + gap - GAP_LEGEND, y_legend, label,
                     fontsize=FONT_SIZE_LEGEND_CAT,
                     color="black",
                     va="center", ha="left",
                     transform=fig.transFigure)


if __name__ == "__main__":

    all_data = {}
    for stem in FILE_LABELS:
        fp = os.path.join(DATA_DIR, stem + ".txt")
        all_data[stem] = parse_txt(fp)

    fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(FIG_WIDTH, FIG_HEIGHT))
    fig.subplots_adjust(left=LEFT_MARGIN, right=RIGHT_MARGIN,
                        top=TOP_MARGIN,   bottom=BOTTOM_MARGIN,
                        hspace=HSPACE,    wspace=WSPACE)

    col_indicators = ["acc", "drink", "ratio"]
    col_ylims      = [YLIM_ACC, YLIM_DRINK, YLIM_RATIO]
    col_hlines     = [HLINE_Y_ACC, HLINE_Y_DRINK, HLINE_Y_RATIO]
    col_ylabels    = [YLABEL_ACC, YLABEL_DRINK, YLABEL_RATIO]

    panel_idx = 0

    for row_idx, file_stems in enumerate(ROW_ORDER):
        for col_idx in range(N_COLS):
            ax        = axes[row_idx, col_idx]
            indicator = col_indicators[col_idx]
            ylim      = col_ylims[col_idx]
            hline_y   = col_hlines[col_idx]
            ylabel    = col_ylabels[col_idx]

            for stem in file_stems:
                color  = GROUP_COLORS.get(stem[:3], "#333333")
                data_d = all_data[stem][indicator]
                plot_series(ax, data_d, color)

            show_xlabel = (row_idx == N_ROWS - 1)
            show_ylabel = (col_idx == 0)
            style_ax(ax, ylim, hline_y, show_xlabel, show_ylabel, ylabel)

            panel_letter = chr(ord("a") + panel_idx)
            ax.text(PANEL_LABEL_X, PANEL_LABEL_Y, panel_letter,
                    transform=ax.transAxes,
                    fontsize=FONT_SIZE_PANEL,
                    fontweight=FONT_WEIGHT_PANEL,
                    va="bottom", ha="left")
            panel_idx += 1

    draw_bottom_legend(fig)

    plt.savefig(SAVE_PATH, dpi=SAVE_DPI, bbox_inches="tight")
    plt.savefig(SAVE_PATH_PDF, dpi=SAVE_DPI, bbox_inches="tight")
    plt.show()
    print(f"Saved: {SAVE_PATH}")
