"""
SI-FIG1: Heterogeneity Analysis of Drunk Driving Trends
Author: Hui Liu
Date: 2026-04-15
Description: This script visualizes the temporal trends of drunk driving cases,
             crashes, and incidence rates. It compares historical data (2016–2019)
             with the COVID-19 pandemic period (2019.10–2020.09) across three
             distinct phases: Scenario I (Pre-lockdown), Scenario II (Lockdown), and
             Scenario III (Recovery).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker
from config import __version, __date

# ═══════════════════════════════════════════════════════════════════════
# File Paths
# ═══════════════════════════════════════════════════════════════════════

SAVE_FP = f"../output/{__date}-{__version}"
LOAD_FP = f"../data/fig1/pop_1year_all_fig0_1211.txt"


# ═══════════════════════════════════════════════════════════════════════
# Font Settings
# ═══════════════════════════════════════════════════════════════════════

FONT_FAMILY       = "Arial"
FONT_SIZE_TICK_X  = 8       # X-axis tick labels
FONT_SIZE_TICK_Y  = 8       # Y-axis tick labels
FONT_SIZE_YLABEL  = 8       # Y-axis title
FONT_SIZE_LEGEND  = 7       # Legend font size
FONT_SIZE_PANEL   = 10      # Subplot labels (a/b/c)
FONT_WEIGHT_PANEL = 'bold'

# ═══════════════════════════════════════════════════════════════════════
# Figure Geometry and Layout
# ═══════════════════════════════════════════════════════════════════════

FIG_SIZE    = (7, 2.8)   # (Width, Height) in inches
AX_Y_POS    = 0.20       # Bottom margin of subplots (figure coordinates)
AX_WIDTH    = 0.265      # Width of each subplot
AX_HEIGHT   = 0.70       # Height of each subplot

# Horizontal starting positions for the three subplots
AX_LEFT_POSITIONS = {
    'a': 0.06,
    'b': 0.39,
    'c': 0.72,
}

# Vertical Reference Line (Intervention Point)
VLINE_X         = 3           # Vertical Reference Line (Intervention Point)
VLINE_COLOR     = "orange"
VLINE_ALPHA     = 0.6
VLINE_LINESTYLE = "--"


# ═══════════════════════════════════════════════════════════════════════
# Scatter Plots: Historical Baseline (2016–2019)
# ═══════════════════════════════════════════════════════════════════════
HIST_COLORS = ['purple', 'y', 'g']
HIST_LABELS = [
    "2016.10–2017.09",
    "2017.10–2018.09",
    "2018.10–2019.09",
]
HIST_MARKER = "v"
HIST_SIZE   = 7       # Marker size
HIST_ALPHA  = 0.5


# ═══════════════════════════════════════════════════════════════════════
# Line Plot: Historical Mean (Gray Dashed Line)
# ═══════════════════════════════════════════════════════════════════════

MEAN_COLOR      = "grey"
MEAN_ALPHA      = 0.4
MEAN_MARKER     = "^"
MEAN_MARKERSIZE = 3
MEAN_LINEWIDTH  = 1
MEAN_LINESTYLE  = '--'
MEAN_LABEL      = "2016.10–2019.09"


# ═══════════════════════════════════════════════════════════════════════
# Line Plots: Pandemic Year (Segmented by Phase)
# ═══════════════════════════════════════════════════════════════════════

# Format: (Start Index, End Index inclusive, Color, Label)
COVID_SEGMENTS = [
    (0,  3,  "grey",  "2019.10–2020.01"),   # Pre-pandemic
    (3,  6,  "red",   "2020.01–2020.04"),   # Lockdown period
    (6,  12, "blue",  "2020.04–2020.09"),   # Recovery period
]
COVID_MARKER     = "o"
COVID_MARKERSIZE = 3
COVID_LINEWIDTH  = 1
COVID_ALPHA      = 0.7
COVID_ZORDER     = 3


# ═══════════════════════════════════════════════════════════════════════
# Legend Configuration (Justified Alignment)
# ═══════════════════════════════════════════════════════════════════════

LEGEND_PANEL   = "b"
LEGEND_FRAMEON = False

# Legend Row 1: Pandemic phases
LEGEND_ROW1_LABELS = [
    "2019.10–2020.01",
    "2020.01–2020.04",
    "2020.04–2020.09",
]
# Legend Row 2: Historical years and baseline mean
LEGEND_ROW2_LABELS = [
    "2016.10–2017.09",
    "2017.10–2018.09",
    "2018.10–2019.09",
    "2016.10–2019.09",
]
# Boundaries for justified legend rows (figure coordinates)
LEGEND_X_LEFT  = 0.25
LEGEND_X_RIGHT = 0.75
LEGEND_Y_ROW1  = 0.08
LEGEND_Y_ROW2  = 0.03


# ═══════════════════════════════════════════════════════════════════════
# Panel Label Positions
# ═══════════════════════════════════════════════════════════════════════

PANEL_TEXT_POS = {
    'a': {'text_x': -2.5, 'text_y': 43500},
    'b': {'text_x': -2.5, 'text_y': 13000},
    'c': {'text_x': -2,   'text_y':  1.09},
}


# ═══════════════════════════════════════════════════════════════════════
# Y-Axis Limits
# ═══════════════════════════════════════════════════════════════════════

YLIMS = {
    'a': [0, 40000],
    'b': [0, 12000],
    'c': [0, 1],
}


# ═══════════════════════════════════════════════════════════════════════
# Y-Axis Labels
# ═══════════════════════════════════════════════════════════════════════
YLABELS = {
    'a': 'Num. of cases',
    'b': 'Num. of crashes',
    'c': 'Drunk Driving Crash Incidence',
}


# ═══════════════════════════════════════════════════════════════════════
# Data Mapping for Subplots
# ═══════════════════════════════════════════════════════════════════════

# Options: 'drunks' / 'accs' / 'incidences'
PANEL_DATA_MAP = {
    'a': 'drunks',
    'b': 'accs',
    'c': 'incidences',
}


# ═══════════════════════════════════════════════════════════════════════
# Export Parameters
# ═══════════════════════════════════════════════════════════════════════

SAVE_DPI    = 500

matplotlib.rcParams['font.sans-serif'] = FONT_FAMILY
matplotlib.rcParams['font.family']     = "sans-serif"

x_label_shunxu = [
    ("10", 'Oct.',  '-3'),
    ('11', 'Nov.',  '-2'),
    ('12', 'Dec.',  '-1'),
    ('01', 'Jan.',   '0'),
    ('02', 'Feb.',   '1'),
    ('03', 'Mar.',   '2'),
    ('04', 'Apr.',   '3'),
    ('05', 'May.',   '4'),
    ('06', 'Jun.',   '5'),
    ('07', 'Jul.',   '6'),
    ('08', 'Aug.',   '7'),
    ('09', "Sep.",   '8'),
]
years = ['2017', '2018', '2019', '2020']


def get_month_data():
    """Parses raw data and calculates descriptive statistics (mean, max, min)."""
    data = {}
    with open(LOAD_FP, 'r', encoding='utf-8') as f:
        for line in f.read().splitlines():
            sp        = line.split(',')
            year_     = sp[0].split('-')[0]
            month_    = sp[0].split('-')[1]
            data[(year_, month_)] = [int(sp[1]), int(sp[2]), float(sp[3])]

    drunks, accs, incidences = [], [], []
    for year in years:
        d, a, i = [], [], []
        for month in x_label_shunxu:
            key = (year, month[0]) if int(month[0]) <= 9 \
                  else (str(int(year) - 1), month[0])
            d.append(data[key][0])
            a.append(data[key][1])
            i.append(data[key][2])
        drunks.append(d)
        accs.append(a)
        incidences.append(i)

    for indicator in [drunks, accs, incidences]:
        arr = np.asarray(indicator)[:-1, :]
        indicator.append(list(np.average(arr, axis=0)))  # mean
        indicator.append(list(np.max(arr,     axis=0)))  # max
        indicator.append(list(np.min(arr,     axis=0)))  # min

    return drunks, accs, incidences


def plot_single(ax, title_name, data, ylim, text_x, text_y, y_label):
    """Renders a single subplot including historical dots, mean line, and pandemic segments."""
    covid_index = len(data) - 4

    ax.axvline(VLINE_X, linestyle=VLINE_LINESTYLE,
               color=VLINE_COLOR, alpha=VLINE_ALPHA, zorder=0)

    x = np.arange(len(x_label_shunxu))

    # Explicitly track handles for the justified legend to prevent duplicates
    legend_handles = {}   # label -> handle

    for i, d in enumerate(data):

        # Plot historical baseline years (Scatter)
        if i < 3:
            sc = ax.scatter(x, d,
                            color=HIST_COLORS[i], alpha=HIST_ALPHA,
                            s=HIST_SIZE, marker=HIST_MARKER)
            legend_handles[HIST_LABELS[i]] = sc
            continue

        # Plot multi-year mean baseline (Line)
        if i == len(data) - 3:
            ln, = ax.plot(x, d,
                          color=MEAN_COLOR, alpha=MEAN_ALPHA,
                          marker=MEAN_MARKER, markersize=MEAN_MARKERSIZE,
                          linewidth=MEAN_LINEWIDTH, linestyle=MEAN_LINESTYLE,
                          zorder=0)
            legend_handles[MEAN_LABEL] = ln
            continue

        # Plot pandemic year in segmented phases (Line)
        if i == covid_index:
            for (x0, x1, color, label) in COVID_SEGMENTS:
                ln, = ax.plot(x[x0:x1+1], d[x0:x1+1],
                              color=color,
                              marker=COVID_MARKER,
                              markersize=COVID_MARKERSIZE,
                              linewidth=COVID_LINEWIDTH,
                              zorder=COVID_ZORDER,
                              alpha=COVID_ALPHA)
                legend_handles[label] = ln
            continue

    # Implement Justified Multi-row Legend
    if title_name == LEGEND_PANEL:
        fig = ax.get_figure()

        def place_row(labels, y_fig):
            """Distributes legend items evenly across a specified width at height y_fig."""
            n = len(labels)
            if n == 1:
                xs = [(LEGEND_X_LEFT + LEGEND_X_RIGHT) / 2]
            else:
                step = (LEGEND_X_RIGHT - LEGEND_X_LEFT) / (n - 1)
                xs = [LEGEND_X_LEFT + i * step for i in range(n)]
            for x_fig, lb in zip(xs, labels):
                if lb not in legend_handles:
                    continue
                fig.legend(
                    [legend_handles[lb]], [lb],
                    loc='center',
                    bbox_to_anchor=(x_fig, y_fig),
                    bbox_transform=fig.transFigure,
                    fontsize=FONT_SIZE_LEGEND,
                    frameon=LEGEND_FRAMEON,
                    handlelength=1.5,
                    handletextpad=0.4,
                )

        place_row(LEGEND_ROW1_LABELS, LEGEND_Y_ROW1)
        place_row(LEGEND_ROW2_LABELS, LEGEND_Y_ROW2)

    # Axis Formatting
    ax.set_ylim(ylim)
    ax.set_xticks(x)
    ax.set_xticklabels([it[2] for it in x_label_shunxu])   # Use relative month indices
    ax.set_ylabel(y_label, fontsize=FONT_SIZE_YLABEL)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.text(x=text_x, y=text_y, s=title_name,
            fontsize=FONT_SIZE_PANEL, fontweight=FONT_WEIGHT_PANEL)
    ax.tick_params(axis='y', labelsize=FONT_SIZE_TICK_Y)
    ax.tick_params(axis='x', labelsize=FONT_SIZE_TICK_X)

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.offsetText.set_fontsize(FONT_SIZE_TICK_Y)

    return ax


if __name__ == "__main__":

    drunks, accs, incidences = get_month_data()
    data_pool = {'drunks': drunks, 'accs': accs, 'incidences': incidences}

    fig = plt.figure(figsize=FIG_SIZE)
    axes = {
        key: fig.add_axes([AX_LEFT_POSITIONS[key], AX_Y_POS, AX_WIDTH, AX_HEIGHT])
        for key in ('a', 'b', 'c')
    }

    for key in ('a', 'b', 'c'):
        plot_single(
            ax         = axes[key],
            title_name = key,
            data       = data_pool[PANEL_DATA_MAP[key]],
            ylim       = YLIMS[key],
            text_x     = PANEL_TEXT_POS[key]['text_x'],
            text_y     = PANEL_TEXT_POS[key]['text_y'],
            y_label    = YLABELS[key],
        )

    save_name = f'SI_fig1.png'
    save_name_pdf = f'SI_fig1.pdf'
    plt.savefig(os.path.join(SAVE_FP, save_name), dpi=SAVE_DPI)
    plt.savefig(os.path.join(SAVE_FP, save_name_pdf), dpi=SAVE_DPI)
    plt.show()
