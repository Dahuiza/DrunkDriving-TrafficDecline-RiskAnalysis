"""
FIG4: Heterogeneity Analysis
Author: Hui Liu
Date: 2026-04-15
Description: This script parses multiple regression result files to perform a
             subgroup (heterogeneity) analysis. It visualizes variations in
             policy impact across regions, road types, vehicle types, and time
             periods for both Short-term (Scenario II) and Medium-term (Scenario III)
             windows.
"""

import os
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from config import __date, __version

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


# ═══════════════════════════════════════════════════════════════════════
# Path & Output Configuration
# ═══════════════════════════════════════════════════════════════════════

DATA_DIR = f"../data/reg_results/{__date}-{__version}/01_main_hetero"
OUTPUT   = f"../output/{__date}-{__version}/fig4.png"
OUTPDF   = f"../output/{__date}-{__version}/fig4.pdf"
DPI      = 500


# ═══════════════════════════════════════════════════════════════════════
# Category and Mappings
# ═══════════════════════════════════════════════════════════════════════

label_name = {
    "All Pool Samples":   "Baseline",
    "Urbanization areas": "Urbanization areas",
    "Suburban areas":     "Suburban areas",
    "Rural areas":        "Rural areas",
    "East":  "East",
    "Middle":"Middle",
    "West":  "West",
    "mainRoads":   "Main Roads",
    "motorway":    "Motorway",
    "residential": "Residential",
    "service":     "Service",
    "B": "Trunk",
    "C":  "Cars",
    "D":  "Motorcycles",
    "afternoon": "Afternoon",
    "midnight":  "Midnight",
    "morning":   "Morning",
    "night":     "Night",
    "higherRate":"High Risk",
    "lowerRate": "Low Risk",
}

category_map = {
    "pop": "Baseline",
    "crl": "Urban–Rural\nClassification",
    "reg": "Region",
    "rak": "Road Type",
    "car": "Vehicle Type",
    "tim": "Time Period",
    "rat": "Prior DUI\nCrash Incidence",
}

group_order = [
    "Baseline",
    "Urban–Rural\nClassification",
    "Region",
    "Road Type",
    "Vehicle Type",
    "Time Period",
    "Prior DUI\nCrash Incidence",
]


# ═══════════════════════════════════════════════════════════════════════
# File parsing
# ═══════════════════════════════════════════════════════════════════════

def parse_file(filepath):
    with open(filepath, encoding="utf-8") as f:
        lines = f.readlines()
    result = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        for var in ("saps_short", "saps_medium"):
            if line.startswith(var):
                parts = re.split(r'\t', lines[i].rstrip('\n'))
                try:
                    acc   = float(parts[2].strip().replace(',', ''))
                    drink = float(parts[3].strip().replace(',', ''))
                    ratio = float(parts[4].strip().replace(',', ''))
                    coef_vals = {"acc": acc, "drink": drink, "ratio": ratio}
                except:
                    coef_vals = {"acc": np.nan, "drink": np.nan, "ratio": np.nan}
                j = i + 1
                while j < len(lines) and lines[j].strip() == '':
                    j += 1
                ci_line  = lines[j].rstrip('\n') if j < len(lines) else ''
                ci_parts = re.split(r'\t', ci_line)
                cis = {}
                for col_idx, col_name in zip([2, 3, 4], ["acc", "drink", "ratio"]):
                    try:
                        ci_str = ci_parts[col_idx].strip().strip('()')
                        lo_hi  = re.split(r'\s+-\s+', ci_str, maxsplit=1)
                        lo, hi = float(lo_hi[0].strip()), float(lo_hi[1].strip())
                        cis[col_name] = (lo, hi)
                    except:
                        cis[col_name] = (np.nan, np.nan)
                result[var] = {"coef": coef_vals, "ci": cis}
        i += 1
    return result


def get_label_and_group(filename):
    """Maps filenames to human-readable categories and labels."""
    name   = os.path.splitext(os.path.basename(filename))[0]
    parts  = name.split('_', 1)
    prefix = parts[0].lower()
    suffix = parts[1] if len(parts) > 1 else ""
    suffix_readable = suffix.replace('_', ' ')
    group  = category_map.get(prefix, prefix)
    label  = label_name.get(suffix_readable, suffix_readable)
    return group, label


def load_all(data_dir):
    """Iterates through data directory and builds the data structure."""
    files = sorted([f for f in os.listdir(data_dir) if f.endswith('.txt')])
    data  = {}
    for fname in files:
        print(fname)
        if "car_B" not in fname:
            fpath = os.path.join(data_dir, fname)
            try:
                parsed = parse_file(fpath)
            except Exception as e:
                print(f"  Skipping {fname}: {e}")
                continue
            if not parsed:
                continue
            group, label = get_label_and_group(fname)
            if group not in data:
                data[group] = []
            data[group].append((label, parsed.get("saps_short"), parsed.get("saps_medium")))
    return data


# ═══════════════════════════════════════════════════════════════════════
# Layout builder
# ═══════════════════════════════════════════════════════════════════════

ROW_H = 1.0   # Height per row
GAP   = 0.5   # Extra gap between groups

def build_layout(data, group_order):
    """Calculates Y-coordinates for every row and grouping header."""
    rows      = []
    sep_lines = []
    y         = 0.0

    for g_idx, grp in enumerate(group_order):
        if grp not in data:
            continue
        is_overall = (grp == "Overall")

        if g_idx > 0 and rows:
            y -= GAP

        if g_idx > 0:
            sep_lines.append(y)   # separator exactly at this header y

        if is_overall:
            # single row: bold label + data plotted here (no DODGE offset needed
            # for the label tick — data still uses DODGE visually)
            rows.append({"type": "overall", "y": y, "label": grp, "group": grp,
                         "short": data[grp][0][1], "medium": data[grp][0][2]})
            y -= ROW_H
        else:
            rows.append({"type": "header", "y": y, "label": grp, "group": grp})
            y -= ROW_H
            for (label, short_d, med_d) in data[grp]:
                rows.append({"type": "data", "y": y, "label": label,
                             "group": grp, "short": short_d, "medium": med_d})
                y -= ROW_H

    MARGIN = ROW_H * 0.55
    y_min  = y + ROW_H - MARGIN
    y_max  = 0.0 + MARGIN

    return rows, sep_lines, y_min, y_max


# ═══════════════════════════════════════════════════════════════════════
# Main plot
# ═══════════════════════════════════════════════════════════════════════

def plot(data, group_order, output, outpdf, dpi):
    DODGE = 0.18

    rows, sep_lines, y_min, y_max = build_layout(data, group_order)
    total_height = abs(y_max - y_min)

    fig_width  = 7
    fig_height = max(7, total_height * 0.28)

    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_height),
                             sharey=True,
                             gridspec_kw={"wspace": 0.06})

    col_keys     = ["drink",  "acc",   "ratio"]
    panel_labels = ["a",      "b",     "c"]
    ref_xvals    = [1.0,      1.0,     0.0]

    color_short  = "#E87461"
    color_medium = "#4F9E8F"
    ms  = 4.4
    lw  = 1.4
    cap = 2.5

    label_lookup = {row["label"]: (row["short"], row["medium"])
                    for row in rows if row["type"] in ("data", "overall")}

    for ax_idx, (ax, col, panel_lbl, ref_x) in enumerate(
            zip(axes, col_keys, panel_labels, ref_xvals)):

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis='y', left=False, labelsize=0)
        ax.tick_params(axis='x', labelsize=9)
        ax.grid(False)

        if col != 'ratio':
            ax.set_xlabel("Estimated Coefficient (IRR)", fontsize=10)
        else:
            ax.set_xlabel("Estimated Coefficient", fontsize=10)

        ax.axvline(ref_x, color="gray", lw=0.9, ls="--", alpha=0.8, zorder=1)

        ax.text(0.01, 1.01, panel_lbl,
                transform=ax.transAxes,
                fontsize=12, fontweight="bold",
                va="bottom", ha="left", color="#111111")

        # error bars
        for row in rows:
            if row["type"] not in ("data", "overall"):
                continue
            y_c  = row["y"]
            pair = label_lookup.get(row["label"])
            if pair is None:
                continue
            short_d, med_d = pair

            for var_d, y_pos, color, marker in [
                    (short_d, y_c + DODGE, color_short,  "s"),
                    (med_d,   y_c - DODGE, color_medium, "D")]:
                if var_d is None:
                    continue
                coef = var_d["coef"].get(col, np.nan)
                lo, hi = var_d["ci"].get(col, (np.nan, np.nan))
                if np.isnan(coef) or np.isnan(lo):
                    continue
                ax.errorbar(coef, y_pos,
                            xerr=[[coef - lo], [hi - coef]],
                            fmt=marker,
                            color=color, ecolor=color,
                            markersize=ms, capsize=cap,
                            lw=lw, elinewidth=lw,
                            markerfacecolor='white',  # 设置标记内部为空白
                            markeredgecolor=color,  # 设置标记边缘颜色（可选）
                            zorder=3)

        # separator lines at header y of non-first groups
        for sep_y in sep_lines:
            ax.axhline(sep_y, color="#bbbbbb", lw=0.7, ls="-", zorder=0)

        ax.set_ylim(y_min, y_max)

        # ── x ticks: ~4 intervals = 5 tick marks, clean round numbers ───────
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        # draw light vertical grid at those ticks
        for xv in ax.get_xticks():
            ax.axvline(xv, color="gray", lw=0.4, alpha=0.35, zorder=0)

    # Y-axis tick labels on axes[0]
    header_rows    = [r for r in rows if r["type"] == "header"]
    all_label_rows = [r for r in rows if r["type"] in ("header", "data", "overall")]
    header_ys      = {r["y"] for r in header_rows}
    overall_ys     = {r["y"] for r in rows if r["type"] == "overall"}

    axes[0].set_yticks([r["y"] for r in all_label_rows])
    axes[0].set_yticklabels([r["label"] for r in all_label_rows],
                             fontsize=9, ha="right")
    axes[0].tick_params(axis='y', left=False)

    for txt, r in zip(axes[0].get_yticklabels(), all_label_rows):
        if r["y"] in header_ys or r["y"] in overall_ys:
            txt.set_fontweight("bold")
            txt.set_fontsize(10)

    # Legend below the plot, with breathing room from xlabel
    legend_elements = [
        Line2D([0], [0], marker="s", color=color_short,  lw=lw, markerfacecolor='white',
               markersize=ms, label="Scenario II"),
        Line2D([0], [0], marker="D", color=color_medium, lw=lw, markerfacecolor='white',
               markersize=ms, label="Scenario III"),
    ]
    fig.legend(handles=legend_elements, loc="lower center",
               ncol=2, fontsize=10, frameon=False,
               bbox_to_anchor=(0.5, 0.0))

    # bottom=0.10 gives room for xlabel + legend without crowding
    fig.subplots_adjust(
        left   = 0.17,
        right  = 0.99,
        top    = 0.96,
        bottom = 0.10,
        wspace = 0.06,
    )

    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    fig.savefig(outpdf, dpi=dpi, bbox_inches="tight")
    print(f"Saved → {output}")
    plt.show()
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Loading data...")
    data = load_all(DATA_DIR)
    print(f"Groups found: {list(data.keys())}")
    plot(data, group_order, OUTPUT, OUTPDF, DPI)
