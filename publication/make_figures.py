#!/usr/bin/env python3
"""
make_figures.py

Figures 1 and 2 for the APOL1-MYH9 research note.

Figure 1 is a six-panel A-F layout matching the six rows of Table 1 exactly,
one panel per locus and reference combination, all at ranking depth 200, on one
common log scale. Only estimates actually computed for a panel are drawn: the
comparison locus was evaluated for the podocyte ranking only, so panels E and F
carry a single estimate each and say so.

BMC Research Notes allows three display items in total. One is Table 1, which
documents every measured reference set, so only two figures are drawn. The
locus map, the raw count tables and the linkage and candidate figures drawn for
the earlier drafts are deposited with the code rather than printed.

GREYSCALE BY DESIGN, as on the sister paper. Identity is never carried by hue:
it is carried by marker shape, line weight and label weight, so these print
correctly in monochrome and cost nothing in colour charges at any journal.

Vector PDF is the submission format; PNG is for reading on screen. Figure text
is capped at 7pt and sized to 88 mm single column or 180 mm double column, so
production does not rescale and push labels below legibility.

Author: Christopher Lawrence
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.use("Agg")

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGS = ROOT / "publication" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
ONE_COL, TWO_COL = 88 * MM, 180 * MM

ACCENT = "#111111"     # the highlighted series: near-black
NEUTRAL = "#8a8a84"    # everything else: mid grey
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d8d8d4"

# APOL1 gene body, GRCh38
APOL1_START, APOL1_END = 36_253_071, 36_267_530

mpl.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK,
    "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def save(fig, name: str) -> None:
    fig.savefig(FIGS / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.pdf (submission) and {name}.png (screen)")


def distance_to_apol1(pos: np.ndarray) -> np.ndarray:
    """Zero inside the gene body, otherwise distance to the nearer edge."""
    d = np.zeros_like(pos, dtype=float)
    left = pos < APOL1_START
    right = pos > APOL1_END
    d[left] = APOL1_START - pos[left]
    d[right] = pos[right] - APOL1_END
    return d


# ---------------------------------------------------------------- figure 1
TISSUES = ["kidney (podocyte)", "kidney (whole)", "lung", "hepatocyte",
           "stomach", "brain", "liver"]
KIDNEY = {"kidney (podocyte)", "kidney (whole)"}
DEPTH = 200

# (panel letter, title, source, endpoint key)
PANELS = [("A", "APOL1-MYH9\nwhole kidney",            "counts", "kidney-only"),
          ("B", "APOL1-MYH9\npodocyte-inclusive union", "counts", "pooled"),
          ("C", "union\nkidney-restricted",             "counts", "kidney-specific"),
          ("D", "union\nshared",                        "counts", "shared"),
          ("E", "beta-globin\nwhole kidney",            "control", "whole"),
          ("F", "beta-globin\npodocyte-inclusive union", "control", "union")]


def figure1() -> None:
    """Six reference sets, one ranking depth, one scale."""
    counts = pd.read_csv(RESULTS / "reconciled_counts.csv")
    counts = counts[counts["depth"] == DEPTH]
    konly = pd.read_csv(RESULTS / "kidney_only_counts.csv")
    counts = pd.concat([counts, konly], ignore_index=True)

    t1 = pd.read_csv(RESULTS / "table_one.csv")
    ctrl = {"whole": t1[(t1.locus == "beta-globin") &
                        (t1.reference_set == "whole kidney only")].iloc[0],
            "union": t1[(t1.locus == "beta-globin") &
                        (t1.reference_set == "kidney + podocyte")].iloc[0]}

    fig, axes = plt.subplots(2, 3, figsize=(TWO_COL, 3.5),
                             sharex=True, sharey=True)
    for ax, (letter, title, src, key) in zip(axes.ravel(), PANELS):
        if src == "control":
            r = ctrl[key]
            rows = [(TISSUES.index("kidney (podocyte)"), "kidney (podocyte)",
                     r.odds_ratio, r.ci_low, r.ci_high)]
            ax.annotate("non-renal rankings\nnot evaluated here",
                        xy=(0.5, 0.60), xycoords="axes fraction", ha="center",
                        va="center", fontsize=6, color=INK2)
        else:
            sub = counts[counts["endpoint"] == key].set_index("ranked_by")
            rows = [(i, t, sub.loc[t]["odds_ratio"], sub.loc[t]["ci_low"],
                     sub.loc[t]["ci_high"])
                    for i, t in enumerate(TISSUES) if t in sub.index]
        for y, lab, o, lo, hi in rows:
            kid = lab in KIDNEY
            c = ACCENT if kid else NEUTRAL
            lo = max(lo, 0.12)          # keep wide lower bounds on the axis
            ax.plot([lo, hi], [y, y], color=c, lw=1.4 if kid else 1.0,
                    solid_capstyle="round", zorder=2)
            ax.plot([o], [y], marker="o" if kid else "s",
                    ms=4.2 if kid else 3.4, color=c, mec=SURFACE, mew=0.7,
                    zorder=3)
        ax.axvline(1, color=INK2, lw=0.7, zorder=1)
        ax.set_xscale("log")
        ax.set_xlim(0.12, 70)
        ax.set_xticks([0.2, 1, 5, 20])
        ax.set_xticklabels(["0.2", "1", "5", "20"])
        ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(title, fontsize=6.8, pad=3)
        ax.annotate(letter, xy=(0.012, 1.10), xycoords="axes fraction",
                    fontsize=8, fontweight="bold", va="top", ha="left")
    for ax in axes[:, 0]:
        ax.set_yticks(range(len(TISSUES)))
        ax.set_yticklabels(TISSUES)
        ax.set_ylim(len(TISSUES) - 0.5, -0.5)
    for ax in axes[1, :]:
        ax.set_xlabel("odds ratio (log scale)")
    fig.text(0.004, 0.55, "variants ranked by predicted accessibility in",
             rotation=90, va="center", ha="left", fontsize=6.5, color=INK2)
    fig.subplots_adjust(left=0.175, right=0.995, top=0.90, bottom=0.115,
                        wspace=0.10, hspace=0.42)
    save(fig, "figure1")


# ---------------------------------------------------------------- figure 2
def figure2() -> None:
    """Change in tissue contrast on partitioning, with numerical limits."""
    d = pd.read_csv(RESULTS / "shared_draw_interaction.csv")
    d = d.sort_values("delta_observed").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(TWO_COL * 0.62, 2.55))
    for i, r in d.iterrows():
        ax.plot([r.delta_ci_low, r.delta_ci_high], [i, i], color=NEUTRAL,
                lw=1.1, solid_capstyle="round", zorder=2)
        ax.plot([r.delta_observed], [i], marker="o", ms=4.2, color=ACCENT,
                mec=SURFACE, mew=0.7, zorder=3)
        ax.annotate(f"{r.delta_observed:+.2f} "
                    f"[{r.delta_ci_low:+.2f}, {r.delta_ci_high:+.2f}]",
                    xy=(1.30, i), xycoords=("axes fraction", "data"),
                    fontsize=6, va="center", ha="right", color=INK2)
    ax.axvline(0, color=INK2, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["comparator"])
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlabel("podocyte-minus-comparator contrast:\n"
                  "kidney-restricted minus union (log odds ratio)\n"
                  "positive = podocyte gains where kidney peaks are restricted\n"
                  "points and 95% percentile intervals, Haldane-corrected\n"
                  "log cross-product ratio, 500 shared draws")
    ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.235, right=0.695, top=0.97, bottom=0.44)
    save(fig, "figure2")


def main() -> None:
    print("building figures")
    figure1()
    figure2()


if __name__ == "__main__":
    main()
