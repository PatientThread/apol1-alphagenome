#!/usr/bin/env python3
"""
make_figures.py

Figures 1 and 2 for the APOL1-MYH9 research note.

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


def figure1() -> None:
    """Enrichment by the tissue output used to rank, across five reference sets.

    The whole argument of the paper is in the gap between the panels: the same
    predictions, the same variants, four measured reference sets, four readings.
    The control panel is the one that matters, because it shows the pipeline
    behaves at beta-globin about as it does at the disease locus.
    """
    counts = pd.read_csv(RESULTS / "reconciled_counts.csv")
    counts = counts[counts["depth"] == DEPTH]
    konly = pd.read_csv(RESULTS / "kidney_only_counts.csv")
    counts = pd.concat([counts, konly], ignore_index=True)
    ctrl = json.loads((RESULTS / "final_reconciliation.json").read_text())
    crow = [r for r in ctrl["control"]["rows"] if r["depth"] == DEPTH][0]

    panels = [("Whole kidney only", "kidney-only"),
              ("Kidney + podocyte", "pooled"),
              ("Kidney-restricted", "kidney-specific"),
              ("Shared", "shared"),
              ("Control locus\n(beta-globin)", None)]

    fig, axes = plt.subplots(1, 5, figsize=(TWO_COL, 2.05),
                             sharex=True, sharey=True)
    for ax, (title, endpoint) in zip(axes, panels):
        if endpoint is None:
            labels = ["kidney (podocyte)"]
            rows = [(crow["odds_ratio"], crow["ci_low"], crow["ci_high"])]
            ys = [TISSUES.index("kidney (podocyte)")]
        else:
            sub = counts[counts["endpoint"] == endpoint].set_index("ranked_by")
            labels, rows, ys = [], [], []
            for i, t in enumerate(TISSUES):
                if t in sub.index:
                    r = sub.loc[t]
                    labels.append(t)
                    rows.append((r["odds_ratio"], r["ci_low"], r["ci_high"]))
                    ys.append(i)
        for y, lab, (o, lo, hi) in zip(ys, labels, rows):
            kid = lab in KIDNEY
            c = ACCENT if kid else NEUTRAL
            lo = max(lo, 0.25)          # keep zero-cell bounds on a log axis
            ax.plot([lo, hi], [y, y], color=c, lw=1.4 if kid else 1.0,
                    solid_capstyle="round", zorder=2)
            ax.plot([o], [y], marker="o" if kid else "s", ms=4.2 if kid else 3.4,
                    color=c, mec=SURFACE, mew=0.7, zorder=3)
        ax.axvline(1, color=GRID, lw=0.7, zorder=1)
        ax.set_xscale("log")
        ax.set_title(title, fontsize=7, pad=4)
        ax.set_xlim(0.25, 60)
        ax.set_xticks([1, 10])
        ax.set_xticklabels(["1", "10"])
        ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
        ax.set_axisbelow(True)
    axes[0].set_yticks(range(len(TISSUES)))
    axes[0].set_yticklabels(TISSUES)
    axes[0].set_ylim(len(TISSUES) - 0.5, -0.5)
    # One shared x label rather than the same words five times.
    axes[2].set_xlabel("odds ratio, top 200 variants against the rest "
                       "(log scale)")
    axes[4].annotate("no non-renal ranking\nevaluated at this locus",
                     xy=(0.5, 0.62), xycoords="axes fraction",
                     ha="center", va="center", fontsize=6, color=INK2)
    fig.text(0.005, 0.5, "variants ranked by predicted accessibility in",
             rotation=90, va="center", ha="left", fontsize=6.5, color=INK2)
    fig.subplots_adjust(left=0.165, right=0.995, top=0.85, bottom=0.185,
                        wspace=0.10)
    save(fig, "figure1")


# ---------------------------------------------------------------- figure 2
def figure2() -> None:
    """Change in tissue contrast on partitioning. No interval excludes zero."""
    d = pd.read_csv(RESULTS / "reconciled_interaction.csv")
    d = d.sort_values("delta_observed").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(ONE_COL, 1.95))
    for i, r in d.iterrows():
        ax.plot([r.delta_ci_low, r.delta_ci_high], [i, i], color=NEUTRAL,
                lw=1.1, solid_capstyle="round", zorder=2)
        ax.plot([r.delta_observed], [i], marker="o", ms=4.2, color=ACCENT,
                mec=SURFACE, mew=0.7, zorder=3)
    ax.axvline(0, color=INK2, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d["comparator"])
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlabel("change in podocyte-minus-comparator contrast\n"
                  "on partitioning (log odds ratio)")
    ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.34, right=0.985, top=0.96, bottom=0.34)
    save(fig, "figure2")


def main() -> None:
    print("building figures")
    figure1()
    figure2()
    for stale in ("figure3", "figure4"):
        for ext in ("pdf", "png"):
            p = FIGS / f"{stale}.{ext}"
            if p.exists():
                p.unlink()
                print(f"  removed superseded {p.name}")


if __name__ == "__main__":
    main()
