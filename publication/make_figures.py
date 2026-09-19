#!/usr/bin/env python3
"""
make_figures.py

Figures 1 to 4 for the APOL1-MYH9 manuscript.

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
# The APOL1 gene body, exactly as the analysis scripts define it.
GENE_START, GENE_END = 36_253_010, 36_267_530
CODING = {"missense_variant", "synonymous_variant", "stop_gained",
          "frameshift_variant", "inframe_deletion", "start_lost", "stop_lost",
          "splice_acceptor_variant", "splice_donor_variant"}


def figure1() -> None:
    """Why the expression channel was discarded.

    Two claims have to be visible: that predicted effect tracks position rather
    than regulation, and that the coding risk variants stay extreme after the
    correction meant to remove that. Panel A carries the first as a
    distance-matched three-group comparison, which is the honest version: an
    earlier draft compared gene-body non-coding variants against ALL coding
    variants at the locus, most of which sit in other genes, and drew the
    opposite conclusion. Panel B carries the second.
    """
    d = pd.concat([pd.read_csv(f, sep="\t")
                   for f in sorted(glob.glob(str(RESULTS / "scored" / "*.tsv.gz")))],
                  ignore_index=True)
    d = d.dropna(subset=["rna_cortex_of_kidney"]).copy()
    d["abs_eff"] = d["rna_cortex_of_kidney"].abs()
    d["in_gene"] = d["pos"].between(GENE_START, GENE_END)
    d["coding"] = d["consequence"].isin(CODING)
    d["dist"] = distance_to_apol1(d["pos"].values)

    groups = [
        ("Coding,\nin gene", d[d.in_gene & d.coding]),
        ("Non-coding,\nin gene", d[d.in_gene & ~d.coding]),
        ("Outside\nthe gene", d[~d.in_gene]),
    ]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(TWO_COL, 2.5),
                                   gridspec_kw={"width_ratios": [1, 1.15]})

    rng = np.random.default_rng(20260919)
    for i, (lab, g) in enumerate(groups):
        y = g["abs_eff"].values
        x = i + rng.uniform(-0.17, 0.17, size=len(y))
        axA.scatter(x, y, s=4, color=NEUTRAL, alpha=0.45, linewidth=0, zorder=3)
        med = float(np.median(y))
        axA.plot([i - 0.30, i + 0.30], [med, med], color=ACCENT,
                 linewidth=1.8, solid_capstyle="round", zorder=5)
        axA.annotate(f"{med:.4f}", (i + 0.32, med), fontsize=6,
                     color=ACCENT, va="center", fontweight="bold")
        # Place n in AXES coordinates: get_ylim() before the log scale is set
        # returns the default 0-1 range and pushes the label off the figure.
        axA.annotate(f"n = {len(y)}", (i, 0), xycoords=("data", "axes fraction"),
                     xytext=(0, -30), textcoords="offset points",
                     fontsize=5.8, color=INK2, ha="center")
    axA.set_yscale("log")
    axA.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=6.2)
    axA.set_xlim(-0.55, len(groups) - 0.25)
    axA.set_ylabel("Predicted absolute effect\non APOL1 expression", linespacing=1.5)
    axA.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
    axA.set_axisbelow(True)
    axA.set_title("A", loc="left", fontweight="bold", fontsize=8)

    # ---- B: the residual, and where the risk variants sit on it ----------
    lg = np.log10(np.maximum(d["dist"].values, 1.0))
    ly = np.log10(d["abs_eff"].values)
    ok = np.isfinite(lg) & np.isfinite(ly)
    b, a = np.polyfit(lg[ok], ly[ok], 1)
    resid = ly - (a + b * lg)
    rmask = d["is_risk_variant"].notna().values

    axB.hist(resid[ok & ~rmask], bins=46, color=NEUTRAL, linewidth=0, zorder=3)
    for r, name in zip(resid[rmask], d.loc[rmask, "is_risk_variant"]):
        axB.axvline(r, color=ACCENT, linewidth=1.3, zorder=5)
    pct = [100 * (resid[ok] < r).mean() for r in resid[rmask]]
    axB.annotate(f"coding risk variants,\nstill {min(pct):.0f}-{max(pct):.0f}th "
                 f"percentile\nafter correcting for distance",
                 (max(resid[rmask]), axB.get_ylim()[1] * 0.92),
                 xytext=(-6, 0), textcoords="offset points",
                 ha="right", va="top", fontsize=6.2, color=ACCENT,
                 fontweight="bold", linespacing=1.5)
    axB.set_xlabel("Residual after regressing effect on distance")
    axB.set_ylabel("Variants")
    axB.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
    axB.set_axisbelow(True)
    axB.set_title("B", loc="left", fontweight="bold", fontsize=8)

    fig.suptitle("The expression channel tracks position, and its negative control fails",
                 x=0.02, ha="left", fontweight="bold", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "figure1")


# ---------------------------------------------------------------- figure 2
def figure2() -> None:
    """The enrichment that the rest of the paper then constrains."""
    e = json.loads((RESULTS / "encode_validation.json").read_text())
    t = e["support_1"]["tests"]
    tops = [25, 50, 100, 200]
    series = {
        "Kidney ATAC": [t[f"atac_kidney_top{n}"]["odds_ratio"] for n in tops],
        "Podocyte DNase": [t[f"dnase_glomerular_visceral_epithelial_cell_top{n}"]
                           ["odds_ratio"] for n in tops],
    }
    fig, ax = plt.subplots(figsize=(ONE_COL, 2.6))
    x = np.arange(len(tops))
    for i, (lab, vals) in enumerate(series.items()):
        ax.plot(x, vals, color=ACCENT if i == 0 else NEUTRAL,
                linewidth=1.4 if i == 0 else 1.1,
                marker="D" if i == 0 else "o", markersize=4.5,
                markeredgecolor=SURFACE, markeredgewidth=0.7, zorder=4)
        ax.annotate(lab, (x[-1], vals[-1]), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=6.4,
                    color=ACCENT if i == 0 else INK2,
                    fontweight="bold" if i == 0 else "normal")
    ax.axhline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
    ax.text(-0.35, 1.06, "no enrichment", fontsize=6, color=INK2)
    ax.set_xticks(x, [f"top {n}" for n in tops])
    ax.set_xlim(-0.45, len(tops) - 0.35)
    ax.set_ylim(0, 9.4)
    ax.set_ylabel("Odds ratio for falling in measured\nkidney open chromatin",
                  linespacing=1.5)
    ax.set_xlabel("Ranking depth")
    ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("Top-ranked variants fall in measured kidney chromatin",
                 loc="left", pad=6, fontweight="bold")
    save(fig, "figure2")


# ---------------------------------------------------------------- figure 3
def figure3() -> None:
    """The tissue-label control, and the point of the paper.

    One row per tissue output used to rank. The four markers are the four
    ranking depths, so both the level and its spread are visible. Kidney rows
    are black and heavier; every other row is grey. Nothing here depends on
    hue, which is the whole idea.
    """
    e = pd.read_csv(RESULTS / "wrong_tissue_enrichment.csv")
    # Clinical readership: ENCODE biosample names are not how nephrologists
    # refer to these cells, and "glomerular visceral epithelial cell" buries
    # the one row the figure exists to highlight.
    PRETTY = {
        "glomerular_visceral_epithelial_cell": "Podocyte (kidney)",
        "kidney": "Whole kidney",
        "hepatocyte": "Hepatocyte",
        "frontal_cortex": "Frontal cortex",
        "heart_left_ventricle": "Heart, left ventricle",
        "left_lung": "Lung (left)",
    }
    e["label"] = e["ranked_by"].map(
        lambda x: PRETTY.get(x, x.replace("_", " ").capitalize()))
    order = (e.groupby("label")["odds_ratio"].median()
             .sort_values().index.tolist())
    kidney = set(e[e["is_kidney"]]["label"])

    fig, ax = plt.subplots(figsize=(TWO_COL * 0.62, 3.1))
    for i, lab in enumerate(order):
        vals = e[e["label"] == lab].sort_values("top_n")["odds_ratio"].values
        is_k = lab in kidney
        c = ACCENT if is_k else NEUTRAL
        ax.plot([vals.min(), vals.max()], [i, i], color=c,
                linewidth=2.0 if is_k else 1.2, solid_capstyle="round", zorder=3)
        ax.scatter(vals, [i] * len(vals), s=15 if is_k else 11, color=c,
                   marker="D" if is_k else "o",
                   edgecolor=SURFACE, linewidth=0.6, zorder=4)
    ax.set_yticks(range(len(order)), order, fontsize=6.4)
    for tick, lab in zip(ax.get_yticklabels(), order):
        if lab in kidney:
            tick.set_fontweight("bold")
            tick.set_color(ACCENT)
    ax.axvline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
    ax.set_xlabel("Odds ratio for falling in measured kidney open chromatin\n"
                  "(four markers = ranking depths 25, 50, 100, 200)",
                  linespacing=1.6)
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlim(0, 7.4)
    ax.set_title("Ranking by the wrong tissue finds kidney chromatin just as well",
                 loc="left", pad=6, fontweight="bold")
    save(fig, "figure3")


# ---------------------------------------------------------------- figure 4
def figure4() -> None:
    """The control locus. Same pipeline, somewhere it should not work."""
    n = json.loads((RESULTS / "negative_control_locus.json").read_text())
    e = json.loads((RESULTS / "encode_validation.json").read_text())
    t = e["support_1"]["tests"]
    tops = [25, 50, 100, 200]
    test = [t[f"dnase_glomerular_visceral_epithelial_cell_top{k}"]["odds_ratio"]
            for k in tops]
    ctrl = {r["top_n"]: r["odds_ratio"] for r in n["enrichment"]}
    ctrl = [ctrl.get(k, 0.0) for k in tops]

    fig, ax = plt.subplots(figsize=(ONE_COL, 2.6))
    x = np.arange(len(tops))
    w = 0.36
    ax.bar(x - w / 2, test, w, color=ACCENT, zorder=3,
           label="APOL1-MYH9 locus")
    ax.bar(x + w / 2, ctrl, w, color=NEUTRAL, zorder=3,
           label="beta-globin control locus")
    for xi, v in zip(x - w / 2, test):
        ax.text(xi, v + 0.18, f"{v:.1f}", ha="center", fontsize=6, color=INK2)
    for xi, v in zip(x + w / 2, ctrl):
        ax.text(xi, v + 0.18, f"{v:.1f}", ha="center", fontsize=6, color=INK2)
    ax.axhline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
    ax.set_xticks(x, [f"top {k}" for k in tops])
    ax.set_ylim(0, 7.6)
    ax.set_ylabel("Odds ratio for falling in measured\nkidney open chromatin",
                  linespacing=1.5)
    ax.set_xlabel("Ranking depth")
    ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=6.2, handlelength=1.1)
    ax.set_title("The same ranking finds nothing at a control locus",
                 loc="left", pad=6, fontweight="bold")
    save(fig, "figure4")


def main() -> None:
    print("building figures")
    figure1(); figure2(); figure3(); figure4()
    print(f"all figures in {FIGS}")


if __name__ == "__main__":
    main()
