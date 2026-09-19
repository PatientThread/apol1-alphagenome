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
def figure1() -> None:
    """The finding: partitioning reverses the pooled result.

    Three panels, one per partition, same seven rankings in the same order in
    each. Reading across a row shows what changes when the regions are separated
    by tissue restriction. Kidney rankings are black; identity is carried by
    weight and marker, never by hue.
    """
    e = pd.read_csv(RESULTS / "specific_vs_shared.csv")
    e = e.dropna(subset=["odds_ratio"])
    parts = ["kidney-specific", "shared", "other-tissue-specific"]
    titles = ["A  Kidney-specific regions", "B  Shared regions",
              "C  Other-tissue-specific regions"]
    order = (e[e.partition == "kidney-specific"]
             .sort_values("odds_ratio")["ranked_by"].tolist())

    fig, axes = plt.subplots(1, 3, figsize=(TWO_COL, 2.7), sharey=True)
    for ax, part, title in zip(axes, parts, titles):
        sub = e[e.partition == part].set_index("ranked_by")
        for i, name in enumerate(order):
            if name not in sub.index:
                continue
            r = sub.loc[name]
            is_k = name.startswith("kidney")
            c = ACCENT if is_k else NEUTRAL
            hi = r.ci_high if pd.notna(r.ci_high) else 20
            ax.plot([r.ci_low, hi], [i, i], color=c,
                    linewidth=1.9 if is_k else 1.1,
                    solid_capstyle="round", zorder=3)
            ax.scatter([r.odds_ratio], [i], s=20 if is_k else 13, color=c,
                       marker="D" if is_k else "o",
                       edgecolor=SURFACE, linewidth=0.6, zorder=4)
        ax.axvline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
        ax.set_xscale("log")
        ax.set_xlim(0.3, 16)
        ax.set_xticks([0.5, 1, 2, 5, 10], ["0.5", "1", "2", "5", "10"])
        ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.set_title(title, loc="left", fontsize=6.8, fontweight="bold", pad=5)
        ax.set_xlabel("Odds ratio (95% CI)", fontsize=6.6)
    axes[0].set_yticks(range(len(order)), order, fontsize=6.4)
    for tick, name in zip(axes[0].get_yticklabels(), order):
        if name.startswith("kidney"):
            tick.set_fontweight("bold"); tick.set_color(ACCENT)
    fig.tight_layout()
    save(fig, "figure1")


# ---------------------------------------------------------------- figure 2
def figure2() -> None:
    """Paired differences, peaks resampled, in the kidney-specific partition."""
    j = json.loads((RESULTS / "specific_vs_shared.json").read_text())
    pr = j.get("paired_kidney_specific", {})
    if not pr:
        print("  figure2 skipped: no paired results")
        return
    items = sorted(pr.items(), key=lambda kv: kv[1]["mean_log_or_diff"])
    fig, ax = plt.subplots(figsize=(ONE_COL, 2.4))
    for i, (name, v) in enumerate(items):
        sig = v["excludes_zero"]
        c = ACCENT if sig else NEUTRAL
        ax.plot([v["lo"], v["hi"]], [i, i], color=c, linewidth=1.6,
                solid_capstyle="round", zorder=3)
        ax.scatter([v["mean_log_or_diff"]], [i], s=17, color=c, marker="D",
                   edgecolor=SURFACE, linewidth=0.6, zorder=4)
    ax.axvline(0, color=INK2, linewidth=0.8, zorder=2)
    ax.set_yticks(range(len(items)), [k for k, _ in items], fontsize=6.4)
    ax.set_xlabel("Difference in log odds ratio\n"
                  "podocyte ranking minus each alternative", linespacing=1.6,
                  fontsize=6.8)
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.text(0.02, 1.02, "favours the alternative  |  favours podocyte",
            transform=ax.transAxes, fontsize=6, color=INK2)
    ax.set_title("Kidney-specific regions only", loc="left",
                 pad=14, fontweight="bold", fontsize=7.4)
    save(fig, "figure2")


# ---------------------------------------------------------------- figure 3
# The APOL1 gene body, exactly as the analysis scripts define it.
GENE_START, GENE_END = 36_253_010, 36_267_530
CODING = {"missense_variant", "synonymous_variant", "stop_gained",
          "frameshift_variant", "inframe_deletion", "start_lost", "stop_lost",
          "splice_acceptor_variant", "splice_donor_variant"}


def figure3() -> None:
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
    save(fig, "figure3")


# ---------------------------------------------------------------- figure 4
def figure4() -> None:
    """The control locus, with intervals. No exclusion is claimed."""
    n = json.loads((RESULTS / "negative_control_locus.json").read_text())
    e = json.loads((RESULTS / "encode_validation.json").read_text())
    t = e["support_1"]["tests"]
    tops = [25, 50, 100, 200]
    test = [t[f"dnase_glomerular_visceral_epithelial_cell_top{k}"]["odds_ratio"]
            for k in tops]
    ctrl = {r["top_n"]: r for r in n["enrichment"]}

    fig, ax = plt.subplots(figsize=(ONE_COL, 2.6))
    y = np.arange(len(tops))
    for i, k in enumerate(tops):
        r = ctrl.get(k)
        if r:
            hi = min(r["ci_high"], 30)
            ax.plot([max(r["ci_low"], 0.05), hi], [i, i], color=NEUTRAL,
                    linewidth=1.5, solid_capstyle="round", zorder=3)
            ax.scatter([max(r["or_conditional"], 0.05)], [i], s=14,
                       color=NEUTRAL, edgecolor=SURFACE, linewidth=0.6, zorder=4)
        ax.scatter([test[i]], [i], s=20, color=ACCENT, marker="D",
                   edgecolor=SURFACE, linewidth=0.7, zorder=5)
    ax.axvline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(0.04, 40)
    ax.set_xticks([0.1, 1, 10], ["0.1", "1", "10"])
    ax.set_yticks(y, [f"top {k}" for k in tops], fontsize=6.6)
    ax.set_xlabel("Odds ratio (log scale)")
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.scatter([], [], s=20, color=ACCENT, marker="D", label="APOL1-MYH9 locus")
    ax.plot([], [], color=NEUTRAL, linewidth=1.5,
            marker="o", markersize=4, label="control locus, 95% CI")
    ax.legend(loc="lower right", fontsize=6.2, handlelength=1.4)
    ax.set_title("Control intervals include the test-locus estimate",
                 loc="left", pad=6, fontweight="bold")
    save(fig, "figure4")


def main() -> None:
    print("building figures")
    figure1(); figure2(); figure3(); figure4()
    print(f"all figures in {FIGS}")


if __name__ == "__main__":
    main()
