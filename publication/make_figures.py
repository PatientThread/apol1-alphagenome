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
# Gene spans, GRCh38, from Ensembl. Drawn to scale.
GENES = [("APOL4", 36_195_000, 36_206_000, "+"),
         ("APOL2", 36_219_000, 36_233_000, "+"),
         ("APOL1", 36_253_010, 36_267_530, "+"),
         ("MYH9",  36_281_280, 36_388_000, "-")]
RISK = [("G1", 36_265_860), ("G2", 36_265_996)]
CANDS = [("rs132708", 36_197_151), ("rs4820232", 36_314_496),
         ("rs6000250", 36_350_756), ("rs713797", 36_368_489)]
LOCUS_S, LOCUS_E = 36_140_330, 36_388_018


def figure1() -> None:
    """Locus map: what the analysis is actually looking at.

    None of the other figures gives the reader any spatial orientation. This one
    shows where the genes are, where the risk variants sit, where the surviving
    candidates fall, and how the measured accessible regions distribute between
    the two partitions. It is the only figure carrying biology rather than
    statistics.
    """
    import gzip
    def peaks_from(tag_terms):
        iv = []
        for term in tag_terms:
            for f in sorted((ROOT / "data" / "comparison_tissue_peaks").glob("*.bed.gz")):
                pass
        return iv
    j = json.loads((RESULTS / "reconciled.json").read_text())

    # Rebuild the partition intervals from the deposited peak files.
    import glob as _g
    PK = ROOT / "data" / "comparison_tissue_peaks"
    def load(terms):
        iv = []
        for term in terms:
            idx = PK / f"{term.replace(' ', '_')}.index.json"
            if not idx.exists():
                continue
            for f in json.loads(idx.read_text())[:4]:
                fp = PK / f"{f['accession']}.bed.gz"
                if not fp.exists():
                    continue
                with gzip.open(fp, "rt") as fh:
                    for line in fh:
                        c = line.split("\t")
                        if c and c[0] == "chr22":
                            a, b = int(c[1]), int(c[2])
                            if b > LOCUS_S and a < LOCUS_E:
                                iv.append((a, b))
        if not iv:
            return []
        iv.sort(); out = [list(iv[0])]
        for a, b in iv[1:]:
            if a <= out[-1][1]:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a, b])
        return [tuple(x) for x in out]

    kid = load(["kidney", "glomerular visceral epithelial cell"])
    oth = load(["liver", "lung", "heart left ventricle", "stomach",
                "spleen", "thyroid gland"])
    def inside(pt, iv):
        return any(a <= pt < b for a, b in iv)
    kspec = [q for q in kid if not inside((q[0]+q[1])//2, oth)]
    shar = [q for q in kid if inside((q[0]+q[1])//2, oth)]

    fig, ax = plt.subplots(figsize=(TWO_COL, 2.5))
    mb = lambda x: x / 1e6

    rows = {"genes": 3.1, "risk": 2.45, "cand": 2.0,
            "kspec": 1.25, "shared": 0.85, "other": 0.45}

    for name, a, b, strand in GENES:
        ax.add_patch(plt.Rectangle((mb(a), rows["genes"] - 0.11), mb(b) - mb(a), 0.22,
                                   facecolor=NEUTRAL, edgecolor="none", zorder=3))
        ax.text(mb((a + b) / 2), rows["genes"] + 0.20,
                f"{name} {'>' if strand == '+' else '<'}",
                ha="center", fontsize=6.2, color=INK, fontweight="bold")

    for lab, pos in RISK:
        ax.plot([mb(pos), mb(pos)], [rows["risk"] - 0.13, rows["risk"] + 0.13],
                color=ACCENT, linewidth=1.5, zorder=4)
    ax.text(mb(RISK[0][1]) + 0.004, rows["risk"], "G1, G2", fontsize=6.2,
            color=ACCENT, va="center", fontweight="bold")

    # rs6000250 and rs713797 are 18 kb apart and their labels collide at this
    # scale, so alternate labels drop to a second line.
    for k, (lab, pos) in enumerate(CANDS):
        ax.plot([mb(pos), mb(pos)], [rows["cand"] - 0.11, rows["cand"] + 0.11],
                color=INK, linewidth=1.1, zorder=4)
        dy = -0.20 if k % 2 == 0 else -0.34
        ax.plot([mb(pos), mb(pos)], [rows["cand"] - 0.11, rows["cand"] + dy + 0.05],
                color=GRID, linewidth=0.5, zorder=3)
        ax.text(mb(pos), rows["cand"] + dy, lab, fontsize=5.4, color=INK2,
                ha="center")

    for key, iv, lab, col in [("kspec", kspec, f"kidney-restricted ({len(kspec)})", ACCENT),
                              ("shared", shar, f"shared ({len(shar)})", NEUTRAL),
                              ("other", oth, f"comparison tissues ({len(oth)})", GRID)]:
        for a, b in iv:
            ax.add_patch(plt.Rectangle((mb(a), rows[key] - 0.09),
                                       max(mb(b) - mb(a), 0.0006), 0.18,
                                       facecolor=col, edgecolor="none", zorder=3))
        ax.text(mb(LOCUS_S) - 0.004, rows[key], lab, fontsize=6, color=INK2,
                ha="right", va="center")

    ax.text(mb(LOCUS_S) - 0.004, rows["genes"], "genes", fontsize=6,
            color=INK2, ha="right", va="center")
    ax.text(mb(LOCUS_S) - 0.004, rows["risk"], "risk variants", fontsize=6,
            color=INK2, ha="right", va="center")
    ax.text(mb(LOCUS_S) - 0.004, rows["cand"], "candidates", fontsize=6,
            color=INK2, ha="right", va="center")

    ax.set_xlim(mb(LOCUS_S) - 0.055, mb(LOCUS_E) + 0.004)
    ax.set_ylim(0.15, 3.55)
    ax.set_yticks([])
    ax.set_xlabel("chromosome 22 position (Mb, GRCh38)")
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    save(fig, "figure1")


# ---------------------------------------------------------------- figure 2
def figure2() -> None:
    """Enrichment by endpoint, from the reconciled counts."""
    c = pd.read_csv(RESULTS / "reconciled_counts.csv")
    c = c[c.depth == 200]
    eps = ["pooled", "kidney-specific", "shared"]
    titles = ["A  Pooled kidney union", "B  Kidney-restricted", "C  Shared"]
    order = (c[c.endpoint == "pooled"].sort_values("odds_ratio")["ranked_by"].tolist())
    fig, axes = plt.subplots(1, 3, figsize=(TWO_COL, 2.5), sharey=True)
    for ax, ep, title in zip(axes, eps, titles):
        sub = c[c.endpoint == ep].set_index("ranked_by")
        for i, name in enumerate(order):
            if name not in sub.index:
                continue
            r = sub.loc[name]
            is_k = name.startswith("kidney")
            col = ACCENT if is_k else NEUTRAL
            hi = r.ci_high if pd.notna(r.ci_high) else 30
            ax.plot([r.ci_low, hi], [i, i], color=col,
                    linewidth=1.9 if is_k else 1.1, solid_capstyle="round", zorder=3)
            ax.scatter([r.odds_ratio], [i], s=20 if is_k else 13, color=col,
                       marker="D" if is_k else "o", edgecolor=SURFACE,
                       linewidth=0.6, zorder=4)
        ax.axvline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
        ax.set_xscale("log"); ax.set_xlim(0.3, 20)
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
    save(fig, "figure2")


# ---------------------------------------------------------------- figure 3
def figure3() -> None:
    """The interaction: does partitioning change the contrast? It does not."""
    it = pd.read_csv(RESULTS / "reconciled_interaction.csv").sort_values("delta_observed")
    fig, ax = plt.subplots(figsize=(ONE_COL, 2.3))
    for i, r in enumerate(it.itertuples()):
        col = ACCENT if r.excludes_zero else NEUTRAL
        ax.plot([r.delta_ci_low, r.delta_ci_high], [i, i], color=col,
                linewidth=1.6, solid_capstyle="round", zorder=3)
        ax.scatter([r.delta_observed], [i], s=17, color=col, marker="D",
                   edgecolor=SURFACE, linewidth=0.6, zorder=4)
    ax.axvline(0, color=INK2, linewidth=0.9, zorder=2)
    ax.set_yticks(range(len(it)), it["comparator"], fontsize=6.4)
    ax.set_xlabel("Change in podocyte-minus-comparator contrast\n"
                  "on partitioning (natural-log odds ratio)",
                  linespacing=1.6, fontsize=6.8)
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("No interval excludes zero", loc="left", pad=6,
                 fontweight="bold", fontsize=7.4)
    save(fig, "figure3")


# ---------------------------------------------------------------- figure 4
def figure4() -> None:
    """Control locus, with zero-event groups labelled rather than plotted."""
    n = json.loads((RESULTS / "negative_control_locus.json").read_text())
    e = json.loads((RESULTS / "encode_validation.json").read_text())
    t = e["support_1"]["tests"]
    tops = [25, 50, 100, 200]
    test = [t[f"dnase_glomerular_visceral_epithelial_cell_top{k}"]["odds_ratio"]
            for k in tops]
    ctrl = {r["top_n"]: r for r in n["enrichment"]}
    fig, ax = plt.subplots(figsize=(ONE_COL, 2.4))
    for i, k in enumerate(tops):
        r = ctrl.get(k)
        if r and r.get("in_peak", 1) == 0:
            ax.text(0.075, i, f"0 of {k} in peak", fontsize=5.8, color=INK2,
                    va="center", ha="left", style="italic")
            ax.annotate("", xy=(min(r["ci_high"], 25), i), xytext=(0.9, i),
                        arrowprops=dict(arrowstyle="->", color=NEUTRAL, lw=1.0))
        elif r:
            ax.plot([r["ci_low"], min(r["ci_high"], 25)], [i, i], color=NEUTRAL,
                    linewidth=1.5, solid_capstyle="round", zorder=3)
            ax.scatter([r["or_conditional"]], [i], s=14, color=NEUTRAL,
                       edgecolor=SURFACE, linewidth=0.6, zorder=4)
        ax.scatter([test[i]], [i], s=20, color=ACCENT, marker="D",
                   edgecolor=SURFACE, linewidth=0.7, zorder=5)
    ax.axvline(1, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
    ax.set_xscale("log"); ax.set_xlim(0.06, 40)
    ax.set_xticks([0.1, 1, 10], ["0.1", "1", "10"])
    ax.set_yticks(range(len(tops)), [f"top {k}" for k in tops], fontsize=6.6)
    ax.set_xlabel("Odds ratio (log scale)")
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.scatter([], [], s=20, color=ACCENT, marker="D", label="APOL1-MYH9 estimate")
    ax.plot([], [], color=NEUTRAL, linewidth=1.5, marker="o", markersize=4,
            label="control locus, 95% CI")
    ax.legend(loc="lower right", fontsize=6, handlelength=1.4)
    ax.set_title("Control intervals include the test-locus estimate",
                 loc="left", pad=6, fontweight="bold")
    save(fig, "figure4")


def main() -> None:
    print("building figures")
    figure1(); figure2(); figure3(); figure4()
    print(f"all figures in {FIGS}")


if __name__ == "__main__":
    main()
