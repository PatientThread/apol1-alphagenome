#!/usr/bin/env python3
"""
03_analyse.py

APOL1 paper, step 3: what the predictions say, and what they do not.

FOUR ANALYSES, in increasing order of how much they claim.

  A. DISTRIBUTION. How large are predicted regulatory effects across the locus,
     and is there structure, or is everything near zero? If the latter, that is
     the result and the paper is short.

  B. NEGATIVE CONTROL, and this one decides whether anything else can be
     believed. G1 and G2 are CODING variants: two missense substitutions and a
     six-base deletion in the last exon. They change the protein, not its
     regulation. A sequence model that assigns them large REGULATORY effects is
     telling us about coding sequence constraint leaking into a regulatory
     prediction, not about expression. So the risk variants should sit in the
     BULK of the regulatory effect distribution, not at its extreme. If they do
     not, every downstream ranking is suspect and must be reported as such.

  C. RANKING. Which non-coding variants are predicted to have the largest effect
     on APOL1 expression in kidney, and where do they sit relative to the gene?
     Clustering at the promoter and 5' region would be mechanistically coherent;
     a flat scatter across 248 kb would not.

  D. PODOCYTE vs BULK KIDNEY. APOL1 kills podocytes. Does the podocyte
     accessibility track rank variants differently from bulk kidney cortex
     expression? A variant that matters in podocyte and not in bulk tissue is
     exactly what a bulk eQTL study would miss, and is the argument for using a
     cell-type-resolved prediction at all.

WHAT THIS SCRIPT DELIBERATELY DOES NOT DO. It does not claim any variant is
causal, and it does not test the G1/G2 haplotype-background hypothesis, which
needs phased reference haplotypes rather than allele frequencies. That is the
next step and is stated as such.

Output: results/analysis.json, results/top_variants.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCORED = RESULTS / "scored"

APOL1_START, APOL1_END = 36_253_010, 36_267_530
APOL1_TSS = APOL1_START            # gene is on the + strand

RISK_RSIDS = {"rs73885319", "rs60910145", "rs71785313", "rs2239785"}

RNA_KIDNEY = "rna_cortex_of_kidney"
RNA_MEDULLA = "rna_outer_medulla_of_kidney"
PODO = "dnase_glomerular_visceral_epithelial_cell"
ATAC_KIDNEY = "atac_kidney"


def load() -> pd.DataFrame:
    files = sorted(SCORED.glob("*.tsv.gz"))
    if not files:
        raise SystemExit("no scored chunks yet")
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    return d.drop_duplicates(subset=["variant_id"])


def pct_rank(s: pd.Series, v: float) -> float:
    return float((s.abs() < abs(v)).mean() * 100)


def main() -> None:
    d = load()
    n = len(d)
    have = [c for c in (RNA_KIDNEY, RNA_MEDULLA, PODO, ATAC_KIDNEY)
            if c in d.columns]
    print(f"loaded {n:,} scored variants; key tracks present: {have}\n")

    coding_terms = ("missense", "synonymous", "stop", "frameshift", "inframe",
                    "start_lost", "splice")
    d["is_coding"] = d["consequence"].fillna("").str.contains(
        "|".join(coding_terms), case=False)
    d["dist_to_tss"] = d["pos"] - APOL1_TSS
    d["in_gene"] = d["pos"].between(APOL1_START, APOL1_END)

    out: dict = {"variants_scored": n,
                 "coding": int(d["is_coding"].sum()),
                 "non_coding": int((~d["is_coding"]).sum())}

    # ---------------------------------------------------------------- A
    print("=" * 74)
    print("A. DISTRIBUTION OF PREDICTED EFFECTS")
    print("=" * 74)
    dist = {}
    for c in have:
        s = d[c].dropna()
        dist[c] = {"n": int(len(s)), "median_abs": float(s.abs().median()),
                   "p95_abs": float(s.abs().quantile(0.95)),
                   "max_abs": float(s.abs().max()),
                   "frac_above_0.01": float((s.abs() > 0.01).mean())}
        print(f"  {c:<46} median|x| {s.abs().median():.5f}  "
              f"p95 {s.abs().quantile(0.95):.4f}  max {s.abs().max():.4f}")
    out["distribution"] = dist

    # ---------------------------------------------------------------- B
    print()
    print("=" * 74)
    print("B. NEGATIVE CONTROL: WHERE DO THE CODING RISK VARIANTS RANK?")
    print("=" * 74)
    print("  G1/G2 are coding. They should NOT be regulatory outliers.")
    risk = d[d["rsid"].isin(RISK_RSIDS)]
    ctrl = {}
    if len(risk):
        for c in have:
            s = d[c].dropna()
            per = {}
            for _, r in risk.iterrows():
                if pd.isna(r.get(c)):
                    continue
                per[r["rsid"]] = {"value": float(r[c]),
                                  "percentile_of_abs": round(pct_rank(s, r[c]), 1)}
            ctrl[c] = per
        for c in have:
            print(f"  {c}")
            for rs, v in ctrl.get(c, {}).items():
                flag = "  <-- EXTREME, investigate" if v["percentile_of_abs"] > 95 else ""
                print(f"    {rs:<13} {v['value']:+.5f}   "
                      f"{v['percentile_of_abs']:>5.1f}th pct of |effect|{flag}")
        extreme = any(v["percentile_of_abs"] > 95
                      for c in ctrl for v in ctrl[c].values())
        out["negative_control"] = {"detail": ctrl, "any_extreme": bool(extreme)}
        print()
        if extreme:
            print("  WARNING: at least one coding risk variant is a regulatory")
            print("  outlier. Coding constraint may be leaking into the")
            print("  regulatory prediction. Report this prominently; it bounds")
            print("  what any ranking below can be said to mean.")
        else:
            print("  PASS: the coding risk variants sit in the bulk of the")
            print("  distribution, as they should. The regulatory channel is")
            print("  not simply re-reading protein-coding constraint.")
    else:
        print("  risk variants not present in the scored set (check AF filter)")

    # ---------------------------------------------------------------- C
    print()
    print("=" * 74)
    print("C. TOP NON-CODING VARIANTS BY PREDICTED APOL1 EXPRESSION EFFECT")
    print("=" * 74)
    if RNA_KIDNEY in d.columns:
        nc = d[~d["is_coding"] & d[RNA_KIDNEY].notna()].copy()
        nc["abs_rna"] = nc[RNA_KIDNEY].abs()
        top = nc.nlargest(20, "abs_rna")
        cols = ["rsid", "pos", "dist_to_tss", "in_gene", "af_afr",
                RNA_KIDNEY] + ([PODO] if PODO in d.columns else [])
        top[cols].to_csv(RESULTS / "top_variants.csv", index=False)
        print(f"  {'rsid':<13}{'pos':>11}{'d(TSS)':>9}{'inGene':>7}"
              f"{'AF_afr':>8}{'RNA kidney':>12}")
        for _, r in top.head(12).iterrows():
            print(f"  {str(r['rsid'])[:12]:<13}{int(r['pos']):>11,}"
                  f"{int(r['dist_to_tss']):>9,}{str(r['in_gene']):>7}"
                  f"{r['af_afr']:>8.3f}{r[RNA_KIDNEY]:>+12.5f}")
        in_gene_frac = float(top["in_gene"].mean())
        out["top20_fraction_in_gene"] = in_gene_frac
        print()
        print(f"  fraction of top 20 inside the APOL1 gene body: {in_gene_frac:.0%}")
        print(f"  (gene is {(APOL1_END-APOL1_START)/1000:.0f} kb of a "
              f"{247688/1000:.0f} kb locus = "
              f"{100*(APOL1_END-APOL1_START)/247688:.0f}% by chance)")

    # ---------------------------------------------------------------- D
    print()
    print("=" * 74)
    print("D. PODOCYTE vs BULK KIDNEY")
    print("=" * 74)
    if PODO in d.columns and RNA_KIDNEY in d.columns:
        sub = d[[PODO, RNA_KIDNEY]].dropna()
        if len(sub) > 10:
            from scipy import stats
            rho, p = stats.spearmanr(sub[PODO].abs(), sub[RNA_KIDNEY].abs())
            out["podocyte_vs_bulk_spearman"] = {"rho": float(rho),
                                                "p": float(p), "n": int(len(sub))}
            print(f"  rank correlation of |podocyte accessibility| with "
                  f"|bulk kidney expression|")
            print(f"    rho = {rho:.3f}  (n = {len(sub):,})")
            print()
            if abs(rho) < 0.3:
                print("  The two rank variants LARGELY DIFFERENTLY. A podocyte-")
                print("  specific signal is not recoverable from bulk kidney,")
                print("  which is the argument for cell-type-resolved prediction")
                print("  and the reason a bulk eQTL study would miss it.")
            else:
                print("  The two are substantially concordant, so the podocyte")
                print("  track adds less independent information than hoped.")

    (RESULTS / "analysis.json").write_text(json.dumps(out, indent=2, default=str))
    print()
    print(f"  Wrote {RESULTS}/analysis.json and top_variants.csv")


if __name__ == "__main__":
    main()
