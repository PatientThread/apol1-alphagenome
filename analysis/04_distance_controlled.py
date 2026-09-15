#!/usr/bin/env python3
"""
04_distance_controlled.py

APOL1 paper, step 4: redo the ranking with the proximity confound removed.

WHY THIS SCRIPT EXISTS. Script 03 found that the gene-masked expression scorer
assigns systematically larger predicted effects to variants inside the target
gene body:

    coding                        n=  65   median |LFC| 0.00109
    non-coding, inside the gene   n= 121   median |LFC| 0.00199
    non-coding, outside the gene  n=1744   median |LFC| 0.00077

Note the ordering. Non-coding variants INSIDE the gene score higher than coding
ones, so this is not coding constraint leaking into a regulatory prediction. It
is position. That is reasonable behaviour for a scorer that masks on a gene, but
it means a raw ranking by predicted expression effect largely rediscovers
proximity to APOL1, and the apparent result that 90% of top-ranked non-coding
variants lie inside the gene body (against 6% expected by chance) is an artefact.

THREE CORRECTIONS APPLIED HERE.

  1. STRATIFY. Rank within gene-body and outside-gene-body strata separately,
     so variants only compete against others at comparable distance.

  2. RESIDUALISE. Regress log|effect| on log distance to the gene and rank on
     the residual. A variant that scores far above what its distance predicts is
     the interesting one; a variant that scores high because it sits in the gene
     is not.

  3. LEAD WITH ACCESSIBILITY. The DNase and ATAC scorers are position-local and
     showed the OPPOSITE pattern (coding variants score slightly lower), so they
     do not carry this confound. They become the primary regulatory readout and
     the expression channel is reported as supporting evidence with its
     limitation stated.

The negative control is re-run after each correction. If the coding risk
variants stop being outliers once distance is controlled, the correction worked
and the residual ranking can be interpreted. If they remain outliers, something
else is wrong and the ranking must not be used at all.

Output: results/distance_controlled.json, results/top_variants_corrected.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

APOL1_START, APOL1_END = 36_253_010, 36_267_530
RISK = {"rs73885319": "G1 S342G", "rs60910145": "G1 I384M",
        "rs71785313": "G2 del", "rs2239785": "N264K protective"}

RNA = "rna_cortex_of_kidney"
PODO = "dnase_glomerular_visceral_epithelial_cell"
ATAC = "atac_kidney"


def dist_to_gene(pos: int) -> int:
    if APOL1_START <= pos <= APOL1_END:
        return 0
    return min(abs(pos - APOL1_START), abs(pos - APOL1_END))


def pct_of_abs(series: pd.Series, value: float) -> float:
    return float((series.abs() < abs(value)).mean() * 100)


def main() -> None:
    files = sorted(glob.glob(str(RESULTS / "scored" / "*.tsv.gz")))
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    d = d.drop_duplicates(subset=["variant_id"])
    ct = ("missense", "synonymous", "stop", "frameshift", "inframe",
          "start_lost", "splice")
    d["coding"] = d["consequence"].fillna("").str.contains("|".join(ct), case=False)
    d["dist"] = d["pos"].map(dist_to_gene)
    d["in_gene"] = d["dist"] == 0
    print(f"loaded {len(d):,} variants\n")

    out: dict = {"n": int(len(d))}

    # ------------------------------------------------------------ 1 stratify
    print("=" * 76)
    print("CORRECTION 1: RANK WITHIN DISTANCE STRATA")
    print("=" * 76)
    strata = {}
    for name, sub in (("inside_gene", d[d["in_gene"]]),
                      ("outside_gene", d[~d["in_gene"]])):
        s = sub[sub[RNA].notna()]
        strata[name] = {"n": int(len(s)),
                        "median_abs": float(s[RNA].abs().median()),
                        "p95_abs": float(s[RNA].abs().quantile(0.95))}
        print(f"  {name:<14} n={len(s):>5}  median|LFC| {s[RNA].abs().median():.6f}"
              f"  p95 {s[RNA].abs().quantile(0.95):.5f}")
    out["strata"] = strata

    print("\n  negative control, WITHIN the gene-body stratum only:")
    ing = d[d["in_gene"] & d[RNA].notna()]
    ctrl_strat = {}
    for rs, label in RISK.items():
        row = ing[ing["rsid"] == rs]
        if len(row):
            v = float(row.iloc[0][RNA])
            p = pct_of_abs(ing[RNA], v)
            ctrl_strat[rs] = {"value": v, "pct_within_stratum": round(p, 1)}
            flag = "  <-- still extreme" if p > 95 else ""
            print(f"    {rs:<13} {label:<18} {v:+.5f}  {p:>5.1f}th pct{flag}")
    out["negative_control_stratified"] = ctrl_strat
    still = [k for k, v in ctrl_strat.items() if v["pct_within_stratum"] > 95]
    print()
    if still:
        print(f"  STILL EXTREME within stratum: {still}")
        print("  Distance alone does not explain it. Do not rank on this channel.")
    else:
        print("  RESOLVED: once variants compete only against others at the same")
        print("  distance, the coding risk variants are unremarkable. The raw")
        print("  ranking was measuring proximity, exactly as suspected.")

    # ---------------------------------------------------------- 2 residualise
    print()
    print("=" * 76)
    print("CORRECTION 2: RESIDUAL AFTER REGRESSING ON DISTANCE")
    print("=" * 76)
    r = d[d[RNA].notna()].copy()
    r["log_abs"] = np.log10(r[RNA].abs() + 1e-9)
    r["log_dist"] = np.log10(r["dist"] + 1)
    sl, ic, rv, pv, se = stats.linregress(r["log_dist"], r["log_abs"])
    r["residual"] = r["log_abs"] - (ic + sl * r["log_dist"])
    print(f"  log|effect| ~ log(distance):  slope {sl:+.4f}, r = {rv:.3f}, "
          f"p = {pv:.2e}")
    print(f"  distance explains {100*rv**2:.1f}% of the variance in log|effect|")
    out["distance_regression"] = {"slope": float(sl), "r": float(rv),
                                  "r2": float(rv ** 2), "p": float(pv)}

    print("\n  negative control on the RESIDUAL:")
    ctrl_res = {}
    for rs, label in RISK.items():
        row = r[r["rsid"] == rs]
        if len(row):
            v = float(row.iloc[0]["residual"])
            p = float((r["residual"] < v).mean() * 100)
            ctrl_res[rs] = {"residual": v, "pct": round(p, 1)}
            flag = "  <-- still extreme" if p > 95 else ""
            print(f"    {rs:<13} {label:<18} residual {v:+.3f}  {p:>5.1f}th pct{flag}")
    out["negative_control_residual"] = ctrl_res

    # ------------------------------------------------- 3 accessibility primary
    print()
    print("=" * 76)
    print("CORRECTION 3: ACCESSIBILITY AS PRIMARY READOUT")
    print("=" * 76)
    for col, lab in ((PODO, "podocyte DNase"), (ATAC, "kidney ATAC")):
        if col not in d.columns:
            continue
        s = d[d[col].notna()]
        sl2, ic2, rv2, pv2, _ = stats.linregress(
            np.log10(s["dist"] + 1), np.log10(s[col].abs() + 1e-9))
        print(f"  {lab:<16} distance explains {100*rv2**2:>4.1f}% of variance "
              f"(r={rv2:+.3f})")
        out[f"distance_r2_{col}"] = float(rv2 ** 2)

    print()
    print("  TOP NON-CODING VARIANTS BY PODOCYTE ACCESSIBILITY EFFECT")
    if PODO in d.columns:
        nc = d[~d["coding"] & d[PODO].notna()].copy()
        nc["abs_podo"] = nc[PODO].abs()
        top = nc.nlargest(15, "abs_podo")
        cols = ["rsid", "pos", "dist", "in_gene", "af_afr", PODO, RNA]
        top[[c for c in cols if c in top.columns]].to_csv(
            RESULTS / "top_variants_corrected.csv", index=False)
        print(f"  {'rsid':<14}{'pos':>11}{'dist':>9}{'AF_afr':>8}"
              f"{'podocyte':>11}{'RNA':>10}")
        for _, x in top.head(10).iterrows():
            print(f"  {str(x['rsid'])[:13]:<14}{int(x['pos']):>11,}"
                  f"{int(x['dist']):>9,}{x['af_afr']:>8.3f}"
                  f"{x[PODO]:>+11.5f}{x[RNA]:>+10.5f}"
                  if pd.notna(x.get(RNA)) else "")
        print()
        print(f"  fraction of top 15 inside the gene body: "
              f"{top['in_gene'].mean():.0%}  (6% expected by chance)")
        out["podocyte_top15_in_gene"] = float(top["in_gene"].mean())

    (RESULTS / "distance_controlled.json").write_text(json.dumps(out, indent=2))
    print()
    print(f"  Wrote {RESULTS}/distance_controlled.json")


if __name__ == "__main__":
    main()
