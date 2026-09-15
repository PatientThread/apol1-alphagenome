#!/usr/bin/env python3
"""
12_analyse_wrong_tissue.py

APOL1 paper, step 12: read out the wrong-tissue control.

Step 10 scored the same 1,930 variants and kept, for each one, the predicted
accessibility effect in kidney AND in eight deliberately non-renal biosamples,
all from the same API call. This script ranks the variants separately by each of
those columns and asks the same question of each ranking:

    how enriched are its top-ranked variants in MEASURED KIDNEY open chromatin?

Only the ranking column changes. Same variants, same measured peaks, same
locus-matched background, same test. So any difference between rows is
attributable to the tissue label and to nothing else.

READING THE RESULT. If the kidney rows sit clearly above the liver, brain, lung,
heart and stomach rows, the label carries information. If they do not, the model
is finding regulatory positions rather than kidney regulatory positions, and the
enrichment reported in step 5 cannot be described as a kidney finding.

Output: results/wrong_tissue_enrichment.csv
        results/wrong_tissue_control.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
SCORED = RESULTS / "scored_wrong_tissue"
PEAKS = DATA / "encode_peaks"

CHROM = "chr22"
LOCUS_START, LOCUS_END = 36_140_330, 36_388_018
TOP_N = (25, 50, 100, 200)
KIDNEY_COLS = ["dnase_glomerular_visceral_epithelial_cell", "dnase_kidney"]


def peak_union() -> list[tuple[int, int]]:
    iv = []
    for p in sorted(PEAKS.glob("*.bed.gz")):
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.split("\t")
                if f and f[0] == CHROM:
                    s, e = int(f[1]), int(f[2])
                    if e > LOCUS_START and s < LOCUS_END:
                        iv.append((s, e))
    iv.sort()
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def in_any(pos: int, iv: list[tuple[int, int]]) -> bool:
    lo, hi = 0, len(iv) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        if pos < iv[m][0]:
            hi = m - 1
        elif pos >= iv[m][1]:
            lo = m + 1
        else:
            return True
    return False


def main() -> None:
    d = pd.concat([pd.read_csv(f, sep="\t")
                   for f in sorted(SCORED.glob("chunk_*.tsv.gz"))],
                  ignore_index=True).drop_duplicates("rsid")
    cols = [c for c in d.columns if c.startswith("dnase_")]
    print(f"variants: {len(d):,}   tissue columns: {len(cols)}")

    peaks = peak_union()
    d["hit"] = d["pos"].map(lambda p: in_any(int(p), peaks))
    print(f"measured kidney peaks on the locus: {len(peaks)}, "
          f"{sum(e - s for s, e in peaks):,} bp\n")

    rows = []
    for c in cols:
        s = d.dropna(subset=[c]).copy()
        s["k"] = s[c].abs()
        s = s.sort_values("k", ascending=False).reset_index(drop=True)
        hit = s["hit"].values
        for n in TOP_N:
            top, rest = hit[:n], hit[n:]
            a, b = int(top.sum()), int(n - top.sum())
            cc, dd = int(rest.sum()), int(len(rest) - rest.sum())
            orr, p = stats.fisher_exact([[a, b], [cc, dd]], alternative="greater")
            rows.append({"ranked_by": c.replace("dnase_", ""),
                         "is_kidney": c in KIDNEY_COLS,
                         "top_n": n, "in_peak": a,
                         "pct": round(100 * a / n, 1),
                         "odds_ratio": round(float(orr), 3), "p": float(p)})

    e = pd.DataFrame(rows)
    e.to_csv(RESULTS / "wrong_tissue_enrichment.csv", index=False)

    piv = e.pivot(index="ranked_by", columns="top_n", values="odds_ratio")
    kid = sorted({c.replace("dnase_", "") for c in KIDNEY_COLS} & set(piv.index))
    wrong = [i for i in piv.index if i not in kid]
    order = kid + sorted(wrong, key=lambda i: -piv.loc[i, TOP_N[-1]])

    print("=" * 74)
    print("ENRICHMENT IN MEASURED KIDNEY CHROMATIN, BY THE TISSUE USED TO RANK")
    print("=" * 74)
    print(f"  {'ranked by':<34}" + "".join(f"{'top'+str(n):>9}" for n in TOP_N))
    for i in order:
        mark = "  <- kidney" if i in kid else ""
        print(f"  {i:<34}" + "".join(f"{piv.loc[i, n]:>9.2f}" for n in TOP_N) + mark)

    best_k = {n: max(piv.loc[i, n] for i in kid) for n in TOP_N}
    best_w = {n: max(piv.loc[i, n] for i in wrong) for n in TOP_N}
    med_w = {n: piv.loc[wrong, n].median() for n in TOP_N}

    print("\n  best kidney vs best wrong tissue:")
    for n in TOP_N:
        print(f"    top {n:<5} kidney {best_k[n]:>6.2f}   "
              f"best wrong {best_w[n]:>6.2f}   median wrong {med_w[n]:>6.2f}")

    kidney_wins = sum(1 for n in TOP_N if best_k[n] > best_w[n])
    print("\n  READING:")
    if kidney_wins == len(TOP_N):
        print("  The kidney ranking beats every wrong tissue at every threshold.")
        print("  The tissue label carries real information and the step 5")
        print("  enrichment can be described as a kidney result.")
    elif kidney_wins == 0:
        print("  A wrong tissue matches or beats kidney at EVERY threshold. The")
        print("  label carries no usable information here: the model is finding")
        print("  regulatory positions, not kidney regulatory positions. Step 5")
        print("  must be reworded and the cell-type language dropped entirely.")
    else:
        print(f"  Kidney leads at {kidney_wins} of {len(TOP_N)} thresholds. Too weak")
        print("  to support a kidney-specific claim; report the whole table.")

    out = {"variants": int(len(d)), "peaks": len(peaks),
           "kidney_columns": kid, "wrong_columns": wrong,
           "kidney_leads_at_n_thresholds": kidney_wins,
           "n_thresholds": len(TOP_N),
           "best_kidney_or": {str(n): float(best_k[n]) for n in TOP_N},
           "best_wrong_or": {str(n): float(best_w[n]) for n in TOP_N},
           "median_wrong_or": {str(n): float(med_w[n]) for n in TOP_N}}
    (RESULTS / "wrong_tissue_control.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/wrong_tissue_control.json")


if __name__ == "__main__":
    main()
