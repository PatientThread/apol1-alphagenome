#!/usr/bin/env python3
"""
20_shared_draw_interaction.py

Rebuilds the change-in-contrast analysis on ONE set of bootstrap draws shared by
every comparator, and recomputes locus coverage from the intersection with the
evaluated interval. Both changes were requested at review.

WHY THE DRAWS ARE NOW SHARED. Step 16 paired the two reference definitions
within a comparator but drew fresh peaks for each comparator, so the six deltas
were not estimated on common draws and were not comparable with each other.
Here a single list of B draws is generated once and every comparator is
evaluated on each draw.

THE ALGORITHM, stated without ambiguity, because the previous wording conflated
two different things (which records exist, and which of them land in the top k):

  for each replicate:
      draw 39 kidney-restricted peaks WITH REPLACEMENT (39 = the number of
          such peaks; a peak drawn k times contributes its allele records k
          times)
      build the replicate allele pool = those copied records
                                      + every non-restricted record, once
      for each tissue output:
          rerank the WHOLE replicate pool by |score| for that tissue
          take the top 200 records of that reranking
          form the 2x2 against membership of the reference set in question
  The fixed part of the pool is fixed in its RECORDS AND SCORES. It is not
  fixed in contingency membership: an allele in the observed top 200 need not
  survive reranking of a replicate pool, and frequently does not.

The observed point estimate and every replicate use the SAME estimator, a
Haldane-corrected log cross-product ratio. That estimator is deliberately not
the conditional maximum likelihood odds ratio used in Table 1; the two are
reported for different purposes and are not interchangeable.

Output: results/shared_draw_interaction.csv
        results/shared_draw_stability.csv
        results/coverage_clipped.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import argparse
import glob
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA, RESULTS = ROOT / "data", ROOT / "results"
TEST = ("chr22", 36_140_330, 36_388_018)
CTRL = ("chr11", 5_150_000, 5_397_688)
PRIMARY, B = 200, 500
SEED = 20260921
PCTS = (2.5, 97.5)

RANKINGS = {
    "kidney (podocyte)": "dnase_glomerular_visceral_epithelial_cell",
    "kidney (whole)": "dnase_kidney", "hepatocyte": "dnase_hepatocyte",
    "liver": "dnase_liver", "lung": "dnase_lung", "brain": "dnase_brain",
    "stomach": "dnase_stomach",
}
REF = "kidney (podocyte)"
KIDNEY_TERMS = ["kidney", "glomerular visceral epithelial cell"]
COMPARISON = ["liver", "lung", "heart_left_ventricle", "stomach", "spleen",
              "thyroid_gland"]
KO = ["ENCFF473YKR", "ENCFF554CKH", "ENCFF350FBW", "ENCFF214HJG",
      "ENCFF851SWR", "ENCFF922WKR", "ENCFF948WNJ", "ENCFF146AYH"]
KU = ["ENCFF473YKR", "ENCFF554CKH", "ENCFF350FBW", "ENCFF214HJG",
      "ENCFF618RWK", "ENCFF325DZQ", "ENCFF019PKU", "ENCFF971ABV"]
PEAK_DIRS = ["encode_peaks", "encode_celltype_peaks", "comparison_tissue_peaks",
             "control_locus_peaks", "control_locus_peaks_matched"]


def merge(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def find(acc):
    for d in PEAK_DIRS:
        p = DATA / d / f"{acc}.bed.gz"
        if p.exists():
            return p
    raise FileNotFoundError(acc)


def load(accs, locus):
    ch, lo, hi = locus
    iv = []
    for a in accs:
        with gzip.open(find(a), "rt") as fh:
            for line in fh:
                c = line.split("\t")
                if c and c[0] == ch:
                    s, e = int(c[1]), int(c[2])
                    if e > lo and s < hi:
                        iv.append((s, e))
    return merge(iv)


def which(pos, iv):
    lo, hi = 0, len(iv) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        if pos < iv[m][0]:
            hi = m - 1
        elif pos >= iv[m][1]:
            lo = m + 1
        else:
            return m
    return None


def logor(a, b, c, d):
    return float(np.log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))))


def clipped_bp(iv, locus):
    _, lo, hi = locus
    return sum(min(e, hi) - max(s, lo) for s, e in iv)


def coverage() -> None:
    """Report peaks classified unclipped, coverage measured clipped."""
    rows = []
    for nm, accs, locus, label in (
            ("APOL1-MYH9", KU, TEST, "podocyte-inclusive union"),
            ("APOL1-MYH9", KO, TEST, "whole kidney"),
            ("beta-globin", KU, CTRL, "podocyte-inclusive union"),
            ("beta-globin", KO, CTRL, "whole kidney")):
        m = load(accs, locus)
        span = locus[2] - locus[1]
        unc, cl = sum(e - s for s, e in m), clipped_bp(m, locus)
        rows.append({"locus": nm, "reference_set": label, "peaks": len(m),
                     "bp_unclipped": unc, "bp_within_locus": cl,
                     "overhang_bp": unc - cl,
                     "peaks_crossing_boundary": sum(
                         1 for s, e in m if s < locus[1] or e > locus[2]),
                     "pct_of_locus": round(100 * cl / span, 1)})
    # partitions of the union at the test locus
    kid = load(KU, TEST)
    comp = merge([iv for t in COMPARISON
                  for iv in load([f["accession"] for f in json.loads(
                      (DATA / "comparison_tissue_peaks" /
                       f"{t}.index.json").read_text())], TEST)])
    ks = [p for p in kid if which((p[0] + p[1]) // 2, comp) is None]
    sh = [p for p in kid if which((p[0] + p[1]) // 2, comp) is not None]
    span = TEST[2] - TEST[1]
    for lbl, part in (("union, kidney-restricted", ks), ("union, shared", sh)):
        unc, cl = sum(e - s for s, e in part), clipped_bp(part, TEST)
        rows.append({"locus": "APOL1-MYH9", "reference_set": lbl,
                     "peaks": len(part), "bp_unclipped": unc,
                     "bp_within_locus": cl, "overhang_bp": unc - cl,
                     "peaks_crossing_boundary": sum(
                         1 for s, e in part if s < TEST[1] or e > TEST[2]),
                     "pct_of_locus": round(100 * cl / span, 1)})
    cv = pd.DataFrame(rows)
    cv.to_csv(RESULTS / "coverage_clipped.csv", index=False)
    print(cv.to_string(index=False))
    tot = cv[cv.reference_set.str.startswith("union,")]["bp_within_locus"].sum()
    par = cv[(cv.locus == "APOL1-MYH9") &
             (cv.reference_set == "podocyte-inclusive union")
             ]["bp_within_locus"].iloc[0]
    print(f"\npartitions sum to parent on clipped bp: {tot} == {par} -> "
          f"{'OK' if tot == par else 'MISMATCH'}")
    return ks, kid, comp


def main(stability: bool = False) -> None:
    ks, kid, comp = coverage()

    d = pd.concat([pd.read_csv(f, sep="\t") for f in
                   sorted(glob.glob(str(RESULTS / "scored_wrong_tissue" /
                                        "*.tsv.gz")))], ignore_index=True)
    d = d.drop_duplicates(["pos", "af_afr"]).reset_index(drop=True)
    d["in_kidney"] = d["pos"].map(lambda p: which(int(p), kid) is not None)
    d["in_kspec"] = d["pos"].map(lambda p: which(int(p), ks) is not None)

    elem = d["pos"].map(lambda p: which(int(p), ks))
    by_elem = {e: g.index.to_numpy() for e, g in d.groupby(elem)}
    outside = d.index[elem.isna()].to_numpy()
    elems = list(by_elem)
    print(f"\nresampling unit: {len(elems)} kidney-restricted peaks; "
          f"{len(outside)} records held fixed; {B} replicates; seed {SEED}")

    def contrast(fr, col_a, col_b, flag, n=PRIMARY):
        out = []
        for col in (col_a, col_b):
            s = fr.dropna(subset=[col]).copy()
            s["k"] = s[col].abs()
            s = s.sort_values("k", ascending=False)
            h = s[flag].values
            a, b = int(h[:n].sum()), int(n - h[:n].sum())
            c, dd = int(h[n:].sum()), int(len(h) - n - h[n:].sum())
            out.append(logor(a, b, c, dd))
        return out[0] - out[1]

    def deltas_on(fr):
        return {r: (contrast(fr, RANKINGS[REF], RANKINGS[r], "in_kspec")
                    - contrast(fr, RANKINGS[REF], RANKINGS[r], "in_kidney"))
                for r in RANKINGS if r != REF}

    obs = deltas_on(d)
    # ONE list of draws, every comparator evaluated on each of them.
    rng = np.random.default_rng(SEED)
    draws = [rng.choice(elems, size=len(elems), replace=True) for _ in range(B)]
    acc = {r: [] for r in obs}
    for j, pick in enumerate(draws):
        idx = np.concatenate([by_elem[e] for e in pick] + [outside])
        for r, v in deltas_on(d.loc[idx]).items():
            acc[r].append(v)
        if (j + 1) % 100 == 0:
            print(f"  {j + 1}/{B}")

    if stability:
        # Nine configurations (three replicate counts x three seeds), each
        # yielding six comparator-specific intervals, so 54 intervals in all.
        srows = []
        for Bn in (500, 2000, 5000):
            for sd in (20260921, 20260922, 20260923):
                r2 = np.random.default_rng(sd)
                a2 = {r: [] for r in obs}
                for _ in range(Bn):
                    pk = r2.choice(elems, size=len(elems), replace=True)
                    ix = np.concatenate([by_elem[e] for e in pk] + [outside])
                    for r, v in deltas_on(d.loc[ix]).items():
                        a2[r].append(v)
                for r, vals in a2.items():
                    lo, hi = np.percentile(vals, PCTS)
                    srows.append({"comparator": r, "replicates": Bn,
                                  "seed": sd, "ci_low": round(float(lo), 4),
                                  "ci_high": round(float(hi), 4),
                                  "excludes_zero": bool(lo > 0 or hi < 0)})
                print(f"  stability B={Bn} seed={sd}", flush=True)
        st = pd.DataFrame(srows)
        st.to_csv(RESULTS / "shared_draw_stability.csv", index=False)
        g = st.groupby("comparator").agg(
            lo_lo=("ci_low", "min"), lo_hi=("ci_low", "max"),
            hi_lo=("ci_high", "min"), hi_hi=("ci_high", "max"))
        g["max_shift"] = np.maximum(g.lo_hi - g.lo_lo, g.hi_hi - g.hi_lo).round(3)
        print(g.round(3).to_string())
        print(f"configurations {st[['replicates','seed']].drop_duplicates().shape[0]}"
              f" | comparator intervals {len(st)}"
              f" | any excludes zero {bool(st.excludes_zero.any())}"
              f" | largest endpoint shift {g.max_shift.max()}")

    rows = []
    for r, vals in acc.items():
        lo, hi = np.percentile(vals, PCTS)
        rows.append({"comparator": r, "delta_observed": round(obs[r], 3),
                     "delta_ci_low": round(float(lo), 3),
                     "delta_ci_high": round(float(hi), 3),
                     "excludes_zero": bool(lo > 0 or hi < 0),
                     "replicates": B, "seed": SEED,
                     "pct_low": PCTS[0], "pct_high": PCTS[1],
                     "estimator": "Haldane-corrected log cross-product ratio",
                     "draws_shared_across_comparators": True})
    it = pd.DataFrame(rows).sort_values("delta_observed", ascending=False)
    it.to_csv(RESULTS / "shared_draw_interaction.csv", index=False)
    print("\nCHANGE IN CONTRAST, common draws across all six comparators")
    print(it[["comparator", "delta_observed", "delta_ci_low", "delta_ci_high",
              "excludes_zero"]].to_string(index=False))
    print(f"\nany interval excluding zero: {bool(it.excludes_zero.any())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stability", action="store_true",
                    help="also run the nine replicate-count and seed "
                         "configurations and write the stability table")
    main(**vars(ap.parse_args()))
