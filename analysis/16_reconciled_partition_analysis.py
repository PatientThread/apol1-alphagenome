#!/usr/bin/env python3
"""
16_reconciled_partition_analysis.py

APOL1 paper, step 16: one frozen universe, and a direct test of the interaction.

WHAT WENT WRONG, AND IT IS FATAL TO THE DRAFT-3 HEADLINE AS WRITTEN.

Review observed that the pooled and partitioned kidney sets cannot be the same
data: pooled gave 30,358 bp and 69 podocyte hits in the top 200, while the
partitions summed to 44,576 bp and 101 hits. A disjoint partition must sum to its
parent. It did not, because the two analyses used DIFFERENT measured peak sets.
Step 5 queried ENCODE for biosample term "kidney" alone. Step 15 used "kidney"
AND "glomerular visceral epithelial cell". So draft 3 compared a pooled result
built on whole-kidney peaks against partitions built on whole-kidney plus
podocyte peaks, and called the difference an effect of partitioning.

Two further errors of the same family. The scored files hold 1,930 rows but only
1,921 unique rsIDs, so every denominator quoted as 1,930 was wrong. And the
claim that pooling hides tissue discrimination rested on ranking depth 100, while
the declared primary depth is 200 where podocyte already led the pooled table.

THIS SCRIPT THEREFORE DOES THREE THINGS.

  1. FREEZES ONE UNIVERSE. One variant set, deduplicated by rsID. One kidney peak
     union. Partitions defined so that kidney-specific and shared are mutually
     exclusive and their union IS the pooled kidney set, by construction rather
     than by hope. Every count is emitted so the arithmetic can be checked.

  2. COMPUTES POOLED AND PARTITIONED FROM THAT SAME OBJECT, at every ranking
     depth, with all four contingency cells reported.

  3. TESTS THE INTERACTION DIRECTLY, which draft 3 never did. The question is not
     whether podocyte beats a comparator in the kidney-specific partition, nor
     whether it does so pooled. It is whether the CONTRAST CHANGES:

         delta = [logOR_podocyte - logOR_comparator] in kidney-specific regions
               - [logOR_podocyte - logOR_comparator] pooled

     estimated on the SAME resampled peaks, so the two contrasts are paired and
     their difference has an interval. Separate significance statements in two
     analyses do not establish that partitioning changed anything.

If delta does not exclude zero, the draft-3 headline is not supported and must be
replaced by the weaker, true statement: the partitions differ in what they show,
but partitioning has not been shown to change tissue discrimination.

Output: results/reconciled_counts.csv      every 2x2, every depth, every ranking
        results/reconciled_interaction.csv the paired interaction estimates
        results/reconciled.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import glob
import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
PEAKS = DATA / "comparison_tissue_peaks"
PEAKS.mkdir(parents=True, exist_ok=True)

CHROM, START, END = "chr22", 36_140_330, 36_388_018
UA = {"User-Agent": "apol1-research/1.0", "Accept": "application/json"}
FILES_PER_TISSUE = 4
B = 2000
RNG = np.random.default_rng(20260920)

# ONE definition of the kidney measured set, used for pooled AND partitioned.
KIDNEY_TERMS = ["kidney", "glomerular visceral epithelial cell"]
COMPARISON_TERMS = ["liver", "lung", "heart left ventricle", "stomach",
                    "spleen", "thyroid gland"]
RANKINGS = {
    "kidney (podocyte)": "dnase_glomerular_visceral_epithelial_cell",
    "kidney (whole)": "dnase_kidney",
    "hepatocyte": "dnase_hepatocyte",
    "liver": "dnase_liver",
    "lung": "dnase_lung",
    "brain": "dnase_brain",
    "stomach": "dnase_stomach",
}
REF = "kidney (podocyte)"
DEPTHS = [25, 50, 100, 200]
PRIMARY = 200


def fetch(term: str) -> list[tuple[int, int]]:
    tag = term.replace(" ", "_")
    idx = PEAKS / f"{tag}.index.json"
    if idx.exists():
        graph = json.loads(idx.read_text())
    else:
        url = ("https://www.encodeproject.org/search/?type=File&file_format=bed"
               "&output_type=peaks&assay_title=DNase-seq"
               f"&biosample_ontology.term_name={urllib.parse.quote(term)}"
               "&assembly=GRCh38&status=released"
               f"&limit={FILES_PER_TISSUE}&format=json&field=accession&field=href")
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                    timeout=180) as fh:
            graph = json.load(fh)["@graph"]
        idx.write_text(json.dumps(graph))
        time.sleep(0.4)
    iv, accs = [], []
    for f in graph[:FILES_PER_TISSUE]:
        accs.append(f["accession"])
        p = PEAKS / f"{f['accession']}.bed.gz"
        if not p.exists():
            with urllib.request.urlopen(
                    urllib.request.Request("https://www.encodeproject.org" + f["href"],
                                           headers=UA), timeout=300) as fh, open(p, "wb") as o:
                o.write(fh.read())
            time.sleep(0.4)
        with gzip.open(p, "rt") as fh:
            for line in fh:
                c = line.split("\t")
                if c and c[0] == CHROM:
                    s, e = int(c[1]), int(c[2])
                    if e > START and s < END:
                        iv.append((s, e))
    return merge(iv), accs


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
    """Haldane-corrected log odds ratio, finite for zero cells."""
    return float(np.log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))))


def main() -> None:
    kid_iv, kid_acc = [], []
    for t in KIDNEY_TERMS:
        iv, acc = fetch(t)
        kid_iv += iv
        kid_acc += acc
    kidney = merge(kid_iv)
    comp_iv, comp_acc = [], {}
    for t in COMPARISON_TERMS:
        iv, acc = fetch(t)
        comp_iv += iv
        comp_acc[t] = acc
    other = merge(comp_iv)

    # Partitions by construction: every kidney peak is in exactly one of them.
    k_spec = [p for p in kidney if which((p[0] + p[1]) // 2, other) is None]
    shared = [p for p in kidney if which((p[0] + p[1]) // 2, other) is not None]
    assert len(k_spec) + len(shared) == len(kidney)
    bp = lambda v: sum(e - s for s, e in v)
    print("FROZEN MEASURED SET")
    print(f"  kidney union        {len(kidney):>4} peaks  {bp(kidney):>7,} bp")
    print(f"    kidney-specific   {len(k_spec):>4} peaks  {bp(k_spec):>7,} bp")
    print(f"    shared            {len(shared):>4} peaks  {bp(shared):>7,} bp")
    print(f"    sum checks        {len(k_spec)+len(shared)} peaks  "
          f"{bp(k_spec)+bp(shared):,} bp  "
          f"{'OK' if bp(k_spec)+bp(shared) == bp(kidney) else 'MISMATCH'}")

    # ---- one variant universe, deduplicated ------------------------------
    d = pd.concat([pd.read_csv(f, sep="\t") for f in
                   sorted(glob.glob(str(RESULTS / "scored_wrong_tissue" / "*.tsv.gz")))],
                  ignore_index=True)
    before = len(d)
    # Deduplicate on the ALLELE, not the identifier. Nine rsIDs appear twice at
    # the same position with different African allele frequencies: these are
    # multi-allelic sites where one identifier covers two alternate alleles, and
    # they are distinct variants. An earlier version deduplicated on rsID and
    # silently discarded nine real variants. The scored file carries no ref/alt
    # column, so (position, allele frequency) serves as the allele key here; the
    # full ledger with ref and alt is deposited.
    d = d.drop_duplicates(["pos", "af_afr"]).reset_index(drop=True)
    print(f"\nFROZEN VARIANT UNIVERSE: {len(d)} distinct variants "
          f"at {d['pos'].nunique()} positions "
          f"({before - len(d)} exact duplicate rows removed)")

    d["in_kidney"] = d["pos"].map(lambda p: which(int(p), kidney) is not None)
    d["in_kspec"] = d["pos"].map(lambda p: which(int(p), k_spec) is not None)
    d["in_shared"] = d["pos"].map(lambda p: which(int(p), shared) is not None)
    assert (d["in_kspec"] & d["in_shared"]).sum() == 0
    assert (d["in_kspec"] | d["in_shared"]).equals(d["in_kidney"])
    print(f"  variants in kidney union {int(d.in_kidney.sum())} = "
          f"{int(d.in_kspec.sum())} specific + {int(d.in_shared.sum())} shared  OK")

    cols = {k: v for k, v in RANKINGS.items() if v in d.columns}
    endpoints = {"pooled": "in_kidney", "kidney-specific": "in_kspec",
                 "shared": "in_shared"}

    rows = []
    for ep, flag in endpoints.items():
        for rname, col in cols.items():
            s = d.dropna(subset=[col]).copy()
            s["k"] = s[col].abs()
            s = s.sort_values("k", ascending=False).reset_index(drop=True)
            h = s[flag].values
            for n in DEPTHS:
                a, b = int(h[:n].sum()), int(n - h[:n].sum())
                c, dd = int(h[n:].sum()), int(len(h) - n - h[n:].sum())
                t = stats.contingency.odds_ratio([[a, b], [c, dd]], kind="conditional")
                lo, hi = t.confidence_interval(confidence_level=0.95)
                rows.append({"endpoint": ep, "ranked_by": rname, "depth": n,
                             "a_top_in": a, "b_top_out": b,
                             "c_bg_in": c, "d_bg_out": dd, "eligible_n": len(s),
                             "odds_ratio": round(float(t.statistic), 3),
                             "ci_low": round(float(lo), 3),
                             "ci_high": None if np.isinf(hi) else round(float(hi), 3)})
    counts = pd.DataFrame(rows)
    counts.to_csv(RESULTS / "reconciled_counts.csv", index=False)

    print(f"\nPOOLED vs PARTITIONED at the primary depth of {PRIMARY}, same data")
    print(f"  {'ranked by':<20}{'pooled':>9}{'k-specific':>12}{'shared':>9}")
    for rname in cols:
        g = counts[(counts.ranked_by == rname) & (counts.depth == PRIMARY)]
        v = {r.endpoint: r.odds_ratio for r in g.itertuples()}
        print(f"  {rname:<20}{v['pooled']:>9.2f}{v['kidney-specific']:>12.2f}"
              f"{v['shared']:>9.2f}")

    # ---- THE INTERACTION, paired on shared resamples ----------------------
    ks_elem = d["pos"].map(lambda p: which(int(p), k_spec))
    by_elem = {e: g.index.to_numpy() for e, g in d.groupby(ks_elem)}
    outside = d.index[ks_elem.isna()].to_numpy()
    elems = list(by_elem)

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

    inter = []
    for rname, col in cols.items():
        if rname == REF:
            continue
        obs_ks = contrast(d, cols[REF], col, "in_kspec")
        obs_pool = contrast(d, cols[REF], col, "in_kidney")
        deltas = []
        for _ in range(B // 4):
            pick = RNG.choice(elems, size=len(elems), replace=True)
            idx = np.concatenate([by_elem[e] for e in pick] + [outside])
            fr = d.loc[idx]
            deltas.append(contrast(fr, cols[REF], col, "in_kspec")
                          - contrast(fr, cols[REF], col, "in_kidney"))
        lo, hi = np.percentile(deltas, [2.5, 97.5])
        inter.append({"comparator": rname,
                      "contrast_kidney_specific": round(obs_ks, 3),
                      "contrast_pooled": round(obs_pool, 3),
                      "delta_observed": round(obs_ks - obs_pool, 3),
                      "delta_ci_low": round(float(lo), 3),
                      "delta_ci_high": round(float(hi), 3),
                      "excludes_zero": bool(lo > 0 or hi < 0)})
    it = pd.DataFrame(inter).sort_values("delta_observed", ascending=False)
    it.to_csv(RESULTS / "reconciled_interaction.csv", index=False)

    print("\n" + "=" * 84)
    print("DOES PARTITIONING CHANGE THE TISSUE CONTRAST? (paired, same resamples)")
    print("=" * 84)
    print(f"  {'vs':<20}{'k-spec':>9}{'pooled':>9}{'delta':>9}   95% CI of delta")
    for r in it.itertuples():
        star = "  *" if r.excludes_zero else ""
        print(f"  {r.comparator:<20}{r.contrast_kidney_specific:>+9.2f}"
              f"{r.contrast_pooled:>+9.2f}{r.delta_observed:>+9.2f}   "
              f"[{r.delta_ci_low:>+6.2f}, {r.delta_ci_high:>+6.2f}]{star}")
    n_sig = int(it["excludes_zero"].sum())
    print(f"\n  comparators whose contrast CHANGED with partitioning: "
          f"{n_sig} of {len(it)}")

    print("\n  READING:")
    if n_sig == 0:
        print("  Partitioning did NOT measurably change the tissue contrast. The")
        print("  draft-3 headline is unsupported. What can be said is that the")
        print("  partitions differ in absolute enrichment, not that partitioning")
        print("  reveals discrimination the pooled endpoint hides.")
    else:
        print(f"  The contrast changed for {n_sig} of {len(it)} comparators, so")
        print("  partitioning does alter what the benchmark shows. State the")
        print("  magnitude and direction per comparator rather than generally.")

    out = {"kidney_accessions": kid_acc, "comparison_accessions": comp_acc,
           "kidney_peaks": len(kidney), "kidney_bp": bp(kidney),
           "kidney_specific_peaks": len(k_spec), "kidney_specific_bp": bp(k_spec),
           "shared_peaks": len(shared), "shared_bp": bp(shared),
           "partition_sums_to_parent": bp(k_spec) + bp(shared) == bp(kidney),
           "variant_universe": int(len(d)),
           "duplicates_removed": int(before - len(d)),
           "primary_depth": PRIMARY, "bootstrap_replicates": B // 4,
           "resampling_unit": "kidney-specific peaks",
           "n_comparators_contrast_changed": n_sig,
           "n_comparators": int(len(it))}
    (RESULTS / "reconciled.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote reconciled_counts.csv, reconciled_interaction.csv, reconciled.json")


if __name__ == "__main__":
    main()
