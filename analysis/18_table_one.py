#!/usr/bin/env python3
"""
18_table_one.py

Regenerates every row of Table 1 under ONE merge convention and ONE
deduplication rule, so the six rows are comparable with each other.

The manuscript's point is that the conclusion changed with the measured
reference set. That point is only legible if the rows differ ONLY in the
reference set. Earlier scripts were written at different times and differ in
incidentals (merge on `<` versus `<=`, dedup on rsID versus on the allele),
which moves peak counts by one or two and n by nine. Those differences are
immaterial to any conclusion but they make the table look unreliable. This
script recomputes all of it in one pass.

Reference sets compared, all DNase-seq, GRCh38, released, from the ENCODE
portal, all merged the same way:

  kidney_only     biosample_ontology.term_name=kidney
  kidney_union    that plus term_name=glomerular visceral epithelial cell
  kspec / shared  kidney_union partitioned on six comparison tissues

Output: results/table_one.csv, results/table_one.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

ROOT = Path(__file__).resolve().parent.parent
DATA, RESULTS = ROOT / "data", ROOT / "results"

TEST = ("chr22", 36_140_330, 36_388_018)
CTRL = ("chr11", 5_150_000, 5_397_688)
DEPTH = 200
PODO = "dnase_glomerular_visceral_epithelial_cell"

# term_name=kidney, the set steps 5, 10 and 12 used.
KIDNEY_ONLY = ["ENCFF473YKR", "ENCFF554CKH", "ENCFF350FBW", "ENCFF214HJG",
               "ENCFF851SWR", "ENCFF922WKR", "ENCFF948WNJ", "ENCFF146AYH"]
# term_name=kidney + glomerular visceral epithelial cell, the reconciled set.
KIDNEY_UNION = ["ENCFF473YKR", "ENCFF554CKH", "ENCFF350FBW", "ENCFF214HJG",
                "ENCFF618RWK", "ENCFF325DZQ", "ENCFF019PKU", "ENCFF971ABV"]
COMPARISON = ["liver", "lung", "heart_left_ventricle", "stomach", "spleen",
              "thyroid_gland"]

PEAK_DIRS = [DATA / "encode_peaks", DATA / "encode_celltype_peaks",
             DATA / "comparison_tissue_peaks", DATA / "control_locus_peaks",
             DATA / "control_locus_peaks_matched"]


def find(acc: str) -> Path:
    for d in PEAK_DIRS:
        p = d / f"{acc}.bed.gz"
        if p.exists():
            return p
    raise FileNotFoundError(acc)


def merge(iv: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """One convention everywhere, and it is step 16's: touching intervals
    merge (s <= end), intervals unclipped."""
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


def load(accs: list[str], locus) -> list[tuple[int, int]]:
    chrom, lo, hi = locus
    iv = []
    for a in accs:
        with gzip.open(find(a), "rt") as fh:
            for line in fh:
                f = line.split("\t")
                if f and f[0] == chrom:
                    s, e = int(f[1]), int(f[2])
                    if e > lo and s < hi:
                        # NOT clipped to the locus bounds. Step 16 merges the
                        # unclipped intervals and the partitions and the
                        # bootstrap are defined on those peaks; clipping here
                        # would silently redefine them and reintroduce exactly
                        # the mismatch this paper is about.
                        iv.append((s, e))
    return merge(iv)


def hit(pos: np.ndarray, iv) -> np.ndarray:
    if not iv:
        return np.zeros(len(pos), bool)
    starts = np.array([s for s, _ in iv])
    ends = np.array([e for _, e in iv])
    j = np.searchsorted(starts, pos, "right") - 1
    ok = j >= 0
    out = np.zeros(len(pos), bool)
    out[ok] = pos[ok] < ends[j[ok]]
    return out


def enrich(d: pd.DataFrame, col: str, inpeak: np.ndarray, depth: int) -> dict:
    order = d[col].abs().values.argsort()[::-1]
    top = np.zeros(len(d), bool)
    top[order[:depth]] = True
    a = int((top & inpeak).sum())
    b = int((top & ~inpeak).sum())
    c = int((~top & inpeak).sum())
    dd = int((~top & ~inpeak).sum())
    t = [[a, b], [c, dd]]
    r = odds_ratio(t, kind="conditional")
    lo, hi = r.confidence_interval(0.95)
    return {"a_top_in": a, "b_top_out": b, "c_bg_in": c, "d_bg_out": dd,
            "eligible_n": a + b + c + dd,
            "odds_ratio": round(float(r.statistic), 3),
            "ci_low": round(float(lo), 3), "ci_high": round(float(hi), 3),
            "p_one_sided": float(fisher_exact(t, "greater")[1])}


def variants(scored_dir: Path) -> pd.DataFrame:
    d = pd.concat([pd.read_csv(p, sep="\t", low_memory=False)
                   for p in sorted(scored_dir.glob("chunk_*.tsv.gz"))],
                  ignore_index=True)
    # Deduplicate on the ALLELE where ref/alt are carried, otherwise on
    # position and frequency, which separates the nine multi-allelic sites.
    keys = ["pos", "ref", "alt"] if "ref" in d.columns else ["pos", "af_afr"]
    return d.drop_duplicates(keys).sort_values("pos").reset_index(drop=True)


def main() -> None:
    test = variants(RESULTS / "scored")
    wrong = variants(RESULTS / "scored_wrong_tissue")
    ctrl = variants(RESULTS / "scored_control_locus")

    k_only = load(KIDNEY_ONLY, TEST)
    k_union = load(KIDNEY_UNION, TEST)
    comp = load([a for t in COMPARISON
                 for a in [x["accession"] for x in json.loads(
                     (DATA / "comparison_tissue_peaks" /
                      f"{t}.index.json").read_text())]
                 ], TEST) if all((DATA / "comparison_tissue_peaks" /
                                  f"{t}.index.json").exists()
                                 for t in COMPARISON) else None
    if comp is None:
        raise SystemExit("comparison tissue index files missing")

    # Step 16's rule, kept verbatim: a kidney peak is SHARED if its midpoint
    # lies inside a comparison peak, otherwise kidney-restricted. The
    # partitions and the Figure 2 bootstrap are defined on this rule, so
    # substituting any-overlap here would put the table and the figure on
    # different peak sets, which is the error this paper reports.
    mid = lambda p: (p[0] + p[1]) // 2
    inside = lambda x: hit(np.array([x], dtype=np.int64), comp)[0]
    kspec = [p for p in k_union if not inside(mid(p))]
    shared = [p for p in k_union if inside(mid(p))]
    assert len(kspec) + len(shared) == len(k_union)
    assert (sum(e - s for s, e in kspec) + sum(e - s for s, e in shared)
            == sum(e - s for s, e in k_union))

    c_only = load(KIDNEY_ONLY, CTRL)
    c_union = load(KIDNEY_UNION, CTRL)

    span_t = TEST[2] - TEST[1]
    span_c = CTRL[2] - CTRL[1]
    rows = []

    def add(label, locus, d, iv, span, best_nonrenal=""):
        pos = d["pos"].values.astype(np.int64)
        r = enrich(d, PODO, hit(pos, iv), DEPTH)
        bp = sum(e - s for s, e in iv)
        rows.append({"reference_set": label, "locus": locus,
                     "peaks": len(iv), "bp": bp,
                     "pct_of_locus": round(100 * bp / span, 1),
                     "best_nonrenal": best_nonrenal, **r})

    # Best non-renal ranking against kidney-only peaks, the comparison that
    # produced the first conclusion. Uses the wrong-tissue scoring run, which
    # is the only run carrying non-renal columns.
    wpos = wrong["pos"].values.astype(np.int64)
    w_in = hit(wpos, k_only)
    nonrenal = [c for c in wrong.columns
                if c.startswith("dnase_") and c not in
                ("dnase_kidney", PODO)]
    best = max(((c, enrich(wrong, c, w_in, DEPTH)) for c in nonrenal),
               key=lambda kv: kv[1]["odds_ratio"])
    best_at_100 = max(((c, enrich(wrong, c, w_in, 100)) for c in nonrenal),
                      key=lambda kv: kv[1]["odds_ratio"])

    add("whole kidney only", "APOL1-MYH9", wrong, k_only, span_t,
        f"{best[1]['odds_ratio']:.2f} ({best[0][6:]})")
    # All four APOL1-MYH9 rows use the step 16 variant universe, so the rows
    # differ only in the measured reference set.
    add("kidney + podocyte", "APOL1-MYH9", wrong, k_union, span_t)
    add("kidney + podocyte, kidney-restricted", "APOL1-MYH9", wrong, kspec,
        span_t)
    add("kidney + podocyte, shared", "APOL1-MYH9", wrong, shared, span_t)
    add("whole kidney only", "beta-globin", ctrl, c_only, span_c)
    add("kidney + podocyte", "beta-globin", ctrl, c_union, span_c)

    # Every ranking against the kidney-only reference, for Figure 1's first
    # panel. This is the comparison that produced the first conclusion, so it
    # belongs in the figure and not only in the table.
    NAME = {"dnase_glomerular_visceral_epithelial_cell": "kidney (podocyte)",
            "dnase_kidney": "kidney (whole)", "dnase_lung": "lung",
            "dnase_hepatocyte": "hepatocyte", "dnase_stomach": "stomach",
            "dnase_brain": "brain", "dnase_liver": "liver"}
    kidney_only_rows = [{"ranked_by": NAME[c], "endpoint": "kidney-only",
                         "depth": DEPTH, **enrich(wrong, c, w_in, DEPTH)}
                        for c in NAME]
    pd.DataFrame(kidney_only_rows).to_csv(
        RESULTS / "kidney_only_counts.csv", index=False)
    print("\nkidney-only reference, every ranking, depth 200")
    print(pd.DataFrame(kidney_only_rows)[
        ["ranked_by", "a_top_in", "odds_ratio", "ci_low", "ci_high"]
    ].to_string(index=False))

    t1 = pd.DataFrame(rows)
    t1.to_csv(RESULTS / "table_one.csv", index=False)
    print(t1.to_string(index=False))
    print(f"\nbest non-renal at depth 100: {best_at_100[0][6:]} "
          f"OR {best_at_100[1]['odds_ratio']}")
    print(f"best non-renal at depth 200: {best[0][6:]} "
          f"OR {best[1]['odds_ratio']}")
    print(f"\npodocyte vs kidney-only peaks, depth 100: "
          f"{enrich(wrong, PODO, w_in, 100)['odds_ratio']}")
    (RESULTS / "table_one.json").write_text(json.dumps(
        {"rows": rows, "depth": DEPTH,
         "best_nonrenal_depth100": {"track": best_at_100[0],
                                    **best_at_100[1]},
         "best_nonrenal_depth200": {"track": best[0], **best[1]},
         "merge_convention": "touching intervals merged (s <= end)",
         "dedup": "pos+ref+alt where carried, else pos+af_afr"}, indent=2))


if __name__ == "__main__":
    main()
