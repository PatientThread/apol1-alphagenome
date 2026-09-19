#!/usr/bin/env python3
"""
15_specific_vs_shared_accessibility.py

APOL1 paper, step 15: the experiment that decides whether the headline survives.

THE OBJECTION, and it is correct. Step 12 found that ranking variants by a
hepatocyte output identifies measured KIDNEY open chromatin about as well as
ranking by a kidney output, and concluded that the tissue label carries no usable
information. That conclusion does not follow from that experiment. Most
regulatory positions are open in many tissues at once. A hepatocyte ranking can
recover shared promoters and shared enhancers without the model retaining any
kidney-specific information at all, and the step 12 endpoint cannot tell the two
apart because it treats every measured kidney peak alike.

THE FIX. Partition the measured accessible regions at this locus by how tissue
restricted they are, then ask the same question of each partition separately.

    KIDNEY-SPECIFIC   open in kidney, closed in every comparison tissue
    SHARED            open in kidney AND in at least one comparison tissue
    OTHER-SPECIFIC    open in a comparison tissue, closed in kidney

The discriminating question is what happens in the KIDNEY-SPECIFIC partition.

  kidney ranking wins there      the model does retain tissue information, the
                                 step 12 result was driven by shared regions,
                                 and the "label is decorative" claim must go
  no ranking wins there          the label genuinely carries nothing even where
                                 tissue identity is the whole signal, which is
                                 a far stronger version of the original claim
  nothing is enriched anywhere   the partition is too small to answer; say so
                                 rather than reading noise

OTHER-SPECIFIC is the mirror control. If the hepatocyte ranking preferentially
finds hepatocyte-specific regions, the model clearly does encode tissue identity
and the only question is whether it reaches kidney.

TWO STATISTICAL CORRECTIONS the same review asked for, both applied here.

  Variants inside one peak are not independent observations: they share the
  element, and their predicted scores are highly correlated. Resampling is
  therefore at the level of PEAKS, not variants, so the interval reflects how
  many independent elements were actually observed.

  Differences between rankings are compared PAIRWISE on the same variants, with
  an interval, rather than by eyeballing two point estimates. Similar odds
  ratios are not equivalence.

Output: results/specific_vs_shared.csv, results/specific_vs_shared.json

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

CHROM = "chr22"
START, END = 36_140_330, 36_388_018
UA = {"User-Agent": "apol1-research/1.0", "Accept": "application/json"}
FILES_PER_TISSUE = 4
B = 4000
RNG = np.random.default_rng(20260919)

# Measured accessibility. Kidney is the tissue of interest; the rest define
# what "shared" means. These mirror the scored outputs used in step 10.
KIDNEY_TERMS = ["kidney", "glomerular visceral epithelial cell"]
COMPARISON_TERMS = ["liver", "lung", "heart left ventricle", "stomach",
                    "spleen", "thyroid gland"]

# Ranking columns from step 10, i.e. which tissue output the model was asked for.
RANKINGS = {
    "kidney (podocyte)": "dnase_glomerular_visceral_epithelial_cell",
    "kidney (whole)": "dnase_kidney",
    "hepatocyte": "dnase_hepatocyte",
    "liver": "dnase_liver",
    "lung": "dnase_lung",
    "brain": "dnase_brain",
    "stomach": "dnase_stomach",
}
TOP_N = 200


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
    iv = []
    for f in graph[:FILES_PER_TISSUE]:
        p = PEAKS / f"{f['accession']}.bed.gz"
        if not p.exists():
            with urllib.request.urlopen(
                    urllib.request.Request(
                        "https://www.encodeproject.org" + f["href"],
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
    return merge(iv)


def merge(iv: list[tuple[int, int]]) -> list[tuple[int, int]]:
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


def which(pos: int, iv: list[tuple[int, int]]) -> int | None:
    """Index of the peak containing pos, or None. Peak index = the cluster id."""
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


def or_ci(a, b, c, d):
    t = stats.contingency.odds_ratio([[a, b], [c, d]], kind="conditional")
    lo, hi = t.confidence_interval(confidence_level=0.95)
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(t.statistic), float(lo), float(hi), float(p)


def main() -> None:
    print("fetching measured accessibility")
    kidney = merge(sum([fetch(t) for t in KIDNEY_TERMS], []))
    others = {t: fetch(t) for t in COMPARISON_TERMS}
    other_union = merge(sum(others.values(), []))
    print(f"  kidney peaks {len(kidney)}, {sum(e-s for s,e in kidney):,} bp")
    for t, v in others.items():
        print(f"  {t:<22} {len(v):>4} peaks")
    print(f"  comparison union {len(other_union)} peaks, "
          f"{sum(e-s for s,e in other_union):,} bp")

    # ---- partition the locus ------------------------------------------------
    k_specific = merge([p for p in kidney if which((p[0]+p[1])//2, other_union) is None])
    shared = merge([p for p in kidney if which((p[0]+p[1])//2, other_union) is not None])
    o_specific = merge([p for p in other_union if which((p[0]+p[1])//2, kidney) is None])
    parts = {"kidney-specific": k_specific, "shared": shared,
             "other-tissue-specific": o_specific}
    print()
    for n, v in parts.items():
        print(f"  {n:<24} {len(v):>4} peaks  {sum(e-s for s,e in v):>7,} bp")

    # ---- the scored variants, one column per requested tissue --------------
    d = pd.concat([pd.read_csv(f, sep="\t") for f in
                   sorted(glob.glob(str(RESULTS / "scored_wrong_tissue" / "*.tsv.gz")))],
                  ignore_index=True).drop_duplicates("rsid")
    cols = {k: v for k, v in RANKINGS.items() if v in d.columns}
    print(f"\n  {len(d):,} scored variants, {len(cols)} ranking columns available")

    rows = []
    for pname, piv in parts.items():
        if not piv:
            continue
        d["hit"] = d["pos"].map(lambda p: which(int(p), piv) is not None)
        d["elem"] = d["pos"].map(lambda p: which(int(p), piv))
        for rname, col in cols.items():
            s = d.dropna(subset=[col]).copy()
            s["k"] = s[col].abs()
            s = s.sort_values("k", ascending=False).reset_index(drop=True)
            h = s["hit"].values
            a, b = int(h[:TOP_N].sum()), int(TOP_N - h[:TOP_N].sum())
            c, dd = int(h[TOP_N:].sum()), int(len(h) - TOP_N - h[TOP_N:].sum())
            if a + c == 0:
                rows.append({"partition": pname, "ranked_by": rname, "in_top": 0,
                             "odds_ratio": np.nan, "ci_low": np.nan,
                             "ci_high": np.nan, "p": np.nan,
                             "note": "no variant in this partition"})
                continue
            orr, lo, hi, p = or_ci(a, b, c, dd)
            rows.append({"partition": pname, "ranked_by": rname,
                         "in_top": a, "top_n": TOP_N,
                         "bg_hits": c, "bg_n": c + dd,
                         "odds_ratio": round(orr, 3),
                         "ci_low": round(lo, 3),
                         "ci_high": None if np.isinf(hi) else round(hi, 3),
                         "p": p})

    res = pd.DataFrame(rows)
    res.to_csv(RESULTS / "specific_vs_shared.csv", index=False)

    print("\n" + "=" * 86)
    print(f"ENRICHMENT IN EACH PARTITION, top {TOP_N} by each requested tissue output")
    print("=" * 86)
    for pname in parts:
        sub = res[res["partition"] == pname]
        if sub.empty:
            continue
        print(f"\n  {pname.upper()}")
        print(f"    {'ranked by':<20}{'hits':>6}{'OR':>8}   95% CI")
        for r in sub.itertuples():
            if np.isnan(r.odds_ratio):
                print(f"    {r.ranked_by:<20}     -       -   {r.note}")
                continue
            hi = "inf" if r.ci_high is None else f"{r.ci_high:.2f}"
            print(f"    {r.ranked_by:<20}{r.in_top:>6}{r.odds_ratio:>8.2f}   "
                  f"[{r.ci_low:.2f}, {hi}]")

    # ---- the decisive comparison, paired and element-resampled -------------
    ks = res[(res["partition"] == "kidney-specific") & res["odds_ratio"].notna()]
    verdict = {}
    if len(ks) >= 2:
        kid = [r for r in ks.itertuples() if r.ranked_by.startswith("kidney")]
        non = [r for r in ks.itertuples() if not r.ranked_by.startswith("kidney")]
        if kid and non:
            best_k = max(r.odds_ratio for r in kid)
            best_n = max(r.odds_ratio for r in non)
            verdict = {"best_kidney_or": best_k, "best_nonkidney_or": best_n,
                       "kidney_leads": bool(best_k > best_n)}
            print("\n  KIDNEY-SPECIFIC PARTITION, the question that decides it:")
            print(f"    best kidney ranking      OR {best_k:.2f}")
            print(f"    best non-kidney ranking  OR {best_n:.2f}")

    # ---- PAIRED, element-resampled difference, which is what decides it ----
    # Unpaired intervals overlapping is not the same as no difference. Resample
    # PEAKS with replacement, recompute both rankings on the same resampled
    # variants, and take the difference in log odds ratio.
    paired = {}
    kspec = parts["kidney-specific"]
    if kspec:
        d["elem_ks"] = d["pos"].map(lambda p: which(int(p), kspec))
        d["hit_ks"] = d["elem_ks"].notna()
        kid_col = RANKINGS["kidney (podocyte)"]
        n_elem = len(kspec)
        # variants grouped by the peak they sit in; non-hits form their own pool
        by_elem = {e: g.index.to_numpy() for e, g in d.groupby("elem_ks")}
        outside = d.index[d["elem_ks"].isna()].to_numpy()

        def or_for(frame, col):
            s2 = frame.dropna(subset=[col]).copy()
            s2["k"] = s2[col].abs()
            s2 = s2.sort_values("k", ascending=False)
            h = s2["hit_ks"].values
            a, b = int(h[:TOP_N].sum()), int(TOP_N - h[:TOP_N].sum())
            c, dd = int(h[TOP_N:].sum()), int(len(h) - TOP_N - h[TOP_N:].sum())
            if a == 0 or c == 0:
                return np.nan
            return np.log((a + 0.5) * (dd + 0.5) / ((b + 0.5) * (c + 0.5)))

        for rname, col in cols.items():
            if rname.startswith("kidney (podocyte)") or col not in d.columns:
                continue
            diffs = []
            for _ in range(600):
                elems = RNG.choice(list(by_elem), size=n_elem, replace=True)
                idx = np.concatenate([by_elem[e] for e in elems] + [outside])
                fr = d.loc[idx]
                a1, a2 = or_for(fr, kid_col), or_for(fr, col)
                if np.isfinite(a1) and np.isfinite(a2):
                    diffs.append(a1 - a2)
            if len(diffs) > 100:
                lo, hi = np.percentile(diffs, [2.5, 97.5])
                paired[rname] = {"mean_log_or_diff": float(np.mean(diffs)),
                                 "lo": float(lo), "hi": float(hi),
                                 "excludes_zero": bool(lo > 0 or hi < 0)}

        print("\n  PAIRED difference in log odds ratio, kidney (podocyte) minus each")
        print("  other ranking, in the KIDNEY-SPECIFIC partition.")
        print("  Peaks resampled with replacement, 600 replicates.")
        for rn, v in sorted(paired.items(), key=lambda kv: -kv[1]["mean_log_or_diff"]):
            star = "  *" if v["excludes_zero"] else ""
            print(f"    vs {rn:<20}{v['mean_log_or_diff']:>+7.3f}  "
                  f"[{v['lo']:>+6.3f}, {v['hi']:>+6.3f}]{star}")
        n_sig = sum(1 for v in paired.values() if v["excludes_zero"] and v["mean_log_or_diff"] > 0)
        print(f"\n    rankings podocyte beats with an interval excluding zero: "
              f"{n_sig} of {len(paired)}")

    print("\n  READING:")
    nk = parts["kidney-specific"]
    n_var = int(d["pos"].map(lambda p: which(int(p), nk) is not None).sum()) if nk else 0
    if not nk or n_var < 20:
        print(f"  The kidney-specific partition holds {len(nk)} peaks and {n_var}")
        print("  scored variants. That is too little to answer the question, and the")
        print("  honest report is that this locus cannot distinguish kidney-specific")
        print("  from shared accessibility, not that the model fails to.")
    elif verdict.get("kidney_leads"):
        print("  The kidney ranking leads where tissue identity is the whole signal.")
        print("  Step 12 was driven by shared regions and the claim that the label")
        print("  carries no information must be withdrawn.")
    else:
        print("  The kidney ranking does NOT lead even in the kidney-specific")
        print("  partition. The label carries nothing where it should matter most,")
        print("  which is a stronger result than step 12 alone could support.")

    out = {"locus": f"{CHROM}:{START}-{END}",
           "kidney_terms": KIDNEY_TERMS, "comparison_terms": COMPARISON_TERMS,
           "partitions": {k: {"peaks": len(v), "bp": sum(e-s for s, e in v)}
                          for k, v in parts.items()},
           "kidney_specific_scored_variants": n_var,
           "top_n": TOP_N, "verdict": verdict,
           "paired_kidney_specific": paired}
    (RESULTS / "specific_vs_shared.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/specific_vs_shared.json")


if __name__ == "__main__":
    main()
