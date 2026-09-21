#!/usr/bin/env python3
"""
19_supplementary_evidence.py

Builds the supplementary evidence the reviewers asked for, all from files
already in the repository. No API calls, so it reruns offline.

  peak_manifest.csv     every ENCODE file used, its locus, term, and the
                        merged peak count and bp it contributes
  (the stability table moved to step 20 --stability, which resamples on ONE
   list of draws shared by every comparator; the version this script used to
   write drew fresh peaks per comparator and is superseded)
  ranking_provenance.json
                        the exact ranking transform, resampling unit and
                        partition rule, written down rather than described

Author: Christopher Lawrence
"""

from __future__ import annotations

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
PRIMARY = 200
B_SETTINGS = [500, 2000, 5000]
SEEDS = [20260920, 20260921, 20260922]

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
KIDNEY_TERMS = ["kidney", "glomerular visceral epithelial cell"]
COMPARISON = ["liver", "lung", "heart_left_ventricle", "stomach", "spleen",
              "thyroid_gland"]
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


def read(path, locus):
    chrom, lo, hi = locus
    iv = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c = line.split("\t")
            if c and c[0] == chrom:
                s, e = int(c[1]), int(c[2])
                if e > lo and s < hi:
                    iv.append((s, e))
    return iv


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
    """Haldane-corrected, so no resample can produce an infinite value."""
    return float(np.log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))))


def manifest() -> None:
    """Every peak file, what it is, and what it contributes at each locus."""
    term_of = {}
    for d in ("encode_celltype_peaks", "comparison_tissue_peaks"):
        for idx in (DATA / d).glob("*.index.json"):
            for f in json.loads(idx.read_text()):
                term_of[f["accession"]] = idx.name.replace(".index.json", "")
    rows = []
    for d in PEAK_DIRS:
        for p in sorted((DATA / d).glob("*.bed.gz")):
            acc = p.stem.replace(".bed", "")
            for name, locus in (("APOL1-MYH9", TEST), ("beta-globin", CTRL)):
                iv = read(p, locus)
                if not iv:
                    continue
                m = merge(iv)
                rows.append({"accession": acc, "directory": d,
                             "biosample_term": term_of.get(acc, "kidney"),
                             "locus": name, "raw_intervals": len(iv),
                             "merged_peaks": len(m),
                             "merged_bp": sum(e - s for s, e in m)})
    mf = pd.DataFrame(rows).sort_values(["locus", "directory", "accession"])
    mf.to_csv(RESULTS / "peak_manifest.csv", index=False)
    print(f"peak_manifest.csv: {mf.accession.nunique()} files, {len(mf)} rows")


def stability() -> None:
    """Recompute the six intervals at three depths of resampling, three seeds."""
    d = pd.concat([pd.read_csv(f, sep="\t") for f in
                   sorted(glob.glob(str(RESULTS / "scored_wrong_tissue" / "*.tsv.gz")))],
                  ignore_index=True)
    d = d.drop_duplicates(["pos", "af_afr"]).reset_index(drop=True)

    kid, comp = [], []
    for t in KIDNEY_TERMS:
        tag = t.replace(" ", "_")
        for f in json.loads((DATA / "comparison_tissue_peaks" /
                             f"{tag}.index.json").read_text()):
            kid += read(DATA / "comparison_tissue_peaks" /
                        f"{f['accession']}.bed.gz", TEST)
    for t in COMPARISON:
        for f in json.loads((DATA / "comparison_tissue_peaks" /
                             f"{t}.index.json").read_text()):
            comp += read(DATA / "comparison_tissue_peaks" /
                         f"{f['accession']}.bed.gz", TEST)
    kidney, other = merge(kid), merge(comp)
    k_spec = [p for p in kidney if which((p[0] + p[1]) // 2, other) is None]

    d["in_kidney"] = d["pos"].map(lambda p: which(int(p), kidney) is not None)
    d["in_kspec"] = d["pos"].map(lambda p: which(int(p), k_spec) is not None)

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

    rows = []
    for rname, col in RANKINGS.items():
        if rname == REF or col not in d.columns:
            continue
        obs = (contrast(d, RANKINGS[REF], col, "in_kspec")
               - contrast(d, RANKINGS[REF], col, "in_kidney"))
        for B in B_SETTINGS:
            for seed in SEEDS:
                rng = np.random.default_rng(seed)
                deltas = []
                for _ in range(B):
                    pick = rng.choice(elems, size=len(elems), replace=True)
                    idx = np.concatenate([by_elem[e] for e in pick] + [outside])
                    fr = d.loc[idx]
                    deltas.append(contrast(fr, RANKINGS[REF], col, "in_kspec")
                                  - contrast(fr, RANKINGS[REF], col, "in_kidney"))
                lo, hi = np.percentile(deltas, [2.5, 97.5])
                rows.append({"comparator": rname, "delta_observed": round(obs, 4),
                             "replicates": B, "seed": seed,
                             "ci_low": round(float(lo), 4),
                             "ci_high": round(float(hi), 4),
                             "excludes_zero": bool(lo > 0 or hi < 0)})
                print(f"  {rname:<18} B={B:<5} seed={seed}  "
                      f"{obs:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    st = pd.DataFrame(rows)
    st.to_csv(RESULTS / "bootstrap_stability.csv", index=False)

    sp = st.groupby("comparator").agg(
        lo_range=("ci_low", lambda x: round(x.max() - x.min(), 3)),
        hi_range=("ci_high", lambda x: round(x.max() - x.min(), 3)),
        any_excludes_zero=("excludes_zero", "any"))
    print("\nspread of interval endpoints across all settings and seeds")
    print(sp.to_string())
    return sp


def provenance() -> None:
    (RESULTS / "ranking_provenance.json").write_text(json.dumps({
        "ranking_transform": "absolute value of the raw differential "
                             "reference-versus-alternate score, sorted "
                             "descending; no quantile transform",
        "tie_handling": "pandas sort_values default (quicksort); ties are "
                        "broken arbitrarily and are not adjudicated",
        "ranked_outputs": list(RANKINGS.values()),
        "n_ranked_outputs": len(RANKINGS),
        "comparison_panel_for_partition": COMPARISON,
        "n_comparison_tissues": len(COMPARISON),
        "note": "the seven ranked outputs are not the six tissues defining "
                "the shared-peak panel; they are separate lists",
        "sequence_context_bp": 1_048_576,
        "context_placement": "variant-centred; start = max(0, pos - 524288)",
        "partition_rule": "a merged kidney peak is SHARED if a comparison peak "
                          "covers its midpoint, otherwise KIDNEY-RESTRICTED",
        "coordinate_convention": "GRCh38, BED half-open [start, end), "
                                 "0-based; intervals unclipped at locus bounds",
        "merge_rule": "touching intervals merged (next.start <= current.end)",
        "resampling_unit": "kidney-restricted merged peak",
        "resampling_fixed_set": "all variants not inside a kidney-restricted "
                                "peak, which includes variants in shared "
                                "kidney peaks",
        "multi_allelic_handling": "alternate alleles at one position share a "
                                  "position and therefore always travel "
                                  "together within one resampling unit",
        "top_set_membership": "recomputed within each resample, not held fixed",
        "repeated_peaks": "a peak drawn k times contributes its variants k "
                          "times to both the top set and the denominator",
        "zero_cells": "Haldane 0.5 correction on the bootstrap log odds "
                      "ratios, so no replicate is infinite; the reported "
                      "point estimates and exact intervals are uncorrected "
                      "conditional maximum likelihood",
        "primary_depth": PRIMARY,
        "multiplicity": "unadjusted across six comparators; exploratory",
    }, indent=2))
    print("ranking_provenance.json written")


if __name__ == "__main__":
    manifest()
    provenance()
    # stability() is superseded by 20_shared_draw_interaction.py --stability
