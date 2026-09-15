#!/usr/bin/env python3
"""
08_annotate_candidates.py

APOL1 paper, step 8: say what the candidates actually are.

THE PROBLEM THIS SOLVES. After steps 4 to 7 the candidates are positions,
frequencies and linkage values. That is enough to defend them statistically and
not nearly enough for a reader to judge whether any of them is biologically
plausible. "A variant 90 kb from APOL1 at 40% frequency in African populations"
tells a clinician nothing. "A variant inside an annotated enhancer that is open
in podocytes but not in whole kidney" tells them something they can act on.

WHAT IS ADDED, all from sources independent of the model.

  1. WHICH KIDNEY CELL TYPES. Each candidate is intersected with the ENCODE peak
     files one at a time rather than the merged set used in step 5, and each file
     is resolved to its biosample and assay through the ENCODE portal. A variant
     open in podocytes only is a different proposition from one open in every
     kidney sample, and the merged analysis in step 5 deliberately could not
     distinguish them.

  2. WHICH REGULATORY ELEMENT, IF ANY. The Ensembl Regulatory Build is queried
     for overlapping annotated features: promoter, enhancer, CTCF binding site,
     open chromatin region. This is an annotation built from many cell types and
     is not kidney-specific, so it answers "is this a regulatory element at all"
     rather than "is this a kidney regulatory element". Both questions matter and
     they are reported separately.

  3. WHAT ELSE IS NEARBY. The nearest gene and the distance to it, because a
     candidate 90 kb from APOL1 may be much closer to something else, and if so
     the assumption that APOL1 is its target needs stating rather than assuming.

WHAT THIS CANNOT DO. None of this is functional evidence. An annotated enhancer
is a prediction from chromatin marks, not a demonstration that the element
regulates APOL1 or that the variant changes its activity. The paper must not
slide from "sits in an annotated enhancer" to "disrupts an enhancer".

Output: results/candidate_annotation.csv
        results/candidate_annotation.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import gzip
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
PEAKS = ROOT / "data" / "encode_peaks"
CACHE = RESULTS / "annot_cache"
CACHE.mkdir(parents=True, exist_ok=True)

ENCODE = "https://www.encodeproject.org"
ENSEMBL = "https://rest.ensembl.org"
INDEPENDENT_R2 = 0.2          # the threshold used in step 7


def get_json(url: str, headers: dict | None = None, tag: str = "") -> dict | list:
    key = CACHE / (tag + ".json")
    if key.exists():
        return json.loads(key.read_text())
    req = urllib.request.Request(url, headers=headers or {"Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                d = json.loads(r.read().decode())
            key.write_text(json.dumps(d))
            time.sleep(0.34)
            return d
        except Exception:                                   # noqa: BLE001
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"failed repeatedly: {url}")


def peak_file_metadata() -> dict[str, dict]:
    """Resolve each ENCODE peak file to its biosample and assay."""
    meta = {}
    for f in sorted(PEAKS.glob("*.bed.gz")):
        acc = f.stem.split(".")[0]
        d = get_json(f"{ENCODE}/files/{acc}/?format=json", tag=f"encode_{acc}")
        bio = d.get("biosample_ontology") or {}
        meta[acc] = {
            "accession": acc,
            "assay": d.get("assay_term_name"),
            "biosample": bio.get("term_name"),
            "classification": bio.get("classification"),
            "assembly": d.get("assembly"),
            "path": f,
        }
    return meta


def load_peaks(path: Path, chrom: str = "chr22") -> list[tuple[int, int]]:
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.split("\t")
            if p and p[0] == chrom:
                out.append((int(p[1]), int(p[2])))
    return sorted(out)


def main() -> None:
    ld = pd.read_csv(RESULTS / "ld_gene_haplotype_per_variant.csv")
    ld["independent"] = ld["max_r2_vs_gene"] < INDEPENDENT_R2
    print(f"candidates: {len(ld)}  "
          f"({int(ld['independent'].sum())} independent of gene haplotype)\n")

    meta = peak_file_metadata()
    print("ENCODE peak files resolved to biosample:")
    for m in meta.values():
        print(f"  {m['accession']}  {m['assay']:<10} {m['classification']:<12} "
              f"{m['biosample']}")
    peaks = {a: load_peaks(m["path"]) for a, m in meta.items()}
    print()

    rows = []
    for c in ld.itertuples():
        pos = int(c.pos)

        # ---- which kidney biosamples is it open in ------------------------
        open_in = []
        for acc, m in meta.items():
            if any(s <= pos < e for s, e in peaks[acc]):
                open_in.append(f"{m['biosample']} ({m['assay']})")

        # ---- annotated regulatory element ---------------------------------
        reg = get_json(
            f"{ENSEMBL}/overlap/region/human/22:{pos-1}-{pos+1}"
            f"?feature=regulatory;content-type=application/json",
            tag=f"ens_reg_{pos}")
        feats = sorted({(r.get("description") or r.get("feature_type"))
                        for r in reg}) if reg else []

        # ---- nearest gene, which need not be APOL1 ------------------------
        genes = get_json(
            f"{ENSEMBL}/overlap/region/human/22:{pos-100000}-{pos+100000}"
            f"?feature=gene;biotype=protein_coding;content-type=application/json",
            tag=f"ens_gene_{pos}")
        near = None
        if genes:
            def dist(g):
                if g["start"] <= pos <= g["end"]:
                    return 0
                return min(abs(pos - g["start"]), abs(pos - g["end"]))
            g = min(genes, key=dist)
            near = (g.get("external_name") or g.get("id"), dist(g))

        rows.append({
            "rsid": c.rsid,
            "pos": pos,
            "af_afr": c.af_afr,
            "max_r2_vs_gene": c.max_r2_vs_gene,
            "independent": bool(c.independent),
            "n_kidney_samples_open": len(open_in),
            "open_in": "; ".join(open_in) if open_in else "none",
            "regulatory_features": "; ".join(feats) if feats else "none annotated",
            "nearest_protein_coding_gene": near[0] if near else None,
            "distance_to_nearest_gene": near[1] if near else None,
        })
        print(f"  {c.rsid:<14} open in {len(open_in)}/{len(meta)} kidney samples"
              f"  |  {', '.join(feats) if feats else 'no annotated element'}")

    d = pd.DataFrame(rows)
    d.to_csv(RESULTS / "candidate_annotation.csv", index=False)

    ind = d[d["independent"]]
    print("\n" + "=" * 78)
    print("THE INDEPENDENT CANDIDATES, ANNOTATED")
    print("=" * 78)
    for r in ind.sort_values("n_kidney_samples_open", ascending=False).itertuples():
        print(f"\n  {r.rsid}   chr22:{r.pos:,}   AF_afr {r.af_afr:.3f}   "
              f"r2 vs gene {r.max_r2_vs_gene:.3f}")
        print(f"    open in            : {r.open_in}")
        print(f"    annotated element  : {r.regulatory_features}")
        print(f"    nearest gene       : {r.nearest_protein_coding_gene} "
              f"({r.distance_to_nearest_gene:,} bp)")

    n_reg = int((ind["regulatory_features"] != "none annotated").sum())
    n_open = int((ind["n_kidney_samples_open"] > 0).sum())
    print("\n" + "-" * 78)
    print(f"  independent candidates in an annotated regulatory element : "
          f"{n_reg} of {len(ind)}")
    print(f"  independent candidates open in >=1 kidney sample          : "
          f"{n_open} of {len(ind)}")
    not_apol1 = ind[ind["nearest_protein_coding_gene"] != "APOL1"]
    if len(not_apol1):
        print(f"\n  NOTE: {len(not_apol1)} of {len(ind)} are nearer a gene other than")
        print("  APOL1. Their assumed target must be stated, not assumed:")
        for r in not_apol1.itertuples():
            print(f"    {r.rsid:<14} nearest {r.nearest_protein_coding_gene} "
                  f"at {r.distance_to_nearest_gene:,} bp")

    out = {
        "candidates": int(len(d)),
        "independent": int(len(ind)),
        "encode_files": {a: {k: v for k, v in m.items() if k != "path"}
                         for a, m in meta.items()},
        "independent_in_annotated_element": n_reg,
        "independent_open_in_any_kidney_sample": n_open,
        "caveat": ("Annotated regulatory elements are predictions from chromatin "
                   "marks across many cell types, not functional evidence that "
                   "the element regulates APOL1 or that the variant alters it."),
    }
    (RESULTS / "candidate_annotation.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/candidate_annotation.csv and .json")


if __name__ == "__main__":
    main()
