#!/usr/bin/env python3
"""
09_cell_type_resolved_validation.py

APOL1 paper, step 9: the validation step 5 should have been.

THE GAP. Step 5 ranked variants by PREDICTED PODOCYTE accessibility and then
tested them against MEASURED WHOLE-KIDNEY accessibility, because its ENCODE query
was pinned to biosample_ontology.term_name=kidney. All eight peak files it
retrieved are whole kidney tissue. The enrichment it found is real, but it cannot
say anything about cell type, and cell type is the whole point: APOL1 nephropathy
is a podocytopathy, and a prediction specific to podocytes is a different claim
from a prediction that something is open somewhere in the kidney.

ENCODE does hold cell-type-resolved kidney accessibility. Step 5 simply did not
ask for it.

THE TEST. Rank by predicted podocyte accessibility, exactly as before, then ask
how enriched the top-ranked variants are in measured peaks from each kidney cell
type separately.

    If the ranking carries genuine cell-type information, enrichment should be
    strongest in measured PODOCYTE peaks, because that is what was predicted.

    If enrichment is the same in podocytes, proximal tubule, kidney epithelium
    and whole kidney, then the accessibility channel is picking out positions
    that are open in kidney generally, not podocyte-specific regulatory
    positions, and the paper must say so.

The second outcome would echo the sister paper on this model, which found that
swapping the tissue output barely changed kidney eQTL predictions. Finding it
again here, by a different route and on a different task, would strengthen that
result rather than weaken this one. Both outcomes are reportable and the wording
above is fixed before the numbers are seen.

BACKGROUND. As in step 5, the comparator is the other common variants at the
same locus, not the genome. Everything here shares sequence composition, gene
density and mappability, so the test isolates the ranking.

Output: results/cell_type_validation.json
        results/cell_type_enrichment.csv
        results/candidate_cell_types.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import glob
import urllib.parse
import gzip
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
CT_PEAKS = DATA / "encode_celltype_peaks"
CT_PEAKS.mkdir(parents=True, exist_ok=True)

CHROM = "chr22"
LOCUS_START, LOCUS_END = 36_140_330, 36_388_018
PODO_SCORE = "dnase_glomerular_visceral_epithelial_cell"
UA = {"User-Agent": "apol1-research/1.0", "Accept": "application/json"}

# Kidney cell types ENCODE exposes as DNase-seq peak files, plus whole kidney as
# the step-5 comparator. Podocyte first: it is what the model predicted.
CELL_TYPES = {
    "podocyte": "glomerular visceral epithelial cell",
    "proximal tubule": "epithelial cell of proximal tubule",
    "kidney epithelial": "kidney epithelial cell",
    "renal cortical epithelial": "renal cortical epithelial cell",
    "kidney tubule": "kidney tubule cell",
    "whole kidney": "kidney",
}
FILES_PER_TYPE = 4
TOP_N = (25, 50, 100, 200)


def fetch_files(term: str, tag: str) -> list[Path]:
    idx = CT_PEAKS / f"{tag}.index.json"
    if idx.exists():
        graph = json.loads(idx.read_text())
    else:
        url = ("https://www.encodeproject.org/search/?type=File&file_format=bed"
               "&output_type=peaks&assay_title=DNase-seq"
               f"&biosample_ontology.term_name={urllib.parse.quote(term)}"
               "&assembly=GRCh38&status=released"
               f"&limit={FILES_PER_TYPE}&format=json&field=accession&field=href")
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=180) as fh:
            graph = json.load(fh)["@graph"]
        idx.write_text(json.dumps(graph))
        time.sleep(0.5)

    paths = []
    for f in graph[:FILES_PER_TYPE]:
        p = CT_PEAKS / f"{f['accession']}.bed.gz"
        if not p.exists():
            req = urllib.request.Request(
                "https://www.encodeproject.org" + f["href"], headers=UA)
            with urllib.request.urlopen(req, timeout=300) as fh, open(p, "wb") as o:
                o.write(fh.read())
            time.sleep(0.5)
        paths.append(p)
    return paths


def peak_union(paths: list[Path]) -> list[tuple[int, int]]:
    """Merged peaks on the locus. Union across replicates of one cell type."""
    iv = []
    for p in paths:
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.split("\t")
                if f and f[0] == CHROM:
                    s, e = int(f[1]), int(f[2])
                    if e > LOCUS_START and s < LOCUS_END:
                        iv.append((s, e))
    if not iv:
        return []
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
    scored = pd.concat(
        [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(RESULTS / "scored" / "*.tsv.gz")))],
        ignore_index=True)
    scored = scored.dropna(subset=[PODO_SCORE]).drop_duplicates("rsid")
    scored["rank_key"] = scored[PODO_SCORE].abs()
    scored = scored.sort_values("rank_key", ascending=False).reset_index(drop=True)
    print(f"scored variants at the locus: {len(scored):,}\n")

    peaks, rows = {}, []
    for label, term in CELL_TYPES.items():
        paths = fetch_files(term, label.replace(" ", "_"))
        iv = peak_union(paths)
        peaks[label] = iv
        bp = sum(e - s for s, e in iv)
        print(f"  {label:<26} {len(paths)} files  {len(iv):>4} peaks  "
              f"{bp:>7,} bp  {100*bp/(LOCUS_END-LOCUS_START):>5.1f}% of locus")

    print()
    for label, iv in peaks.items():
        if not iv:
            print(f"  {label}: no peaks on this locus, skipped")
            continue
        hit = scored["pos"].map(lambda p: in_any(int(p), iv)).values
        for n in TOP_N:
            top = hit[:n]
            rest = hit[n:]
            a, b = int(top.sum()), int(n - top.sum())
            c, d = int(rest.sum()), int(len(rest) - rest.sum())
            orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({"cell_type": label, "top_n": n,
                         "in_peak": a, "pct": round(100 * a / n, 1),
                         "background_pct": round(100 * c / (c + d), 1),
                         "odds_ratio": round(float(orr), 3), "p": float(p)})

    e = pd.DataFrame(rows)
    e.to_csv(RESULTS / "cell_type_enrichment.csv", index=False)

    print("=" * 78)
    print("ENRICHMENT OF PODOCYTE-PREDICTED TOP VARIANTS, BY MEASURED CELL TYPE")
    print("=" * 78)
    print(f"  {'cell type':<26}" + "".join(f"{'top'+str(n):>11}" for n in TOP_N))
    for label in peaks:
        sub = e[e["cell_type"] == label]
        if sub.empty:
            continue
        cells = "".join(f"{sub[sub.top_n==n]['odds_ratio'].iloc[0]:>11.2f}"
                        for n in TOP_N)
        print(f"  {label:<26}{cells}")
    print("\n  (odds ratio against the other common variants at the same locus)")

    # ---- the comparison that decides the claim -----------------------------
    podo = e[(e.cell_type == "podocyte")].set_index("top_n")["odds_ratio"]
    others = (e[(e.cell_type != "podocyte") & (e.cell_type != "whole kidney")]
              .groupby("top_n")["odds_ratio"].median())
    print("\n  podocyte vs the median of other kidney cell types:")
    for n in TOP_N:
        if n in podo.index and n in others.index:
            print(f"    top {n:<4} podocyte {podo[n]:>6.2f}   others {others[n]:>6.2f}"
                  f"   ratio {podo[n]/others[n]:>5.2f}" if others[n] else "")
    podo_better = sum(1 for n in TOP_N
                      if n in podo.index and n in others.index and podo[n] > others[n])

    print("\n  READING:")
    if podo_better == len(TOP_N):
        print("  Enrichment is strongest in measured podocyte chromatin at every")
        print("  threshold. The ranking carries podocyte-specific information.")
    elif podo_better == 0:
        print("  Podocyte enrichment is NOT stronger than other kidney cell types.")
        print("  The channel selects positions open in kidney generally, not")
        print("  podocyte-specific ones. Say so, and drop any cell-type claim.")
    else:
        print(f"  Mixed: podocyte leads at {podo_better} of {len(TOP_N)} thresholds.")
        print("  Too weak to claim cell-type specificity. Report the table.")

    # ---- which cell types are the independent candidates open in -----------
    ann = pd.read_csv(RESULTS / "ld_gene_haplotype_per_variant.csv")
    ind = ann[ann["max_r2_vs_gene"] < 0.2]
    crows = []
    for c in ind.itertuples():
        openin = [lab for lab, iv in peaks.items() if iv and in_any(int(c.pos), iv)]
        crows.append({"rsid": c.rsid, "pos": c.pos, "af_afr": c.af_afr,
                      "max_r2_vs_gene": c.max_r2_vs_gene,
                      "n_cell_types_open": len(openin),
                      "open_in_cell_types": "; ".join(openin) or "none"})
    cdf = pd.DataFrame(crows).sort_values("n_cell_types_open", ascending=False)
    cdf.to_csv(RESULTS / "candidate_cell_types.csv", index=False)

    print("\n" + "=" * 78)
    print("THE INDEPENDENT CANDIDATES, BY MEASURED CELL TYPE")
    print("=" * 78)
    for r in cdf.itertuples():
        print(f"  {r.rsid:<14} AF {r.af_afr:.3f}  open in: {r.open_in_cell_types}")

    out = {"scored_variants": int(len(scored)),
           "cell_types": {k: {"files": FILES_PER_TYPE, "peaks": len(v),
                              "bp": sum(e2 - s for s, e2 in v)}
                          for k, v in peaks.items()},
           "ranking_channel": PODO_SCORE,
           "background": "other common variants at the same locus",
           "podocyte_leads_at_n_thresholds": podo_better,
           "n_thresholds": len(TOP_N)}
    (RESULTS / "cell_type_validation.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/cell_type_validation.json")


if __name__ == "__main__":
    import urllib.parse  # noqa: E402  (used in fetch_files)
    main()
