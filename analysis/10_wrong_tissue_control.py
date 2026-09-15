#!/usr/bin/env python3
"""
10_wrong_tissue_control.py

APOL1 paper, step 10: the control that decides whether the kidney label matters.

THE QUESTION. Step 9 showed the ranking is enriched in every kidney cell type at
roughly equal strength, and that podocyte, which is what was actually predicted,
is not the strongest. That already suggests the accessibility channel finds
positions open in kidney generally rather than in one cell type. It does not test
the harder version of the question.

    If we rank the same variants using a DELIBERATELY WRONG tissue, liver or
    brain or lung, are the top-ranked variants still enriched in measured
    KIDNEY chromatin?

If they are, then the enrichment reported in step 5 is not evidence that the
model knows anything about kidney. It is evidence that the model can find
regulatory positions, and the tissue label on the request is close to decorative.
That is the direct analogue of the tissue-swap control in the sister paper, which
found that scoring kidney eQTLs with other tissues' outputs barely changed
performance.

    H0  wrong-tissue rankings are no less enriched in kidney chromatin than the
        kidney ranking, i.e. the label carries no information
    H1  the kidney ranking is more enriched, i.e. the label carries information

Unlike step 7, here the NULL is the uncomfortable outcome, and it is the one the
sister paper's result predicts. Writing that down before running.

WHAT EACH OUTCOME MEANS, decided in advance.

  kidney clearly strongest     the accessibility channel carries genuine
                               tissue information and the step 5 result stands
                               as a kidney finding
  all tissues similar          the channel finds open chromatin, not kidney
                               chromatin; step 5 must be reported as "selects
                               regulatory positions", never as "kidney-specific"
  wrong tissue stronger        the ranking is not measuring what it claims at
                               all and the accessibility result must be withdrawn

DESIGN NOTE. The same variants, the same 1 Mb windows, the same API call, the
same measured kidney peaks and the same locus-matched background. Only the
column selected from the model's output changes. That isolates the tissue label
and nothing else.

Output: results/scored_wrong_tissue/*.tsv.gz
        results/wrong_tissue_control.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
OUT = RESULTS / "scored_wrong_tissue"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

CHROM = "chr22"
SEQ_1MB = 1_048_576
SLEEP = 0.15
MAX_RETRIES = 3
CHUNK = 200

# The kidney reference, carried through so both come from one run and cannot
# differ through some incidental change between runs.
KIDNEY = {"glomerular visceral epithelial cell", "kidney"}

# Deliberately wrong tissues. Chosen to be unambiguously non-renal, to span
# different germ layers, and to be well represented in the model's DNase output.
WRONG = {
    "liver", "hepatocyte",
    "brain", "frontal cortex",
    "lung", "left lung",
    "heart left ventricle",
    "stomach",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    # MUST reproduce script 02's variant set EXACTLY. This control compares one
    # column of the model's output against another on the same variants; if the
    # variant sets differ the comparison is meaningless. Script 02 applies an
    # SNV-only filter and de-duplicates, which takes 2,429 common variants down
    # to 1,930. Omitting that here silently scored a different 2,429.
    v = pd.read_csv(DATA / "apol1_locus_variants.tsv.gz", sep="\t",
                    low_memory=False)
    v = v[v["af_afr"].fillna(0) >= 0.01].copy()
    v = v[(v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
    v = v.drop_duplicates(subset=["pos", "ref", "alt"]).sort_values("pos")
    v = v.reset_index(drop=True)
    if args.limit:
        v = v.head(args.limit)
    print(f"variants to score: {len(v):,}")
    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    dnase_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]

    want = KIDNEY | WRONG
    rows, failures, t0, chunk_no = [], 0, time.time(), 0

    for i, r in enumerate(v.itertuples()):
        pos = int(r.pos)
        start = max(0, pos - SEQ_1MB // 2)
        got = None
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                got = client.score_variant(
                    interval=genome.Interval(chromosome=CHROM, start=start,
                                             end=start + SEQ_1MB),
                    variant=genome.Variant(chromosome=CHROM, position=pos,
                                           reference_bases=r.ref,
                                           alternate_bases=r.alt),
                    variant_scorers=[dnase_scorer])
                break
            except Exception as exc:                          # noqa: BLE001
                if attempt == MAX_RETRIES:
                    failures += 1
                    print(f"  FAIL {r.variant_id}: {type(exc).__name__}")
                else:
                    time.sleep(SLEEP * 4 * attempt)
        if got is None:
            continue

        ad = got[0]
        if ad is None or "biosample_name" not in ad.var.columns:
            raise RuntimeError("no biosample_name in DNase result")
        X = np.asarray(ad.X)
        rec = {"rsid": r.rsid, "pos": pos, "af_afr": r.af_afr}
        for bs in want:
            cols = np.where(ad.var["biosample_name"].values == bs)[0]
            if len(cols):
                rec[f"dnase_{bs.replace(' ', '_')}"] = float(np.nanmean(X[:, cols]))
        rows.append(rec)
        time.sleep(SLEEP)

        # Monotonic counter, never derived from len(rows): rows is emptied on
        # each write, so a derived name collapses every chunk onto one file.
        # That bug destroyed 1,600 of 1,930 rows on the first run of script 02.
        if len(rows) >= CHUNK:
            chunk_no += 1
            p = OUT / f"chunk_{chunk_no:03d}.tsv.gz"
            if p.exists():
                raise SystemExit(f"refusing to overwrite {p}")
            pd.DataFrame(rows).to_csv(p, sep="\t", index=False,
                                      compression="gzip")
            rows = []
        if (i + 1) % 100 == 0:
            el = time.time() - t0
            print(f"  {i+1}/{len(v)}  {el:.0f}s  "
                  f"eta {(el/(i+1))*(len(v)-i-1)/60:.0f}m", flush=True)

    if rows:
        chunk_no += 1
        p = OUT / f"chunk_{chunk_no:03d}.tsv.gz"
        if p.exists():
            raise SystemExit(f"refusing to overwrite {p}")
        pd.DataFrame(rows).to_csv(p, sep="\t", index=False, compression="gzip")

    # Verify by counting what landed on disk, not by reaching the end of a loop.
    written = sum(len(pd.read_csv(f, sep="\t")) for f in OUT.glob("chunk_*.tsv.gz"))
    expected = len(v) - failures
    print(f"\nscored {written}/{len(v)}, {failures} failures")
    if written != expected:
        raise SystemExit(f"row count mismatch: {written} on disk, "
                         f"{expected} expected. Do not use this output.")

    (RESULTS / "wrong_tissue_scoring_log.json").write_text(json.dumps(
        {"variants": int(len(v)), "written": int(written),
         "failures": int(failures), "seconds": round(time.time() - t0, 1),
         "kidney_biosamples": sorted(KIDNEY),
         "wrong_biosamples": sorted(WRONG)}, indent=2))
    print(f"  Wrote {OUT} and wrong_tissue_scoring_log.json")


if __name__ == "__main__":
    main()
