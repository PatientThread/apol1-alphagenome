#!/usr/bin/env python3
"""
21_track_manifest.py

Records which model output tracks each named tissue ranking actually used, and
how several tracks were combined into one number. Asked for at review, because
CenterMaskScorer returns one value per output track and the manuscript named
seven "tissue outputs" without saying what each one was.

THE RULE, as implemented in steps 02 and 10: tracks are selected by EXACT match
on `biosample_name`, and where more than one track matches, the score is their
`numpy.nanmean`. This is not one track per ranking.

One variant is scored to enumerate the tracks. The scores are discarded; only
the track metadata is kept.

Output: results/track_manifest.csv, results/track_manifest.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

CHROM, SEQ_1MB = "chr22", 1_048_576
PROBE_POS, PROBE_REF, PROBE_ALT = 36_265_860, "A", "G"   # within the locus
BIOSAMPLES = ["glomerular visceral epithelial cell", "kidney", "hepatocyte",
              "liver", "lung", "brain", "stomach"]
RANK_NAME = {"glomerular visceral epithelial cell": "kidney (podocyte)",
             "kidney": "kidney (whole)", "hepatocyte": "hepatocyte",
             "liver": "liver", "lung": "lung", "brain": "brain",
             "stomach": "stomach"}


def main() -> None:
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    client = dna_client.create(get_key())
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]
    start = max(0, PROBE_POS - SEQ_1MB // 2)
    out = client.score_variant(
        interval=genome.Interval(chromosome=CHROM, start=start,
                                 end=start + SEQ_1MB),
        variant=genome.Variant(chromosome=CHROM, position=PROBE_POS,
                               reference_bases=PROBE_REF,
                               alternate_bases=PROBE_ALT),
        variant_scorers=[scorer])
    ad = out[0]
    var = ad.var
    print(f"DNase output tracks returned: {var.shape[0]}")
    print(f"metadata columns: {list(var.columns)}")

    rows = []
    for bs in BIOSAMPLES:
        cols = np.where(var["biosample_name"].values == bs)[0]
        for c in cols:
            r = {"ranking": RANK_NAME[bs], "biosample_name": bs,
                 "track_index": int(c)}
            for f in ("name", "track_name", "ontology_curie", "biosample_type",
                      "strand", "assay", "data_source"):
                if f in var.columns:
                    r[f] = var.iloc[c][f]
            rows.append(r)
    tm = pd.DataFrame(rows)
    tm.to_csv(RESULTS / "track_manifest.csv", index=False)

    summary = (tm.groupby("ranking").size().rename("n_tracks")
               .reset_index().sort_values("n_tracks", ascending=False))
    print("\ntracks per ranking (combined by numpy.nanmean where n > 1)")
    print(summary.to_string(index=False))

    (RESULTS / "track_manifest.json").write_text(json.dumps({
        "selection_rule": "exact match on biosample_name in the DNase output "
                          "metadata",
        "combination_rule": "numpy.nanmean across all matching tracks; a "
                            "ranking with n>1 is a mean, not a single track",
        "total_dnase_tracks_returned": int(var.shape[0]),
        "tracks_per_ranking": {r["ranking"]: int(r["n_tracks"])
                               for _, r in summary.iterrows()},
        "probe_variant": f"{CHROM}:{PROBE_POS}{PROBE_REF}>{PROBE_ALT}",
        "note": "scores from the probe variant are discarded; only track "
                "metadata is retained",
    }, indent=2))
    print("\nwrote track_manifest.csv and track_manifest.json")


if __name__ == "__main__":
    main()
