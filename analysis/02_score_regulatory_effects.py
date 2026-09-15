#!/usr/bin/env python3
"""
02_score_regulatory_effects.py

APOL1 paper, step 2: predict the regulatory effect of every common
African-ancestry variant in the APOL1 cluster, in kidney.

WHAT IS SCORED AND WHY THOSE TRACKS. APOL1 nephropathy is a podocytopathy, and
the toxicity is dose-dependent. So the quantities that matter are APOL1
expression in kidney and chromatin accessibility in podocyte. Both exist in the
model, but only just:

    RNA-seq      Kidney_Cortex (adult), Kidney_Medulla (adult), kidney (embryonic)
    DNase        glomerular visceral epithelial cell = PODOCYTE, one track,
                 from a paediatric donor
                 plus proximal tubule, kidney epithelial, tubule cell, whole kidney
    ATAC         whole kidney (adult), one track

**That single podocyte track is the entire podocyte representation in the model.**
It is the sharpest illustration of the companion paper's finding, and it must be
stated as a limitation here rather than buried: the cell type in which APOL1
kills is represented by one accessibility experiment from one child.

TWO SCORERS, because the question has two halves:
    GeneMaskLFCScorer on RNA_SEQ   does the variant change APOL1 expression?
    CenterMaskScorer on DNASE/ATAC does it change accessibility at that point,
                                   which is the mechanism a regulatory variant
                                   would act through?

VARIANT SET. Variants with African-ancestry allele frequency at or above 1% in
the cluster. Rare variants are excluded deliberately: a modifier that explains
why 80% of two-risk-allele carriers stay well has to be common enough to be
present in most of them.

Output: results/scored/<chunk>.tsv.gz, results/apol1_scoring_log.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
OUT = RESULTS / "scored"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

CHROM = "chr22"
SEQ_1MB = 1_048_576          # exact; 1_000_000 is rejected by the model
APOL1_ENSG = "ENSG00000100342"

# Curated kidney tracks. Named explicitly, NOT substring-matched: matching on
# "renal" would pull in adrenal gland, which is a different organ.
KIDNEY_BIOSAMPLES = {
    "cortex of kidney", "outer medulla of kidney", "kidney", "left kidney",
    "right kidney", "kidney epithelial cell", "kidney tubule cell",
    "epithelial cell of proximal tubule", "renal cortical epithelial cell",
    "glomerular visceral epithelial cell",          # podocyte
    "kidney capillary endothelial cell", "nephron progenitor cell",
    "renal cortex interstitium", "left renal cortex interstitium",
    "right renal cortex interstitium",
}
PODOCYTE = "glomerular visceral epithelial cell"

SLEEP = 0.15
MAX_RETRIES = 3
CHUNK = 200


def build_targets(om) -> dict[str, pd.DataFrame]:
    out = {}
    for mod in ("rna_seq", "dnase", "atac"):
        df = getattr(om, mod).copy()
        df["idx"] = np.arange(len(df))
        out[mod] = df[df["biosample_name"].isin(KIDNEY_BIOSAMPLES)]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-af", type=float, default=0.01)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    v = pd.read_csv(DATA / "apol1_locus_variants.tsv.gz", sep="\t")
    v = v[v["af_afr"].fillna(0) >= args.min_af].copy()
    # single-nucleotide substitutions only; indels need separate handling
    v = v[(v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
    v = v.drop_duplicates(subset=["pos", "ref", "alt"]).sort_values("pos")
    if args.limit:
        v = v.head(args.limit)
    print(f"variants to score: {len(v):,} (AFR AF >= {args.min_af})")

    # Never log the key, not even masked. masked() reveals 10 of 39
    # characters and that fragment reached a committed artefact once on
    # the sister project. Confirm authentication, print nothing about it.
    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    om = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    targets = build_targets(om)
    for mod, df in targets.items():
        print(f"  {mod:<9} {len(df)} kidney tracks"
              + ("  (incl. PODOCYTE)" if PODOCYTE in set(df["biosample_name"]) else ""))

    rna_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
    dnase_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]
    atac_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"]

    # chunk_no MUST be an independent monotonic counter. Deriving the filename
    # from len(rows) is wrong because rows is emptied after each write, so every
    # chunk resolves to the same name and silently overwrites its predecessor.
    # That bug destroyed 1,600 of 1,930 rows on the first full run.
    rows, failures, t0, chunk_no = [], 0, time.time(), 0
    for i, r in enumerate(v.itertuples()):
        pos = int(r.pos)
        start = max(0, pos - SEQ_1MB // 2)
        interval = genome.Interval(chromosome=CHROM, start=start,
                                   end=start + SEQ_1MB)
        variant = genome.Variant(chromosome=CHROM, position=pos,
                                 reference_bases=r.ref, alternate_bases=r.alt)
        got = None
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                got = client.score_variant(
                    interval=interval, variant=variant,
                    variant_scorers=[rna_scorer, dnase_scorer, atac_scorer])
                break
            except Exception as exc:                          # noqa: BLE001
                if attempt == MAX_RETRIES:
                    failures += 1
                    print(f"  FAIL {r.variant_id}: {type(exc).__name__}: "
                          f"{str(exc)[:100]}")
                else:
                    time.sleep(SLEEP * 4 * attempt)
        if got is None:
            continue

        rec = {"variant_id": r.variant_id, "rsid": r.rsid, "pos": pos,
               "ref": r.ref, "alt": r.alt, "af_afr": r.af_afr,
               "consequence": r.consequence, "gene": r.gene,
               "is_risk_variant": r.is_risk_variant}

        # The returned AnnData indexes columns POSITIONALLY ('0','1','2'...),
        # not by track name, and returns a SUBSET of tracks. Select columns from
        # the result's own .var by biosample_name. Looking up an index taken
        # from a different DataFrame silently matches nothing, which is how the
        # first version produced empty score columns.
        def collect(ad, prefix: str, want: set[str]) -> None:
            if ad is None or "biosample_name" not in ad.var.columns:
                raise RuntimeError(f"{prefix}: no biosample_name in result var")
            X = np.asarray(ad.X)
            for bs in want:
                cols = np.where(ad.var["biosample_name"].values == bs)[0]
                if len(cols) == 0:
                    continue
                key = f"{prefix}_{bs.replace(' ', '_')}"
                rec[key] = float(np.nanmean(X[:, cols]))

        # RNA: restrict rows to APOL1 before averaging over kidney tracks
        ad = got[0]
        gcol = next((c for c in ("gene_id", "gene_name", "gene")
                     if c in ad.obs.columns), None)
        if gcol is None:
            raise RuntimeError("RNA result has no gene column")
        gm = (ad.obs[gcol].astype(str).str.split(".").str[0] == APOL1_ENSG).values
        if not gm.any():
            failures += 1
            if failures <= 3:
                print(f"  APOL1 not in returned genes for {r.variant_id}; "
                      f"genes present: {list(ad.obs[gcol][:3])}")
            continue
        Xr = np.asarray(ad.X)
        for bs in KIDNEY_BIOSAMPLES:
            cols = np.where(ad.var["biosample_name"].values == bs)[0]
            if len(cols) == 0:
                continue
            rec[f"rna_{bs.replace(' ', '_')}"] = float(
                np.nanmean(Xr[np.where(gm)[0][:, None], cols]))

        collect(got[1], "dnase", KIDNEY_BIOSAMPLES)
        collect(got[2], "atac", KIDNEY_BIOSAMPLES)

        rows.append(rec)
        time.sleep(SLEEP)

        if len(rows) >= CHUNK:
            chunk_no += 1
            dest = OUT / f"chunk_{chunk_no:03d}.tsv.gz"
            if dest.exists():
                raise RuntimeError(f"refusing to overwrite {dest}")
            pd.DataFrame(rows).to_csv(dest, sep="\t", index=False,
                                      compression="gzip")
            rows = []
        if (i + 1) % 100 == 0:
            el = time.time() - t0
            print(f"  {i+1}/{len(v)}  {el:.0f}s  "
                  f"eta {(len(v)-i-1)*el/(i+1)/60:.0f}m")

    if rows:
        chunk_no += 1
        pd.DataFrame(rows).to_csv(OUT / f"chunk_{chunk_no:03d}.tsv.gz",
                                  sep="\t", index=False, compression="gzip")

    # Verify nothing was lost before declaring success.
    written = sum(len(pd.read_csv(f, sep="\t")) for f in OUT.glob("*.tsv.gz"))
    expected = len(v) - failures
    print(f"rows written {written:,} / expected {expected:,}")
    if written != expected:
        raise SystemExit(f"ROW LOSS: {expected - written} rows missing")

    dt = time.time() - t0
    pd.DataFrame([{"variants": len(v), "failures": failures,
                   "seconds": round(dt, 1)}]).to_csv(
        RESULTS / "apol1_scoring_log.csv", index=False)
    print(f"\ndone: {len(v) - failures}/{len(v)} scored, {failures} failures, "
          f"{dt/60:.0f} min")


if __name__ == "__main__":
    main()
