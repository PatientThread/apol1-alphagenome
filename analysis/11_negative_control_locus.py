#!/usr/bin/env python3
"""
11_negative_control_locus.py

APOL1 paper, step 11: run the whole pipeline somewhere it should not work.

THE QUESTION. Step 5 found that variants ranked by predicted kidney accessibility
fall inside measured kidney open chromatin four to eight times more often than
the other common variants at the APOL1 locus. That is presented as evidence the
ranking finds real regulatory positions. It could instead be an artefact of the
pipeline: of how peaks are merged, of how the locus-matched background is built,
or of the model simply preferring positions that are open in every tissue and
therefore open in kidney too.

The way to tell is to run the identical pipeline at a locus where the same
enrichment has no business appearing, and see what comes out.

THE CONTROL LOCUS. The beta-globin cluster on chromosome 11, matched to the APOL1
locus in the ways that matter:

    same physical length, 247,688 bp, so peak density and background are
      comparable rather than eyeballed
    a paralogue cluster, HBB, HBD, HBG1, HBG2, HBE1, structurally parallel to
      the APOL paralogues at the test locus
    strong, well characterised regulatory architecture in the locus control
      region, so this is NOT a gene desert; it is a locus with real enhancers
      that belong to a different lineage entirely
    erythroid, not renal, and not expressed in kidney

The last two points together make this a fair test. A gene desert would fail
trivially and prove nothing. A locus with powerful enhancers for the wrong
tissue asks the pipeline the actual question: can you tell the difference?

WHAT EACH OUTCOME MEANS, decided before running.

  little or no enrichment      the APOL1 enrichment is specific and the step 5
                               result survives this challenge
  similar enrichment           the pipeline produces enrichment anywhere, and
                               the APOL1 number says nothing about APOL1; report
                               both and drop any locus-specific language
  stronger enrichment          something in the pipeline manufactures the
                               result; step 5 must be withdrawn

Note that the middle outcome is not fatal to the paper. It would mean the
accessibility channel finds generally-open regulatory positions, which is a
reportable and useful thing for it to do, just a much weaker claim than the one
step 5 currently implies. The paper would have to say so.

Output: data/control_locus_variants.tsv.gz
        results/scored_control_locus/*.tsv.gz
        results/negative_control_locus.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
OUT = RESULTS / "scored_control_locus"
PEAKS = DATA / "control_locus_peaks"
for d in (OUT, PEAKS):
    d.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
from ag_auth import get_key  # noqa: E402

# Matched to the APOL1 locus length exactly: 36,388,018 - 36,140,330 = 247,688.
CHROM = "chr11"
START, END = 5_150_000, 5_397_688
CENTRE = (START + END) // 2
SEQ_1MB = 1_048_576
MIN_AF_AFR = 0.01
SLEEP, MAX_RETRIES, CHUNK = 0.15, 3, 200

KIDNEY = {"glomerular visceral epithelial cell", "kidney"}
PODO_COL = "dnase_glomerular_visceral_epithelial_cell"
TOP_N = (25, 50, 100, 200)
UA = {"User-Agent": "apol1-research/1.0", "Accept": "application/json"}

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
QUERY = """
query RegionVariants($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id rsid pos ref alt
      genome { populations { id ac an } }
    }
  }
}
"""


def gql(retries: int = 4) -> dict:
    body = json.dumps({"query": QUERY,
                       "variables": {"chrom": CHROM.replace("chr", ""),
                                     "start": START, "stop": END}}).encode()
    for a in range(retries):
        try:
            req = urllib.request.Request(
                GNOMAD_API, data=body,
                headers={"Content-Type": "application/json", **UA})
            with urllib.request.urlopen(req, timeout=600) as r:
                return json.loads(r.read().decode())
        except Exception:                                     # noqa: BLE001
            time.sleep(10 * (a + 1))
    raise SystemExit("gnomAD failed repeatedly")


def fetch_variants() -> pd.DataFrame:
    out = DATA / "control_locus_variants.tsv.gz"
    if out.exists():
        return pd.read_csv(out, sep="\t", low_memory=False)
    print(f"querying gnomAD for {CHROM}:{START:,}-{END:,} ...")
    d = gql()
    vs = d["data"]["region"]["variants"]
    rows = []
    for v in vs:
        g = v.get("genome") or {}
        af = None
        for p in (g.get("populations") or []):
            if p["id"] == "afr" and p["an"]:
                af = p["ac"] / p["an"]
        rows.append({"variant_id": v["variant_id"], "rsid": v.get("rsid"),
                     "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
                     "af_afr": af})
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"  {len(df):,} variants; {(df['af_afr'] >= MIN_AF_AFR).sum():,} "
          f"common in AFR")
    return df


def fetch_peaks() -> list[tuple[int, int]]:
    """Measured kidney DNase peaks over the CONTROL locus. Same source as step 5."""
    idx = PEAKS / "index.json"
    if idx.exists():
        graph = json.loads(idx.read_text())
    else:
        url = ("https://www.encodeproject.org/search/?type=File&file_format=bed"
               "&output_type=peaks&assay_title=DNase-seq"
               "&biosample_ontology.term_name=kidney"
               "&assembly=GRCh38&status=released"
               "&limit=8&format=json&field=accession&field=href")
        with urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=180) as fh:
            graph = json.load(fh)["@graph"]
        idx.write_text(json.dumps(graph))
    iv = []
    for f in graph[:8]:
        p = PEAKS / f"{f['accession']}.bed.gz"
        if not p.exists():
            with urllib.request.urlopen(
                    urllib.request.Request(
                        "https://www.encodeproject.org" + f["href"],
                        headers=UA), timeout=300) as fh, open(p, "wb") as o:
                o.write(fh.read())
            time.sleep(0.5)
        with gzip.open(p, "rt") as fh:
            for line in fh:
                c = line.split("\t")
                if c and c[0] == CHROM:
                    s, e = int(c[1]), int(c[2])
                    if e > START and s < END:
                        iv.append((s, e))
    if not iv:
        return []
    iv.sort()
    merged = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


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


def score(v: pd.DataFrame) -> None:
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome
    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    sc = variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"]
    rows, failures, t0, chunk_no = [], 0, time.time(), 0
    for i, r in enumerate(v.itertuples()):
        pos = int(r.pos)
        start = max(0, CENTRE - SEQ_1MB // 2)
        got = None
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                got = client.score_variant(
                    interval=genome.Interval(chromosome=CHROM, start=start,
                                             end=start + SEQ_1MB),
                    variant=genome.Variant(chromosome=CHROM, position=pos,
                                           reference_bases=r.ref,
                                           alternate_bases=r.alt),
                    variant_scorers=[sc])
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
        X = np.asarray(ad.X)
        rec = {"rsid": r.rsid, "pos": pos, "af_afr": r.af_afr}
        for bs in KIDNEY:
            cols = np.where(ad.var["biosample_name"].values == bs)[0]
            if len(cols):
                rec[f"dnase_{bs.replace(' ', '_')}"] = float(np.nanmean(X[:, cols]))
        rows.append(rec)
        time.sleep(SLEEP)
        if len(rows) >= CHUNK:                    # monotonic counter, see step 10
            chunk_no += 1
            p = OUT / f"chunk_{chunk_no:03d}.tsv.gz"
            if p.exists():
                raise SystemExit(f"refusing to overwrite {p}")
            pd.DataFrame(rows).to_csv(p, sep="\t", index=False, compression="gzip")
            rows = []
        if (i + 1) % 100 == 0:
            el = time.time() - t0
            print(f"  {i+1}/{len(v)}  {el:.0f}s  "
                  f"eta {(el/(i+1))*(len(v)-i-1)/60:.0f}m", flush=True)
    if rows:
        chunk_no += 1
        pd.DataFrame(rows).to_csv(OUT / f"chunk_{chunk_no:03d}.tsv.gz",
                                  sep="\t", index=False, compression="gzip")
    written = sum(len(pd.read_csv(f, sep="\t")) for f in OUT.glob("chunk_*.tsv.gz"))
    print(f"\nscored {written}/{len(v)}, {failures} failures")
    if written != len(v) - failures:
        raise SystemExit("row count mismatch; do not use this output")


def analyse(peaks: list[tuple[int, int]]) -> dict:
    d = pd.concat([pd.read_csv(f, sep="\t")
                   for f in sorted(OUT.glob("chunk_*.tsv.gz"))], ignore_index=True)
    d = d.dropna(subset=[PODO_COL]).drop_duplicates("rsid")
    d["k"] = d[PODO_COL].abs()
    d = d.sort_values("k", ascending=False).reset_index(drop=True)
    hit = d["pos"].map(lambda p: in_any(int(p), peaks)).values
    bp = sum(e - s for s, e in peaks)
    print(f"\n  control locus: {len(d):,} scored variants, {len(peaks)} kidney "
          f"peaks covering {bp:,} bp ({100*bp/(END-START):.1f}% of locus)")
    rows = []
    for n in TOP_N:
        top, rest = hit[:n], hit[n:]
        a, b = int(top.sum()), int(n - top.sum())
        c, e2 = int(rest.sum()), int(len(rest) - rest.sum())
        if (c + e2) == 0 or n >= len(hit):
            continue
        orr, p = stats.fisher_exact([[a, b], [c, e2]], alternative="greater")
        rows.append({"top_n": n, "in_peak": a, "pct": round(100 * a / n, 1),
                     "background_pct": round(100 * c / (c + e2), 1),
                     "odds_ratio": round(float(orr), 3), "p": float(p)})
        print(f"    top {n:<4} {a}/{n} in peak ({100*a/n:.1f}%)  "
              f"background {100*c/(c+e2):.1f}%  OR {orr:.2f}  p={p:.2g}")
    return {"locus": f"{CHROM}:{START}-{END}", "bp": END - START,
            "scored": int(len(d)), "peaks": len(peaks), "peak_bp": bp,
            "enrichment": rows}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--analyse-only", action="store_true")
    args = ap.parse_args()

    v = fetch_variants()
    # Same filter as script 02, applied to BOTH loci: common in Africa, single
    # nucleotide substitutions only, de-duplicated. Indels are excluded because
    # script 02 excluded them, not because they do not matter.
    v = v[v["af_afr"].fillna(0) >= MIN_AF_AFR].copy()
    v = v[(v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
    v = v.drop_duplicates(subset=["pos", "ref", "alt"]).sort_values("pos")
    v = v.reset_index(drop=True)

    # MATCH THE TEST LOCUS. The control locus carries far more common variants
    # than APOL1 (3,492 vs 1,930), so "top 25" would mean a different percentile
    # at each locus and the comparison would not be like for like. Subsample to
    # the same count, stratified on African allele frequency so the frequency
    # spectra match too, since rare and common variants differ systematically in
    # how often they fall in open chromatin.
    target = pd.read_csv(DATA / "apol1_locus_variants.tsv.gz", sep="\t",
                         low_memory=False)
    target = target[target["af_afr"].fillna(0) >= MIN_AF_AFR]
    target = target[(target["ref"].str.len() == 1)
                    & (target["alt"].str.len() == 1)]
    target = target.drop_duplicates(subset=["pos", "ref", "alt"])
    n_target = len(target)
    if len(v) > n_target:
        bins = [0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.01]
        tb = pd.cut(target["af_afr"], bins).value_counts()
        rng = np.random.default_rng(20260915)
        keep = []
        vb = pd.cut(v["af_afr"], bins)
        for b, want in tb.items():
            pool = v.index[vb == b].to_numpy()
            take = min(int(want), len(pool))
            if take:
                keep.extend(rng.choice(pool, size=take, replace=False))
        v = v.loc[sorted(keep)].reset_index(drop=True)
        print(f"  frequency-matched subsample to {len(v):,} "
              f"(APOL1 locus has {n_target:,})")
    if args.limit:
        v = v.head(args.limit)
    print(f"variants to score at the control locus: {len(v):,}")

    peaks = fetch_peaks()
    if not args.analyse_only:
        score(v)
    res = analyse(peaks)

    apol1 = json.loads((RESULTS / "encode_validation.json").read_text())
    a_top = apol1["support_1"]["tests"]
    res["apol1_comparison"] = {
        f"top{n}": a_top.get(f"dnase_glomerular_visceral_epithelial_cell_top{n}",
                             {}).get("odds_ratio") for n in TOP_N}
    print("\n" + "=" * 70)
    print("CONTROL LOCUS vs APOL1 LOCUS, same pipeline, same peak source")
    print("=" * 70)
    print(f"  {'':<10}{'control OR':>13}{'APOL1 OR':>11}")
    for r in res["enrichment"]:
        n = r["top_n"]
        a = res["apol1_comparison"].get(f"top{n}")
        print(f"  top {n:<6}{r['odds_ratio']:>13.2f}{(a if a else float('nan')):>11.2f}")
    (RESULTS / "negative_control_locus.json").write_text(json.dumps(res, indent=2))
    print(f"\n  Wrote {RESULTS}/negative_control_locus.json")


if __name__ == "__main__":
    main()
