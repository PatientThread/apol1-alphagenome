#!/usr/bin/env python3
"""
01_fetch_locus_variants.py

APOL1 paper, step 1: assemble the variant set for the APOL1 locus.

THE QUESTION. Only about 15 to 20% of people carrying two APOL1 high-risk
alleles develop kidney disease. The field agrees the mechanism is expression
DOSE: risk-variant APOL1 is cytotoxic in a dose-dependent way in cells and in
podocyte-specific mice. The field has characterised the STIMULI that raise
expression (interferon above all) in detail. It has never asked whether the
RESPONSE ELEMENTS THEMSELVES differ between carriers.

That is the gap. The KDIGO consensus report on APOL1 mentions expression 22
times and does not once mention an eQTL, a regulatory variant, an enhancer, or a
genome-wide association study. It then explicitly invites the work: "Genetic or
environmental modifiers of APOL1 expression or activity could be used as
components of biomarker panels."

WHY A SEQUENCE MODEL IS THE RIGHT TOOL HERE, and this is not special pleading.
GTEx holds 133 significant APOL1 expression QTLs and NOT ONE is in kidney. The
tissue-QTL route is underpowered at this locus precisely because kidney sample
sizes are the smallest in the resource, which is the finding of the companion
paper. A model that predicts from sequence does not need eQTL sample size.

LOCUS. APOL1 sits at chr22:36,253,010-36,267,530 (GRCh38). The functional unit is
larger: APOL3 through MYH9 spans roughly chr22:36,140,330-36,388,018, about
248 kb, and MYH9 matters because the original kidney-disease association was
mapped there before being resolved to APOL1. The whole cluster plus generous
flanking regulatory space fits inside ONE 1,048,576 bp model window centred near
chr22:36,260,000, so the entire analysis is one context.

WHAT THIS SCRIPT DOES. Pulls every variant gnomAD records in the cluster, with
African/African-American allele frequencies, and tags the three variants that
define the risk haplotypes:
    G1  rs73885319 (p.S342G) and rs60910145 (p.I384M), on one haplotype
    G2  rs71785313 (6 bp deletion removing N388 and Y389)
plus rs2239785 (p.N264K), the protective modifier that abolishes G2 risk.

Output: data/apol1_locus_variants.tsv.gz, results/locus_summary.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DATA.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

# GRCh38
CHROM = "chr22"
APOL1_START, APOL1_END = 36_253_010, 36_267_530
CLUSTER_START, CLUSTER_END = 36_140_330, 36_388_018
WINDOW_CENTRE = 36_260_000
WINDOW = 1_048_576

RISK_VARIANTS = {
    "rs73885319": "G1 (p.S342G)",
    "rs60910145": "G1 (p.I384M)",
    "rs71785313": "G2 (6bp del, p.N388_Y389del)",
    "rs2239785": "N264K (protective modifier)",
}

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

QUERY = """
query RegionVariants($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      rsids
      pos
      ref
      alt
      consequence
      transcript_consequence { gene_symbol major_consequence }
      genome { ac an populations { id ac an } }
      exome   { ac an populations { id ac an } }
    }
  }
}
"""


def gql(query: str, variables: dict, retries: int = 4) -> dict:
    body = json.dumps({"query": query, "variables": variables}).encode()
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                GNOMAD_API, data=body,
                headers={"Content-Type": "application/json",
                         "User-Agent": "apol1-research/1.0"})
            with urllib.request.urlopen(req, timeout=300) as fh:
                out = json.load(fh)
            if "errors" in out and out["errors"]:
                raise RuntimeError(str(out["errors"])[:300])
            return out
        except Exception as exc:                                # noqa: BLE001
            if attempt == retries:
                raise
            wait = 4 * attempt
            print(f"  {type(exc).__name__}, retry {attempt}/{retries} in {wait}s")
            time.sleep(wait)
    raise RuntimeError("unreachable")


def af(block: dict | None, pop_prefix: str) -> float | None:
    """Allele frequency for a population prefix, from a genome/exome block."""
    if not block:
        return None
    pops = block.get("populations") or []
    ac = an = 0
    for p in pops:
        if p["id"].lower() == pop_prefix:
            ac += p["ac"] or 0
            an += p["an"] or 0
    return (ac / an) if an else None


def main() -> None:
    cache = DATA / "gnomad_region_raw.json"
    if cache.exists():
        payload = json.loads(cache.read_text())
        print("using cached gnomAD response")
    else:
        print(f"querying gnomAD for {CHROM}:{CLUSTER_START}-{CLUSTER_END} ...")
        payload = gql(QUERY, {"chrom": CHROM.replace("chr", ""),
                              "start": CLUSTER_START, "stop": CLUSTER_END})
        cache.write_text(json.dumps(payload))
    variants = payload["data"]["region"]["variants"]
    print(f"gnomAD returned {len(variants):,} variants")

    rows = []
    for v in variants:
        tc = v.get("transcript_consequence") or {}
        g, e = v.get("genome"), v.get("exome")
        afr = af(g, "afr")
        if afr is None:
            afr = af(e, "afr")
        tot_ac = (g or {}).get("ac", 0) or 0
        tot_an = (g or {}).get("an", 0) or 0
        if not tot_an and e:
            tot_ac, tot_an = e.get("ac", 0) or 0, e.get("an", 0) or 0
        rsids = v.get("rsids") or []
        rows.append({
            "variant_id": v["variant_id"],
            "rsid": rsids[0] if rsids else None,
            "pos": v["pos"],
            "ref": v["ref"],
            "alt": v["alt"],
            "consequence": v.get("consequence"),
            "gene": tc.get("gene_symbol"),
            "af_global": (tot_ac / tot_an) if tot_an else None,
            "af_afr": afr,
            "in_apol1_gene": APOL1_START <= v["pos"] <= APOL1_END,
        })
    df = pd.DataFrame(rows)
    df["is_risk_variant"] = df["rsid"].map(RISK_VARIANTS)

    df.to_csv(DATA / "apol1_locus_variants.tsv.gz", sep="\t", index=False,
              compression="gzip")

    coding_terms = ("missense", "synonymous", "stop", "frameshift",
                    "inframe", "start_lost", "splice")
    is_coding = df["consequence"].fillna("").str.contains(
        "|".join(coding_terms), case=False)

    found = df[df["is_risk_variant"].notna()]
    summary = {
        "region": f"{CHROM}:{CLUSTER_START}-{CLUSTER_END}",
        "region_bp": CLUSTER_END - CLUSTER_START,
        "model_window": {"centre": WINDOW_CENTRE, "width": WINDOW,
                         "covers_cluster": (WINDOW // 2) >= max(
                             WINDOW_CENTRE - CLUSTER_START,
                             CLUSTER_END - WINDOW_CENTRE)},
        "variants_total": int(len(df)),
        "variants_in_apol1_gene": int(df["in_apol1_gene"].sum()),
        "coding": int(is_coding.sum()),
        "non_coding": int((~is_coding).sum()),
        "non_coding_pct": round(100 * (~is_coding).sum() / len(df), 1),
        "common_in_afr_1pct": int((df["af_afr"].fillna(0) >= 0.01).sum()),
        "common_in_afr_5pct": int((df["af_afr"].fillna(0) >= 0.05).sum()),
        "risk_variants_found": found[["rsid", "is_risk_variant", "pos",
                                      "af_afr"]].to_dict(orient="records"),
    }
    (RESULTS / "locus_summary.json").write_text(json.dumps(summary, indent=2))

    print()
    print("=" * 72)
    print("APOL1 LOCUS VARIANT SET")
    print("=" * 72)
    print(f"  region                     : {summary['region']} "
          f"({summary['region_bp']:,} bp)")
    print(f"  fits one model window      : {summary['model_window']['covers_cluster']}")
    print(f"  variants in gnomAD         : {summary['variants_total']:,}")
    print(f"    within APOL1 gene body   : {summary['variants_in_apol1_gene']:,}")
    print(f"    coding                   : {summary['coding']:,}")
    print(f"    NON-CODING               : {summary['non_coding']:,} "
          f"({summary['non_coding_pct']}%)")
    print(f"  common in African ancestry : "
          f"{summary['common_in_afr_1pct']:,} at >=1%, "
          f"{summary['common_in_afr_5pct']:,} at >=5%")
    print()
    print("  RISK-DEFINING VARIANTS:")
    for r in summary["risk_variants_found"]:
        a = f"{r['af_afr']:.3f}" if r["af_afr"] else "n/a"
        print(f"    {r['rsid']:<14} {r['is_risk_variant']:<32} "
              f"pos {r['pos']:,}  AFR AF {a}")
    missing = set(RISK_VARIANTS) - {r["rsid"] for r in summary["risk_variants_found"]}
    if missing:
        print(f"    NOT FOUND (check manually): {sorted(missing)}")
    print()
    print(f"  Wrote {DATA}/apol1_locus_variants.tsv.gz")


if __name__ == "__main__":
    main()
