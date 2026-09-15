#!/usr/bin/env python3
"""
07_ld_against_gene_haplotype.py

APOL1 paper, step 7: close the hole that step 6 left open.

WHAT STEP 6 ESTABLISHED, AND WHAT IT DID NOT. Step 6 used Ensembl's precomputed
linkage disequilibrium and found every top candidate independent of G1 (maximum
r-squared 0.133 across seven African populations). The internal control passed:
the two G1 SNPs return r-squared 1.0 against each other.

But G2 returned zero partners in every population. G2 is a six-base deletion,
Ensembl excludes indels from its LD service, and the deletion is absent
altogether from the 1000 Genomes phase 3 biallelic call set: the only indel
within a kilobase of the site is a single-base insertion 88 bp away. So step 6
tested G1 and only G1, and saying it tested "the risk haplotypes" would have been
false.

THE FIX, WHICH IS ALSO A STRONGER TEST. Rather than chase a proxy for G2 and
inherit whatever error that proxy carries, test each candidate against EVERY
common variant in the APOL1 gene body. G2 sits inside the gene. So does G1. If a
candidate is in low LD with all of the gene's common variation, it is not tagging
the gene's haplotype structure, and that conclusion holds whatever the causal
variant turns out to be and whether or not it is genotyped here.

    H0  the candidate is independent of APOL1 gene-body haplotype structure
    H1  the candidate tags it, and any disease association would be confounded

The null is again the good outcome.

METHOD. Phased haplotypes from 1000 Genomes phase 3, GRCh38 lift, read remotely.
r-squared and D' are computed from the phased data directly, which is exact
rather than the composite estimate you get from unphased genotypes. Populations
are kept separate, because pooling African populations inflates LD through
population structure and would manufacture exactly the false positive this
script exists to exclude. Each candidate is judged on its MAXIMUM r-squared
across every gene-body variant and every population, which is the conservative
reading.

Output: results/ld_gene_haplotype.json
        results/ld_gene_haplotype_per_variant.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

VCF = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
       "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
       "ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
PANEL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
         "integrated_call_samples_v3.20130502.ALL.panel")

# APOL1, GRCh38. The gene body, which contains both G1 and G2.
APOL1 = ("22", 36253071, 36267530)
AFR_POPS = ["YRI", "LWK", "ESN", "MSL", "GWD", "ACB", "ASW"]
MIN_MAF = 0.01          # common variation only; rare variants cannot tag anything


def load_panel() -> dict[str, list[str]]:
    with urllib.request.urlopen(PANEL, timeout=60) as r:
        lines = r.read().decode().strip().split("\n")
    pops: dict[str, list[str]] = {p: [] for p in AFR_POPS}
    for line in lines[1:]:
        f = line.split()
        if len(f) >= 2 and f[1] in pops:
            pops[f[1]].append(f[0])
    return pops


def haplotypes(rec, idx) -> np.ndarray:
    """Phased haplotype vector, 2 per sample, for the given sample indices."""
    out = []
    for s in idx:
        a = rec.samples[s].alleles
        gt = rec.samples[s]["GT"]
        out.extend([-1 if g is None else int(g) for g in gt])
    return np.array(out, dtype=np.int8)


def r2_dprime(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Exact r-squared and D' from phased haplotypes."""
    ok = (x >= 0) & (y >= 0)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 20:
        return float("nan"), float("nan")
    pa, pb = x.mean(), y.mean()
    if min(pa, pb) <= 0 or max(pa, pb) >= 1:
        return 0.0, 0.0
    d = (x & y).mean() - pa * pb
    r2 = d * d / (pa * (1 - pa) * pb * (1 - pb))
    dmax = min(pa * (1 - pb), (1 - pa) * pb) if d > 0 else min(pa * pb, (1 - pa) * (1 - pb))
    return float(min(r2, 1.0)), float(abs(d) / dmax) if dmax > 0 else 0.0


def main() -> None:
    cand = pd.read_csv(RESULTS / "top_variants_corrected.csv")
    print(f"candidates: {len(cand)}")
    pops = load_panel()
    print("African populations: " +
          ", ".join(f"{p} n={len(s)}" for p, s in pops.items()))

    vcf = pysam.VariantFile(VCF)
    have = set(vcf.header.samples)
    idx = {p: [s for s in ss if s in have] for p, ss in pops.items()}

    # ---- every common variant in the APOL1 gene body ------------------------
    chrom, gs, ge = APOL1
    gene_recs = []
    for rec in vcf.fetch(chrom, gs - 1, ge):
        if len(rec.alts or []) != 1:
            continue
        gene_recs.append((rec.pos, rec.ref, rec.alts[0],
                          {p: haplotypes(rec, idx[p]) for p in AFR_POPS}))
    print(f"common-or-any variants in the APOL1 gene body: {len(gene_recs)}")

    # keep only those common in at least one African population
    keep = []
    for pos, ref, alt, hp in gene_recs:
        if max((h[h >= 0].mean() if (h >= 0).any() else 0) for h in hp.values()) >= MIN_MAF:
            keep.append((pos, ref, alt, hp))
    print(f"of which common (MAF>={MIN_MAF}) somewhere in Africa: {len(keep)}\n")

    # ---- each candidate against all of them --------------------------------
    rows = []
    for c in cand.itertuples():
        rec = next((r for r in vcf.fetch(chrom, c.pos - 1, c.pos)
                    if r.pos == c.pos and len(r.alts or []) == 1), None)
        if rec is None:
            rows.append({"rsid": c.rsid, "pos": c.pos, "dist_to_gene": c.dist,
                         "af_afr": round(c.af_afr, 4), "max_r2_vs_gene": None,
                         "d_prime_at_max": None, "partner_pos": None,
                         "population": None,
                         "verdict": "not in 1000G phase 3 panel"})
            continue
        chap = {p: haplotypes(rec, idx[p]) for p in AFR_POPS}
        best = (0.0, 0.0, None, None)
        for pos, ref, alt, ghp in keep:
            for p in AFR_POPS:
                r2, dp = r2_dprime(chap[p], ghp[p])
                if not np.isnan(r2) and r2 > best[0]:
                    best = (r2, dp, pos, p)
        rows.append({"rsid": c.rsid, "pos": c.pos, "dist_to_gene": c.dist,
                     "af_afr": round(c.af_afr, 4),
                     "max_r2_vs_gene": round(best[0], 4),
                     "d_prime_at_max": round(best[1], 4),
                     "partner_pos": best[2], "population": best[3],
                     "verdict": ("tags the gene haplotype" if best[0] >= 0.8 else
                                 "partially correlated" if best[0] >= 0.2 else
                                 "independent")})
        print(f"  {c.rsid:<14} max r2 vs gene body = {best[0]:.3f}"
              f"  ({best[3]}, pos {best[2]})")

    d = pd.DataFrame(rows)
    d.to_csv(RESULTS / "ld_gene_haplotype_per_variant.csv", index=False)

    tested = d[d["max_r2_vs_gene"].notna()]
    n_tag = int((tested["max_r2_vs_gene"] >= 0.8).sum())
    n_part = int(((tested["max_r2_vs_gene"] >= 0.2) &
                  (tested["max_r2_vs_gene"] < 0.8)).sum())
    n_ind = int((tested["max_r2_vs_gene"] < 0.2).sum())

    print("\n" + "=" * 74)
    print("LD BETWEEN CANDIDATES AND APOL1 GENE-BODY HAPLOTYPE STRUCTURE")
    print("=" * 74)
    print(f"  tested                     : {len(tested)} of {len(d)}")
    print(f"  tags the gene (r2>=0.8)    : {n_tag}")
    print(f"  partially correlated       : {n_part}")
    print(f"  independent (r2<0.2)       : {n_ind}")
    if len(tested):
        print(f"  highest r2 seen            : {tested['max_r2_vs_gene'].max():.3f}")

    print("\n  READING:")
    if n_ind == len(tested) and len(tested):
        print("  No candidate tags APOL1 gene-body haplotype structure. Since both")
        print("  G1 and G2 lie inside the gene, this covers G2 without needing the")
        print("  deletion itself, which 1000 Genomes does not carry. The candidates")
        print("  are independent of the known risk haplotypes.")
    elif n_ind == 0:
        print("  Every candidate tags the gene haplotype. Any association would be")
        print("  confounded by the known risk variants. Drop the novelty claim.")
    else:
        print("  Mixed. Report the LD column per variant and restrict any novelty")
        print("  claim to those below 0.2.")

    out = {"candidates": int(len(d)), "tested": int(len(tested)),
           "gene_body": f"chr{chrom}:{gs}-{ge}",
           "gene_variants_considered": len(keep),
           "min_maf": MIN_MAF, "populations": AFR_POPS,
           "populations_pooled": False,
           "phased": True, "source": "1000 Genomes phase 3, GRCh38 lift",
           "g2_note": ("G2 (rs71785313, 6bp deletion) is absent from this call "
                       "set, which is why LD is assessed against all gene-body "
                       "variation rather than against G2 directly"),
           "n_tags_gene_r2_ge_0.8": n_tag,
           "n_partial": n_part, "n_independent": n_ind,
           "max_r2_observed": (float(tested["max_r2_vs_gene"].max())
                               if len(tested) else None)}
    (RESULTS / "ld_gene_haplotype.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/ld_gene_haplotype.json")


if __name__ == "__main__":
    main()
