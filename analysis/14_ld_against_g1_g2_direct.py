#!/usr/bin/env python3
"""
14_ld_against_g1_g2_direct.py

APOL1 paper, step 14: the G1/G2 analysis that steps 6 and 7 said was impossible.

WHAT WAS CLAIMED, AND WHY IT WAS WRONG. Step 6 reported that G2 could not be
tested because Ensembl excludes indels from its linkage service. Step 7 went
further and stated that the six-base deletion is "absent from the 1000 Genomes
phase 3 call set entirely", and substituted a test against all common gene-body
variation. Editorial review challenged that, citing a published nomenclature
clarification, and the challenge is correct.

G2 is absent from the file step 7 used, ALL.chr22.shapeit2_integrated_snvindels
_v2a_27022019.GRCh38.phased.vcf.gz, which is a biallelic-filtered GRCh38 lift. It
is PRESENT in the original phase 3 release on GRCh37, at chr22:36,662,041 as
AATAATT -> A, carrying no rsID, which is why every identifier-based search
failed. Its allele frequency confirms the identification beyond doubt: 19.5% in
Gambian, 18.2% in Mende, 11.6% in Esan, 13.0% in African Caribbean, and 0% in
CEU, GBR and CHB. That is the G2 distribution and nothing else looks like it.

So the limitation was a property of the file chosen, not of the resource. This
script does the analysis properly.

THREE THINGS IT COMPUTES, and the third is the one that matters most.

  1. r-squared against G1 and against G2 SEPARATELY, per African population,
     from phased haplotypes.

  2. D' alongside r-squared. Review made the point that these answer different
     questions: D' near 1 with low r-squared means one allele sits entirely on
     the other's background but at a different frequency, which is haplotypic
     dependence that r-squared alone hides.

  3. THE MAXIMUM r-squared ATTAINABLE given the two allele frequencies. This is
     the argument that most damages the earlier analysis. With frequencies of
     0.026 and 0.20, r-squared cannot exceed about 0.107 even under complete
     nesting, so observing 0.115 and calling it independence is meaningless. For
     each pair we report observed r-squared as a FRACTION of its attainable
     maximum, which is the only version of the statistic that can support a
     claim of independence.

Output: results/ld_g1_g2_direct.csv, results/ld_g1_g2_direct.json

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

# The ORIGINAL phase 3 release on GRCh37, which retains the G2 deletion.
VCF = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
       "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
PANEL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
         "integrated_call_samples_v3.20130502.ALL.panel")

AFR_POPS = ["YRI", "LWK", "ESN", "MSL", "GWD", "ACB", "ASW"]

# GRCh37 positions, from Ensembl's GRCh37 service. G1 verified against the
# published coordinate; G2 identified by position and alleles, not by rsID.
RISK = {
    "G1 rs73885319": (36_661_906, None),
    "G1 rs60910145": (36_662_034, None),
    "G2 6bp del":    (36_662_041, ("AATAATT", "A")),
}
CANDIDATES = {
    "rs136204": 36_754_161, "rs183925240": 36_754_162, "rs6000250": 36_746_801,
    "rs4820232": 36_710_541, "rs713797": 36_764_534, "rs132708": 36_593_197,
    "rs5750234": 36_580_879, "rs132754": 36_616_891, "rs114518164": 36_576_518,
    "rs132695": 36_580_965, "rs150645785": 36_758_565, "rs77773194": 36_707_934,
    "rs78440576": 36_616_761, "rs8138011": 36_733_423, "rs148651564": 36_746_984,
}


def load_panel() -> dict[str, list[str]]:
    lines = urllib.request.urlopen(PANEL, timeout=90).read().decode().strip().split("\n")
    pops: dict[str, list[str]] = {p: [] for p in AFR_POPS}
    for line in lines[1:]:
        f = line.split()
        if len(f) >= 2 and f[1] in pops:
            pops[f[1]].append(f[0])
    return pops


def haps(rec, samples) -> np.ndarray:
    out = []
    for s in samples:
        gt = rec.samples[s]["GT"]
        out.extend([-1 if g is None else (1 if g else 0) for g in gt])
    return np.array(out, dtype=np.int8)


def ld_stats(x: np.ndarray, y: np.ndarray) -> dict:
    """r-squared, D', and the maximum r-squared the frequencies permit."""
    ok = (x >= 0) & (y >= 0)
    x, y = x[ok].astype(float), y[ok].astype(float)
    if len(x) < 20:
        return {}
    pa, pb = x.mean(), y.mean()
    if min(pa, pb) <= 0 or max(pa, pb) >= 1:
        return {"r2": 0.0, "dprime": 0.0, "r2_max": 0.0, "frac_of_max": 0.0,
                "pa": pa, "pb": pb, "n_hap": len(x)}
    d = (x * y).mean() - pa * pb
    denom = pa * (1 - pa) * pb * (1 - pb)
    r2 = min(d * d / denom, 1.0)
    dmax = min(pa * (1 - pb), (1 - pa) * pb) if d > 0 else min(pa * pb, (1 - pa) * (1 - pb))
    dprime = abs(d) / dmax if dmax > 0 else 0.0
    # Maximum attainable r-squared: D at its own maximum, same frequencies.
    dmax_pos = min(pa * (1 - pb), (1 - pa) * pb)
    r2_max = min(dmax_pos ** 2 / denom, 1.0)
    return {"r2": r2, "dprime": dprime, "r2_max": r2_max,
            "frac_of_max": r2 / r2_max if r2_max > 0 else np.nan,
            "pa": pa, "pb": pb, "n_hap": len(x)}


def main() -> None:
    pops = load_panel()
    vcf = pysam.VariantFile(VCF)
    have = set(vcf.header.samples)
    idx = {p: [s for s in pops[p] if s in have] for p in AFR_POPS}
    print("populations:", ", ".join(f"{p} n={len(idx[p])}" for p in AFR_POPS), "\n")

    def get(pos, alleles=None):
        for rec in vcf.fetch("22", pos - 2, pos + 2):
            if rec.pos != pos:
                continue
            if alleles and not (rec.ref == alleles[0] and alleles[1] in (rec.alts or [])):
                continue
            return rec
        return None

    risk_recs = {}
    for name, (pos, alleles) in RISK.items():
        r = get(pos, alleles)
        if r is None:
            raise SystemExit(f"{name} not found at {pos}")
        risk_recs[name] = {p: haps(r, idx[p]) for p in AFR_POPS}
        af = {p: round(float(v[v >= 0].mean()), 4) for p, v in risk_recs[name].items()}
        print(f"  {name:<16} pos {pos}  {r.ref}->{','.join(r.alts or [])}  AF {af}")

    rows = []
    for rs, pos in CANDIDATES.items():
        rec = get(pos)
        if rec is None:
            rows.append({"rsid": rs, "risk_variant": None,
                         "note": "not found in phase 3 GRCh37"})
            continue
        ch = {p: haps(rec, idx[p]) for p in AFR_POPS}
        for name in RISK:
            best = None
            for p in AFR_POPS:
                s = ld_stats(ch[p], risk_recs[name][p])
                if not s:
                    continue
                if best is None or s["r2"] > best["r2"]:
                    best = {**s, "population": p}
            if best:
                rows.append({"rsid": rs, "risk_variant": name,
                             "population": best["population"],
                             "af_candidate": round(best["pa"], 4),
                             "af_risk": round(best["pb"], 4),
                             "r2": round(best["r2"], 4),
                             "d_prime": round(best["dprime"], 4),
                             "r2_max_possible": round(best["r2_max"], 4),
                             "frac_of_max": round(best["frac_of_max"], 3)})

    d = pd.DataFrame(rows)
    d.to_csv(RESULTS / "ld_g1_g2_direct.csv", index=False)

    print("\n" + "=" * 92)
    print("LINKAGE WITH G1 AND G2 DIRECTLY, maximum across seven African populations")
    print("=" * 92)
    print(f"  {'variant':<13}{'vs':<16}{'r2':>7}{'D-prime':>9}"
          f"{'r2 max':>8}{'% of max':>10}  population")
    for r in d.dropna(subset=["risk_variant"]).itertuples():
        print(f"  {r.rsid:<13}{r.risk_variant:<16}{r.r2:>7.3f}{r.d_prime:>9.3f}"
              f"{r.r2_max_possible:>8.3f}{100*r.frac_of_max:>9.0f}%  {r.population}")

    g2 = d[d["risk_variant"] == "G2 6bp del"]
    hi_dp = d[(d["d_prime"] >= 0.8) & (d["r2"] < 0.2)]
    print("\n  READING:")
    print(f"  G2 was testable after all: {len(g2)} candidate comparisons computed.")
    print(f"  Highest r-squared against any risk variant: {d['r2'].max():.3f}")
    print(f"  Highest as a fraction of what the frequencies permit: "
          f"{d['frac_of_max'].max():.2f}")
    if len(hi_dp):
        print(f"\n  {len(hi_dp)} comparison(s) have D' >= 0.8 with r-squared < 0.2.")
        print("  Low r-squared there reflects a frequency difference, NOT haplotypic")
        print("  independence, and those candidates cannot be called independent:")
        for r in hi_dp.itertuples():
            print(f"    {r.rsid:<13} vs {r.risk_variant:<16} r2 {r.r2:.3f} "
                  f"D' {r.d_prime:.3f}  (max r2 possible {r.r2_max_possible:.3f})")
    else:
        print("\n  No comparison combines high D' with low r-squared.")

    out = {"source": "1000 Genomes phase 3 original release, GRCh37",
           "why_not_grch38_file": ("the biallelic GRCh38 lift used in step 7 omits "
                                   "the G2 deletion; the original release retains it "
                                   "at chr22:36,662,041 AATAATT>A with no rsID"),
           "populations": AFR_POPS, "populations_pooled": False,
           "comparisons": int(len(d.dropna(subset=["risk_variant"]))),
           "max_r2": float(d["r2"].max()),
           "max_frac_of_attainable": float(d["frac_of_max"].max()),
           "n_high_dprime_low_r2": int(len(hi_dp))}
    (RESULTS / "ld_g1_g2_direct.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/ld_g1_g2_direct.csv and .json")


if __name__ == "__main__":
    main()
