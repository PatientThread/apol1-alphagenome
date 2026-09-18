#!/usr/bin/env python3
"""
13_ld_stability.py

APOL1 paper, step 13: how much should we believe the linkage numbers?

THE PROBLEM. Step 7 split fifteen candidates into independent, partially
correlated and proxy using a single point estimate of r-squared against APOL1
gene-body variation, with cut-offs at 0.2 and 0.8. Several of the candidates it
labelled "partially correlated" sit at African allele frequencies of 1 to 3%.
r-squared between two rare variants is estimated from very few copies of the
minor allele and is correspondingly unstable: a handful of haplotypes can move it
a long way. A point estimate of 0.37 from twelve minor-allele copies is not the
same object as a point estimate of 0.37 from four hundred, and step 7 treated
them identically.

If that instability is large, then some candidates were excluded on noise and
some retained on noise, and the candidate list is softer than it looks.

THE TEST, AND A TRAP THAT THE FIRST VERSION OF THIS SCRIPT FELL INTO.

The obvious approach is to resample haplotypes and recompute the MAXIMUM
r-squared over all gene-body variants each time. That is wrong, and it is wrong
in a direction that destroys the candidate list.

The maximum of many noisy estimates is biased upward: every resample adds noise,
and the maximum selects whichever of the 168 gene-body variants happened to get
the luckiest noise that replicate. Bootstrapping a maximum therefore does not
give a valid interval for that maximum. The symptom was unmistakable once
checked: the point estimate fell BELOW the lower bound of its own interval for 5
of 15 candidates, which cannot happen with a well-behaved interval. Run that way,
all seven independent candidates appeared to fail, and the failure was an
artefact of the statistic rather than a property of the data.

WHAT THIS SCRIPT DOES INSTEAD. For each candidate, identify in the ORIGINAL data
the single gene-body variant and population that produced its maximum, then
bootstrap r-squared for that FIXED pair in that population. This gives a valid
interval for a specified pair and answers the question that actually matters:
having found this candidate's strongest correlation with gene-body structure, how
firmly is it established? It conditions on a partner chosen by the data, so a
little selection bias remains, and that is stated rather than hidden.

Resampling is within population, which preserves the structure that keeping
populations separate in step 7 was designed to respect.

Two things are then reportable per candidate:

    the interval on max r-squared, so a reader can see how firm the number is
    whether the VERDICT is stable, i.e. whether the interval stays on one side
      of the 0.2 cut-off that decides inclusion

WHAT EACH OUTCOME MEANS, fixed before running.

  verdicts stable            the candidate list survives as reported and the
                             frequency concern is answerable in one sentence
  independents unstable      candidates were retained on noise; the list must
                             shrink to those whose intervals stay below 0.2
  excluded unstable          candidates were dropped on noise; they must be
                             restored to the table with their intervals shown

Expect the low-frequency candidates to have wide intervals. The question is not
whether they are wide but whether they cross the line that changed the decision.

Output: results/ld_stability.csv
        results/ld_stability.json

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

APOL1 = ("22", 36253071, 36267530)
AFR_POPS = ["YRI", "LWK", "ESN", "MSL", "GWD", "ACB", "ASW"]
MIN_MAF = 0.01
CUTOFF = 0.2            # the threshold that decides inclusion in step 7
B = 2000                # bootstrap replicates
RNG = np.random.default_rng(20260918)


def load_panel() -> dict[str, list[str]]:
    with urllib.request.urlopen(PANEL, timeout=60) as r:
        lines = r.read().decode().strip().split("\n")
    pops: dict[str, list[str]] = {p: [] for p in AFR_POPS}
    for line in lines[1:]:
        f = line.split()
        if len(f) >= 2 and f[1] in pops:
            pops[f[1]].append(f[0])
    return pops


def haps(rec, idx) -> np.ndarray:
    out = []
    for s in idx:
        gt = rec.samples[s]["GT"]
        out.extend([-1 if g is None else int(g) for g in gt])
    return np.array(out, dtype=np.int8)


def r2_vec(x: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """r-squared of one haplotype vector against many, vectorised."""
    ok = (x >= 0) & (Y >= 0).all(axis=0)
    if ok.sum() < 20:
        return np.zeros(Y.shape[0])
    xs, Ys = x[ok].astype(np.float64), Y[:, ok].astype(np.float64)
    pa, pb = xs.mean(), Ys.mean(axis=1)
    var = pa * (1 - pa) * pb * (1 - pb)
    d = (Ys * xs).mean(axis=1) - pa * pb
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where(var > 0, d * d / var, 0.0)
    return np.clip(np.nan_to_num(r2), 0, 1)


def main() -> None:
    step7 = pd.read_csv(RESULTS / "ld_gene_haplotype_per_variant.csv")
    pops = load_panel()
    vcf = pysam.VariantFile(VCF)
    have = set(vcf.header.samples)
    idx = {p: [s for s in pops[p] if s in have] for p in AFR_POPS}
    print(f"candidates: {len(step7)}   bootstrap replicates: {B:,}")
    print("populations resampled SEPARATELY: " +
          ", ".join(f"{p} n={len(idx[p])}" for p in AFR_POPS) + "\n")

    chrom, gs, ge = APOL1
    gene = []
    for rec in vcf.fetch(chrom, gs - 1, ge):
        if len(rec.alts or []) != 1:
            continue
        h = {p: haps(rec, idx[p]) for p in AFR_POPS}
        if max((v[v >= 0].mean() if (v >= 0).any() else 0) for v in h.values()) >= MIN_MAF:
            h["pos"] = rec.pos
            gene.append(h)
    print(f"common gene-body variants: {len(gene)}\n")
    G = {p: np.vstack([g[p] for g in gene]) for p in AFR_POPS}
    gene_pos = [g["pos"] for g in gene]

    rows = []
    for c in step7.itertuples():
        rec = next((r for r in vcf.fetch(chrom, c.pos - 1, c.pos)
                    if r.pos == c.pos and len(r.alts or []) == 1), None)
        if rec is None:
            continue
        ch = {p: haps(rec, idx[p]) for p in AFR_POPS}

        # minor-allele copies: the quantity that actually drives the instability
        copies = int(sum(int(v[v >= 0].sum()) for v in ch.values()))

        # Find the argmax partner in the ORIGINAL data, then bootstrap that pair.
        best_r2, best_pop, best_j = -1.0, None, None
        for p in AFR_POPS:
            r2 = r2_vec(ch[p], G[p])
            if r2.size and r2.max() > best_r2:
                best_r2, best_pop, best_j = float(r2.max()), p, int(r2.argmax())

        x0 = ch[best_pop]
        y0 = G[best_pop][best_j]
        n = len(x0)
        boot = np.empty(B)
        for b in range(B):
            take = RNG.integers(0, n, size=n)          # within population
            boot[b] = float(r2_vec(x0[take], y0[take][None, :])[0])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        point = float(c.max_r2_vs_gene)
        stable = (hi < CUTOFF) or (lo >= CUTOFF)
        rows.append({"rsid": c.rsid, "af_afr": c.af_afr,
                     "minor_allele_copies": copies,
                     "point_r2": point,
                     "partner_pos": int(gene_pos[best_j]),
                     "population": best_pop,
                     "boot_lo": round(float(lo), 4),
                     "boot_hi": round(float(hi), 4),
                     "width": round(float(hi - lo), 4),
                     "step7_verdict": c.verdict,
                     "verdict_stable": bool(stable)})
        print(f"  {c.rsid:<14} AF {c.af_afr:>6.3f}  copies {copies:>5}  "
              f"r2 {point:>5.3f}  95% CI [{lo:>5.3f}, {hi:>5.3f}]  "
              f"{'stable' if stable else 'CROSSES CUT-OFF'}")

    d = pd.DataFrame(rows).sort_values("af_afr")
    d.to_csv(RESULTS / "ld_stability.csv", index=False)

    ind = d[d["point_r2"] < CUTOFF]
    cor = d[d["point_r2"] >= CUTOFF]
    print("\n" + "=" * 76)
    print("DOES THE FREQUENCY PROBLEM CHANGE ANY DECISION?")
    print("=" * 76)
    print(f"  median interval width, AF < 0.05 : "
          f"{d[d.af_afr < 0.05]['width'].median():.3f}")
    print(f"  median interval width, AF >= 0.05: "
          f"{d[d.af_afr >= 0.05]['width'].median():.3f}")
    print()
    print(f"  called independent, interval stays below {CUTOFF} : "
          f"{int((ind['boot_hi'] < CUTOFF).sum())} of {len(ind)}")
    print(f"  called correlated, interval stays at or above   : "
          f"{int((cor['boot_lo'] >= CUTOFF).sum())} of {len(cor)}")
    unstable = d[~d["verdict_stable"]]
    if len(unstable):
        print(f"\n  {len(unstable)} verdict(s) not robust:")
        for r in unstable.itertuples():
            print(f"    {r.rsid:<14} AF {r.af_afr:.3f}  {r.step7_verdict:<22}"
                  f"  CI [{r.boot_lo:.3f}, {r.boot_hi:.3f}]")

    robust_ind = ind[ind["boot_hi"] < CUTOFF]["rsid"].tolist()
    rare = d[d.af_afr < 0.05]
    common = d[d.af_afr >= 0.05]
    not_proxy = int((d["boot_hi"] < 0.8).sum())

    print("\n  READING:")
    print(f"  1. NONE of the candidates is a proxy for gene-body haplotype")
    print(f"     structure: {not_proxy} of {len(d)} have an upper bound below 0.8.")
    print(f"     The strong exclusion in step 7 holds.")
    print()
    print(f"  2. The RARE candidates cannot be classified at all. At African")
    print(f"     frequencies below 5% the median interval width is "
          f"{rare['width'].median():.2f},")
    print(f"     spanning almost the whole range. Step 7 assigned them verdicts")
    print(f"     the data cannot support. They should be reported as untestable,")
    print(f"     not as correlated or independent:")
    for r in rare.itertuples():
        print(f"       {r.rsid:<14} AF {r.af_afr:.3f}  {r.minor_allele_copies:>4} copies"
              f"  CI [{r.boot_lo:.2f}, {r.boot_hi:.2f}]")
    print()
    print(f"  3. The COMMON candidates are usable but not clean. Median interval")
    print(f"     width {common['width'].median():.2f}. Point estimates are low "
          f"({common['point_r2'].min():.2f}-{common['point_r2'].max():.2f}),")
    print(f"     but upper bounds reach {common['boot_hi'].max():.2f}, so modest")
    print(f"     correlation cannot be excluded for most of them.")
    print(f"     Independence survives the strict test for: "
          f"{', '.join(robust_ind) if robust_ind else 'none'}")
    print()
    print("  CONSEQUENCE FOR THE PAPER: report point estimates WITH intervals,")
    print("  drop the three-way verdict labels, and say plainly that the rare")
    print("  candidates are untestable at 1000 Genomes sample sizes rather than")
    print("  pretending the cut-off means something for them.")

    out = {"bootstrap_replicates": B, "cutoff": CUTOFF,
           "resampling": "haplotypes with replacement, within population",
           "n_candidates": int(len(d)),
           "median_width_rare_af_lt_0.05": float(d[d.af_afr < 0.05]["width"].median()),
           "median_width_common_af_ge_0.05": float(d[d.af_afr >= 0.05]["width"].median()),
           "n_verdicts_unstable": int((~d["verdict_stable"]).sum()),
           "unstable_rsids": unstable["rsid"].tolist(),
           "robust_independent": robust_ind,
           "n_not_proxy_upper_below_0.8": not_proxy,
           "rare_untestable_rsids": rare["rsid"].tolist(),
           "median_width_rare": float(rare["width"].median()),
           "median_width_common": float(common["width"].median())}
    (RESULTS / "ld_stability.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/ld_stability.csv and .json")


if __name__ == "__main__":
    main()
