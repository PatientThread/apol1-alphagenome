#!/usr/bin/env python3
"""
06_ld_with_risk_haplotypes.py

APOL1 paper, step 6: the analysis that decides whether there is a paper.

THE QUESTION. Steps 4 and 5 produced a ranked list of non-coding candidates and
showed that the top of that list falls inside independently measured kidney open
chromatin far more often than the other common variants at this locus. That is a
real result about the model. It is NOT yet a result about APOL1 biology, because
of one thing that has not been tested.

G1 and G2 sit on extended risk haplotypes that reach well beyond the coding
change. Any variant in linkage disequilibrium with them will track kidney disease
risk in an association study for reasons that have nothing to do with its own
regulatory function. If our candidates are simply in LD with G1 or G2, we have
rediscovered the risk haplotype in a new coordinate system and dressed it up as
a set of novel regulatory variants. Script 03 flagged this explicitly and
deferred it, noting it needs phased reference haplotypes rather than the allele
frequencies used there. This script closes that gap.

    H0  top-ranked candidates are no more correlated with G1/G2 than the other
        common variants at this locus
    H1  they are enriched for LD with the risk haplotypes, i.e. confounded

Note the direction. Here the NULL is the good outcome. We want the candidates to
be independent of the risk haplotype, because only then are they new.

THE DESIGN. r-squared and D' are taken from 1000 Genomes phase 3 for every
African population separately, not from a pooled AFR superpopulation. Pooling
African populations inflates LD through population structure, which is exactly
the artefact that would manufacture a false positive here. A candidate is judged
against its MAXIMUM r-squared across populations, which is the conservative
choice: one population showing correlation is enough to make a candidate
suspect.

    r2 >= 0.8   effectively a proxy for the risk haplotype. Not independent.
    0.2 - 0.8   partially correlated. Reportable but cannot be called novel.
    r2 <  0.2   independent of G1/G2 by the usual convention.

D' is reported alongside because high D' with low r-squared means the variants
sit on a shared haplotype background but differ in frequency. That still matters
for interpretation even when r-squared is low, and reporting r-squared alone
would hide it.

WHAT EACH OUTCOME MEANS, decided before running.

  candidates independent      the ranking nominates variants that are NOT
                              explained by the known risk haplotype, and the
                              paper can report them as novel candidates
  candidates in high LD       the ranking recovers the risk haplotype; the
                              "novel candidate" framing must be dropped and the
                              paper becomes a methods caution only
  mixed                       report per-variant, name which survive, and make
                              the LD column part of the main results table

Output: results/ld_with_risk.json
        results/ld_per_variant.csv

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
RESULTS = ROOT / "results"
CACHE = RESULTS / "ld_cache"
CACHE.mkdir(parents=True, exist_ok=True)

SERVER = "https://rest.ensembl.org"

# The risk alleles. Positions verified against Ensembl on GRCh38.
RISK = {
    "rs73885319": "G1 p.S342G  chr22:36,265,860",
    "rs60910145": "G1 p.I384M  chr22:36,265,988",
    "rs71785313": "G2 6bp del  chr22:36,265,996",
}

# African populations, kept SEPARATE. Pooling them inflates LD via structure.
POPS = [
    "1000GENOMES:phase_3:YRI",   # Yoruba, Nigeria
    "1000GENOMES:phase_3:LWK",   # Luhya, Kenya
    "1000GENOMES:phase_3:ESN",   # Esan, Nigeria
    "1000GENOMES:phase_3:MSL",   # Mende, Sierra Leone
    "1000GENOMES:phase_3:GWD",   # Gambian
    "1000GENOMES:phase_3:ACB",   # African Caribbean, Barbados
    "1000GENOMES:phase_3:ASW",   # African ancestry, southwest USA
]

WINDOW_KB = 500          # Ensembl maximum; the locus spans 248 kb
R2_FLOOR = 0.01          # fetch low so we can report true independence, not absence


def fetch_ld(rsid: str, pop: str) -> list[dict]:
    """Every variant in LD with `rsid` in `pop`. Cached, so a re-run is free."""
    key = CACHE / f"{rsid}__{pop.replace(':', '_')}.json"
    if key.exists():
        return json.loads(key.read_text())
    url = (f"{SERVER}/ld/human/{rsid}/{pop}"
           f"?window_size={WINDOW_KB};r2={R2_FLOOR};content-type=application/json")
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = json.loads(r.read().decode())
            key.write_text(json.dumps(data))
            time.sleep(0.4)                      # be polite to Ensembl
            return data
        except urllib.error.HTTPError as e:
            if e.code == 429:                    # rate limited
                time.sleep(int(e.headers.get("Retry-After", 5)) + 1)
                continue
            if e.code in (400, 404):             # variant absent from this panel
                key.write_text("[]")
                return []
            raise
        except Exception:                        # noqa: BLE001
            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"Ensembl failed repeatedly for {rsid} in {pop}")


def main() -> None:
    top = pd.read_csv(RESULTS / "top_variants_corrected.csv")
    print(f"candidates to test: {len(top)}")
    print(f"risk variants: {', '.join(RISK)}")
    print(f"populations: {len(POPS)}, kept separate\n")

    # Build risk-variant -> population -> {partner: (r2, dprime)}
    ld: dict[str, dict[str, dict[str, tuple[float, float]]]] = {}
    control = {}
    for rs in RISK:
        ld[rs] = {}
        for pop in POPS:
            rows = fetch_ld(rs, pop)
            ld[rs][pop] = {r["variation2"]: (float(r["r2"]), float(r["d_prime"]))
                           for r in rows}
            print(f"  {rs:<12} {pop.split(':')[-1]:<4} "
                  f"{len(rows):>5} partners at r2>={R2_FLOOR}")
        # sanity check: the two G1 SNPs must be in near-perfect LD with each other
        if rs == "rs73885319":
            for pop in POPS:
                v = ld[rs][pop].get("rs60910145")
                if v:
                    control[pop.split(":")[-1]] = v[0]

    print("\nCONTROL, the two G1 SNPs against each other (must be near 1.0):")
    print("  " + "  ".join(f"{k} {v:.3f}" for k, v in control.items()))
    if control and min(control.values()) < 0.8:
        raise SystemExit("Control failed: the G1 SNPs are not in LD with each "
                         "other. The LD source is wrong; do not use these numbers.")

    # ---------------------------------------------------------------- per variant
    out_rows = []
    for t in top.itertuples():
        rec = {"rsid": t.rsid, "pos": t.pos, "dist_to_gene": t.dist,
               "af_afr": round(t.af_afr, 4),
               "podocyte_dnase": t.dnase_glomerular_visceral_epithelial_cell}
        best_r2, best_dp, best_where = 0.0, 0.0, ""
        for rs in RISK:
            for pop in POPS:
                hit = ld[rs][pop].get(t.rsid)
                if hit and hit[0] > best_r2:
                    best_r2, best_dp = hit
                    best_where = f"{rs}/{pop.split(':')[-1]}"
        rec["max_r2"] = round(best_r2, 4)
        rec["d_prime_at_max"] = round(best_dp, 4)
        rec["max_against"] = best_where or "none above floor"
        rec["verdict"] = ("proxy for risk haplotype" if best_r2 >= 0.8 else
                          "partially correlated" if best_r2 >= 0.2 else
                          "independent")
        out_rows.append(rec)

    d = pd.DataFrame(out_rows).sort_values("max_r2", ascending=False)
    d.to_csv(RESULTS / "ld_per_variant.csv", index=False)

    print("\n" + "=" * 78)
    print("LD BETWEEN TOP CANDIDATES AND THE G1/G2 RISK HAPLOTYPES")
    print("=" * 78)
    print(f"  {'rsid':<14}{'dist':>8}{'AF_afr':>8}{'max r2':>8}{'D-prime':>9}  verdict")
    for r in d.itertuples():
        print(f"  {r.rsid:<14}{r.dist_to_gene:>8,}{r.af_afr:>8.3f}"
              f"{r.max_r2:>8.3f}{r.d_prime_at_max:>9.3f}  {r.verdict}")

    n_proxy = int((d["max_r2"] >= 0.8).sum())
    n_part = int(((d["max_r2"] >= 0.2) & (d["max_r2"] < 0.8)).sum())
    n_indep = int((d["max_r2"] < 0.2).sum())

    print()
    print(f"  proxy for the risk haplotype (r2>=0.8) : {n_proxy} of {len(d)}")
    print(f"  partially correlated (0.2-0.8)         : {n_part} of {len(d)}")
    print(f"  independent (r2<0.2)                   : {n_indep} of {len(d)}")

    print("\n  READING:")
    if n_indep == len(d):
        print("  Every candidate is independent of G1 and G2. The ranking is not")
        print("  rediscovering the risk haplotype, and these can be reported as")
        print("  novel regulatory candidates.")
    elif n_indep == 0:
        print("  No candidate is independent. The ranking recovers the risk")
        print("  haplotype. The novel-candidate framing must be DROPPED and the")
        print("  paper reduced to the methodological finding.")
    else:
        print(f"  Mixed. {n_indep} of {len(d)} are independent and can be reported")
        print("  as candidates; the rest must be labelled as correlated with the")
        print("  risk haplotype. The LD column belongs in the main results table.")

    out = {
        "candidates_tested": int(len(d)),
        "risk_variants": RISK,
        "populations": POPS,
        "populations_pooled": False,
        "window_kb": WINDOW_KB,
        "r2_floor_fetched": R2_FLOOR,
        "control_g1_internal_r2": control,
        "n_proxy_r2_ge_0.8": n_proxy,
        "n_partial_0.2_to_0.8": n_part,
        "n_independent_r2_lt_0.2": n_indep,
        "max_r2_observed": float(d["max_r2"].max()),
        "median_r2": float(d["max_r2"].median()),
    }
    (RESULTS / "ld_with_risk.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/ld_with_risk.json and ld_per_variant.csv")


if __name__ == "__main__":
    main()
