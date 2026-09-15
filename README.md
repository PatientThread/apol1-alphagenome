# Regulatory variant prioritisation at the APOL1 locus

Analysis code and data for a study in preparation.

**Author:** Christopher Lawrence, Consultant Nephrologist, London, United Kingdom.
ORCID [0000-0002-8159-0879](https://orcid.org/0000-0002-8159-0879).

---

## The question

G1 and G2 explain a large part of the excess kidney disease risk carried on
APOL1, and both are coding. The locus is overwhelmingly non-coding: 98,875 of
107,737 catalogued variants here, 91.8%. Whether any of that non-coding
variation contributes regulatory effects is unresolved. This study asks whether
a sequence-to-function model can prioritise candidates, and then tests the
answer against data the model never saw.

## What is established so far

1. **The obvious readout is confounded and is not used.** Ranking by the
   gene-masked expression scorer largely rediscovers proximity to APOL1.
   Non-coding variants inside the gene body outscore coding ones, so this is
   position rather than constraint. The pre-registered negative control, the
   coding risk variants, stayed at the 97th and 98th percentile even after
   residualising on distance. That channel is reported as confounded and is not
   used for the ranking.
2. **The accessibility channel is not confounded.** Distance explains 0.3% of
   the variance for podocyte DNase and essentially none for kidney ATAC. None of
   the top 15 candidates sits inside the gene body, against 6% expected.
3. **Its top-ranked variants are enriched in independently measured kidney open
   chromatin.** Against a background of the other common variants at the same
   locus, not the genome: odds ratio 8.1 for the top 25 by kidney ATAC
   (p = 2e-6), 6.1 for the top 200 by podocyte DNase (p = 1e-22). The
   enrichment survives a stricter peak definition.
4. **Half the candidates are correlated with APOL1 haplotype structure and
   cannot be called novel.** Tested against every common gene-body variant in
   seven African populations, kept separate, from phased haplotypes: 7 of 15 are
   independent (r² < 0.2), 7 are partially correlated, and one is a perfect
   proxy. Six of the seven independent ones are common enough for the estimate
   to be stable, and those are the defensible candidate set.

## What this repository does *not* claim

No candidate has been functionally validated. Accessibility enrichment shows the
ranking selects positions that are open in kidney; it does not show that any
variant changes expression, and no disease association is tested here.

---

## Reproducing the analysis

Scripts are numbered in dependency order.

| Script | Needs network | Needs API key | Produces |
|---|---|---|---|
| `01_fetch_locus_variants.py` | yes (gnomAD) | no | the frozen locus variant set |
| `02_score_regulatory_effects.py` | yes | **yes** | per-variant predictions |
| `03_analyse.py` | no | no | distributions, first ranking |
| `04_distance_controlled.py` | no | no | **the proximity correction** |
| `05_validate_against_encode.py` | yes (ENCODE) | no | **the independent validation** |
| `06_ld_with_risk_haplotypes.py` | yes (Ensembl) | no | linkage against G1 |
| `07_ld_against_gene_haplotype.py` | yes (1000G) | no | **linkage against the gene** |

Scripts 03, 04 and 06 run offline from what is committed here. Script 07 streams
a remote 1000 Genomes index and needs network but no credentials.

### Why there are two linkage scripts

Script 06 tests against G1 and finds every candidate independent. That result is
real but incomplete: G2 is a six-base deletion, Ensembl excludes indels from its
linkage service, and the deletion is absent from the 1000 Genomes phase 3 call
set entirely. Script 07 therefore tests against all common gene-body variation,
which contains both G1 and G2, and reaches a materially less comfortable
conclusion. **Script 07 supersedes script 06.** Script 06 is kept because the
difference between them is itself worth reporting.

### The API key

Script 02 requires an AlphaGenome API key, obtained under the developers'
non-commercial research terms. **The key is never stored in this repository.**
It is read from `~/.alphagenome/api_key` (mode 600) or `ALPHAGENOME_API_KEY` by
`analysis/ag_auth.py`, which refuses a key whose permissions are looser than 600.

### Not committed

`data/gnomad_region_raw.json` (69 MB) is the raw API response and is regenerable
with script 01. The derived, frozen variant set is committed.

---

## Licence: read this before reusing anything

**The code and the data are under different licences, and one forbids commercial
use.** The MIT licence in `LICENSE` covers the analysis code only. The
AlphaGenome-derived predictions in `results/` are available for **non-commercial
use only**; the binding restrictions and the list of affected files are in
**`LEGALLY_BINDING_TERMS_OF_USE.txt`**.

ENCODE peak calls and 1000 Genomes genotypes remain under their own terms.
