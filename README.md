# Regulatory variant prioritisation at the APOL1-MYH9 locus

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

**This is the APOL1-MYH9 locus, not APOL1.** All seven independent candidates are
nearer another gene: five inside MYH9, two in or beside APOL4. The accessibility
channel used for the ranking is position-local and gene-agnostic, so it never
targeted APOL1.

1. **The obvious readout is confounded and is not used.** Ranking by the
   gene-masked expression scorer rediscovers proximity to the gene. The
   pre-registered negative control, the coding risk variants, stayed at the 97th
   and 96th percentile even after residualising on distance.
2. **The accessibility channel is not confounded** and its top-ranked variants
   are enriched in independently measured kidney open chromatin, odds ratios 8.1
   and 6.1 against a locus-matched background.
3. **Linkage disequilibrium removes over half the candidates.** Tested against
   all common APOL1 gene-body variation, which contains both G1 and G2: 7 of 15
   independent, 7 partially correlated, 1 a perfect proxy. Only 4 of the 7 sit in
   measured open chromatin in any kidney cell type.
4. **There is no cell-type specificity.** The ranking predicts podocyte
   accessibility; podocyte leads at one of four thresholds against other kidney
   cell types.
5. **There is no tissue specificity.** Ranking by hepatocyte finds kidney
   chromatin as well as ranking by kidney does. Kidney leads at one of four
   thresholds against eight non-renal tissues.
6. **The pipeline does not manufacture enrichment.** At a matched beta-globin
   control locus the same ranking finds nothing, odds ratio 1.44 at top 200,
   p = 0.44. Power-checked: an APOL1-sized effect would have been detected.

**The primary finding is methodological.** The model detects regulatory positions
competently and detects them where they exist, but does not use the tissue label
to do it. The candidate list is the secondary output.

## What this repository does *not* claim

No candidate has been functionally validated, and all seven are nearer a gene
other than APOL1. The enrichment shows the ranking selects positions that are
open; findings 4 and 5 show it does not select positions that are open
*specifically in kidney*. No variant is shown to change expression and no disease
association is tested. See `FINDINGS.md` for the full list of claims this work
does not support.

---

## Reproducing the analysis

Scripts are numbered in dependency order.

| Script | Needs network | Needs API key | Produces |
|---|---|---|---|
| `01_fetch_locus_variants.py` | yes (gnomAD) | no | the frozen locus variant set |
| `02_score_regulatory_effects.py` | yes | **yes** | per-variant predictions, kidney tracks |
| `03_analyse.py` | no | no | distributions, first ranking |
| `04_distance_controlled.py` | no | no | **the proximity correction** |
| `05_validate_against_encode.py` | yes (ENCODE) | no | enrichment vs whole-kidney peaks |
| `06_ld_with_risk_haplotypes.py` | yes (Ensembl) | no | linkage against G1 only, superseded |
| `07_ld_against_gene_haplotype.py` | yes (1000G) | no | **linkage against the whole gene** |
| `08_annotate_candidates.py` | yes (ENCODE, Ensembl) | no | regulatory elements, nearest gene |
| `09_cell_type_resolved_validation.py` | yes (ENCODE) | no | **enrichment by kidney cell type** |
| `10_wrong_tissue_control.py` | yes | **yes** | scores with eight non-renal tissues |
| `11_negative_control_locus.py` | yes | **yes** | **the negative control locus** |
| `12_analyse_wrong_tissue.py` | no | no | **the wrong-tissue control readout** |

Scripts 03, 04, 06 and 12 run offline from what is committed here. Scripts 07 to
09 and 11 need network but no credentials. Only 02, 10 and 11 need an API key.

**Read `FINDINGS.md` first.** It states what the paper can and cannot claim, with
every number sourced from `results/`.

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
