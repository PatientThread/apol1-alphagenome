# Reference choice affects positional enrichment of AlphaGenome variant rankings at APOL1-MYH9

**Christopher Lawrence**
Consultant Nephrologist, London, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net


## Abstract

**Objective.** We examined how the choice of measured chromatin-accessibility
reference set affects interpretation of AlphaGenome variant rankings at the
APOL1-MYH9 locus. This exploratory evaluation arose from a regulatory-variant
prioritisation project and included a beta-globin comparison locus.

**Results.** Among 1,930 common alternate alleles at 1,921 positions, the top
200 variants ranked by the absolute podocyte-labelled differential accessibility
score were enriched in whole-kidney peaks (odds ratio 5.98, 95% CI 4.18-8.52)
and in a podocyte-inclusive kidney union (8.76, 6.30-12.20). The podocyte
estimate exceeded the strongest non-renal estimate numerically under both
references, 5.98 against 5.46 for lung and 8.76 against 5.61, with greater
separation under the union. We detected no change in the
podocyte-minus-comparator enrichment contrast for any of six comparators after
partitioning that union, under a conditional peak bootstrap. At the comparison
locus the corresponding estimates were 1.44 (0.16-6.56) and 7.54 (3.66-15.31).
These observations illustrate reference-dependent interpretation of positional
enrichment. They do not establish tissue-specific variant-effect accuracy,
equivalent performance between loci, or generalisation beyond
training-overlapping data.


## Introduction

Two coding haplotypes at APOL1, G1 and G2, account for much of the excess risk
of non-diabetic kidney disease in people of recent African genetic ancestry [1],
an association first mapped to the adjacent gene MYH9 [2]. Penetrance is
incomplete, and transgenic expression of the risk variants in podocytes causes
dose-dependent kidney disease in mice [3], so regulatory modifiers of APOL1
expression are a plausible contributor to who develops disease [4]. The locus is
overwhelmingly non-coding and that variation is largely uncharacterised.

Sequence-to-function models predict regulatory activity from DNA and are
proposed for prioritising such variants [5, 6]. Earlier models of this family
perform less well in regions that distinguish cell types, in an evaluation that
stratified by cell-type specificity [7].

We set out to prioritise regulatory candidates at this locus. That attempt did
not yield a defensible candidate list, but in the course of it we found that the
apparent answer to a simple question, whether a tissue-labelled ranking
concentrates variants in measured accessible chromatin, changed with the
reference set used to evaluate it. We report that sensitivity, because the
reference sets involved are all defensible choices and the differences between
them are easy to overlook.

## Main text

### Methods

**Variants.** The interval spanning APOL3 to MYH9, chr22:36,140,330-36,388,018
(GRCh38, 247,688 bp), fits within one model sequence context. Variants were
retrieved from the gnomAD v4 GraphQL API [9] (dataset `gnomad_r4`, region query on
GRCh38), taking allele counts from the `afr` population block of the genome
callset and falling back to the exome block where the genome block was absent,
retrieved 13 September 2026; the raw API response is deposited.
Retaining single-nucleotide substitutions at African-ancestry allele frequency
of at least 1% gave **1,930 distinct alternate alleles at 1,921 positions**. The
unit of analysis is a unique build, chromosome, position, REF and ALT
combination; nine positions are multi-allelic and contribute two alleles each. A
deduplication ledger is deposited.

**Predictions.** Predictions came from the AlphaGenome public API [6] using the
`ALL_FOLDS` model version, the default when none is specified, with the
recommended DNase variant scorer. Each variant was scored in its own
1,048,576 bp context placed so that the variant is central
(`start = max(0, pos - 524288)`), using `alphagenome` client 0.9.0. Predictions
were obtained on 13 September 2026 at the test locus and 15 September 2026 for
the comparator rankings and the comparison locus. Rankings use the **absolute value of the raw
differential reference-versus-alternate score**, sorted descending, with no
quantile transformation; ties are broken arbitrarily by the sort. Seven tissue
outputs were ranked: podocyte (glomerular visceral epithelial cell), whole
kidney, hepatocyte, liver, lung, brain and stomach. These seven are a different
list from the six tissues used below to define shared peaks. Because the scorer
is differential, the endpoint tests whether predictions concentrate at
accessible positions, not which allele changes accessibility.

**This evaluation is not independent of model training.** AlphaGenome is trained
on ENCODE DNase-seq coverage, and every kidney biosample used here corresponds
to a named training experiment (ENCSR206OJJ, ENCSR785BDQ, ENCSR000EPW,
ENCSR000EOL, ENCSR000EOK, ENCSR175IWT). The served `ALL_FOLDS` version is
described as a distilled model whose teachers were trained across all folds, so
no genomic region is held out from it; the API also exposes `FOLD_0` to
`FOLD_3`. We did not establish which released checkpoint, if any, holds this
locus out, so we make no claim about a held-out evaluation in either direction.
What follows are consistency checks between predictions and the data underlying
them.

**Measured reference sets.** ENCODE DNase-seq peak files (GRCh38, released, peak
output) [10] were merged over each locus. Coordinates are GRCh38 BED half-open
intervals; touching intervals merge (next start not greater than current end);
intervals are not clipped at the locus boundary. Five reference sets were
evaluated (Table 1): whole kidney alone (`term_name=kidney`); a
**podocyte-inclusive kidney union** adding
`glomerular visceral epithelial cell`; that union partitioned into
**kidney-restricted** and **shared**; and both whole-kidney and
podocyte-inclusive definitions at the comparison locus. A merged kidney peak is
classified **shared** if a peak from any of six comparison tissues (liver, lung,
heart left ventricle, stomach, spleen, thyroid gland; four files per tissue)
**covers its midpoint**, and **kidney-restricted** otherwise. This is a
midpoint rule, so part of a kidney-restricted peak may still overlap a
comparison peak. Partitions are mutually exclusive and sum to the union in
peaks, base pairs and alleles, asserted in code. Every row of Table 1 uses the
same allele set and the same podocyte ranking, held fixed across reference
definitions; comparator rankings necessarily use different tissue outputs. A
peak manifest listing every file, term and contribution is deposited.

**Statistics.** Enrichment of the top-ranked alleles was tested against the
remaining common alleles at the same locus by Fisher exact test, with
conditional maximum-likelihood odds ratios and exact intervals; all four
contingency cells and the denominator are deposited for every endpoint. Depth
200 is primary; 25, 50 and 100 are nested sensitivity analyses. These intervals
assume independence between alleles, which clustering within peaks violates, so
they are descriptive summaries rather than tests of ranking differences.

**Resampling.** Comparisons between rankings resample **kidney-restricted merged
peaks** with replacement (500 draws, seed 20260920, `numpy` PCG64). All alleles
outside a kidney-restricted peak are held fixed, which includes alleles in
shared kidney peaks. Alternate alleles at one position share a position and
therefore always travel together within a resampling unit. A peak drawn k times
contributes its alleles k times to both the top set and the denominator, and
top-set membership is recomputed within each resample rather than held fixed.
Bootstrap log odds ratios carry a Haldane 0.5 correction so no replicate is
infinite. Intervals are percentile, unadjusted across six comparators, and
conditional on the fixed background. Interval endpoints were rechecked at 500,
2,000 and 5,000 draws across three seeds. Across all 54 resulting settings no
interval excluded zero and the largest shift in any endpoint was 0.325, in the
widest interval (liver); the full stability table is deposited.

**Comparison locus.** The beta-globin cluster, chr11:5,150,000-5,397,688,
matched on length and paralogue-cluster structure and **approximately
count-matched** by subsampling stratified on African-ancestry allele frequency
in eight bins. This gives 1,929 alleles against 1,930: in the 0.10 to 0.20
stratum the test locus required 463 and the control contained 462. No exact
matching was attempted at any stage.

### Results

**Whole-kidney reference.** Against whole-kidney peaks alone (67 peaks,
30,358 bp, 12.3% of the locus) the podocyte ranking gave an odds ratio of 5.98
(95% CI 4.18-8.52) at depth 200. The strongest non-renal ranking at that depth
was lung, 5.46 (3.81-7.79). The podocyte estimate is therefore numerically above
the strongest non-renal estimate, by a small margin. In a nested sensitivity
analysis at depth 100 that ordering reversed: podocyte 5.25 and hepatocyte 6.72.

**Podocyte-inclusive union.** Against the union (110 peaks, 44,576 bp, 18.0% of
the locus) the podocyte ranking gave 8.76 (6.30-12.20) against 5.61 (4.02-7.81)
for lung (Figure 1). The separation is larger than under the whole-kidney
reference. That is a difference in the size of a numerical gap, not a
demonstration that discrimination is absent under one reference and established
under the other.

**Partitions.** The union partitions into 39 kidney-restricted peaks
(14,686 bp) and 71 shared (29,890 bp), containing 84 and 193 of the 277 alleles
in the union. The podocyte ranking led in both, 4.30 (2.53-7.15) and 7.94
(5.55-11.33). The change in the podocyte-minus-comparator contrast on
partitioning ranged from -0.07 to +0.42 in natural-log odds ratio and no
interval excluded zero for any of six comparators (Figure 2, all six values
deposited). We detected no change in relative enrichment on partitioning under
this resampling scheme. This is not evidence of equivalence, and it answers a
different question from the whole-kidney to union comparison above.

**Comparison locus.** Against whole-kidney peaks the beta-globin cluster gave
1.44 (0.16-6.56, one-sided p = 0.43; 10 peaks, 2,598 bp). Against the
podocyte-inclusive union it gave 7.54 (3.66-15.31; 18 peaks, 5,851 bp). Positional
enrichment is therefore **also observed at the comparison locus** under the
union definition. We do not read this as equivalent performance between loci.
The loci differ substantially in the architecture being evaluated, with union
coverage of 18.0% at APOL1-MYH9 against 2.4% at beta-globin, and overlapping
intervals are not an equivalence test. Equally, the wide whole-kidney interval
at the comparison locus, 0.16 to 6.56, is too imprecise to establish an absence
of enrichment there.

Table 1 sets out each reference set with its denominator and named comparator.

### Discussion

Three interpretations we held during this work did not survive the reconciled
analysis, and they are related but not identical failures.

The first two are endpoint definitions. Whole-kidney and podocyte-inclusive
references describe different positional endpoints: a cell-type ranking
evaluated against a tissue-level measurement is not the same comparison as that
ranking evaluated against a reference containing the cell type. Adding podocyte
files changes reference coverage and peak composition as well as cell-type
resolution, so the union is better described as a broader, explicitly defined
reference than as the uniquely matched one.

The third was an internal consistency error: an earlier analysis compared
partitions built on one reference against a pooled set built on another, so the
partitions did not sum to their stated parent. We have not reproduced it here
and readers should not need it. The forward-looking point is that asserting in
code that partitions sum to their parent is a cheap safeguard, and we now do
so.

A fourth interpretation concerned the comparison locus. A beta-globin interval
with measured renal accessibility is not a zero-accessibility control, and an
accessibility model may reasonably enrich variants at accessible positions
there. Reading a wide, imprecise null interval under one reference as evidence
of locus specificity was an inferential error; that reference was not directly
comparable for the question.

What survives is narrow: a podocyte-labelled ranking concentrates alleles in
measured podocyte-containing accessible regions at both loci examined, on data
the model was trained on. Three practical points follow. Match the measured
reference to the predicted output and report both explicitly. Assert that any
partition sums to its parent before interpreting differences between them. Use
the same reference definition at test and comparison loci, or the comparison is
not informative for the question asked.

## Limitations

This is a single-locus exploratory analysis with no pre-registered plan, and the
comparisons reported were arrived at iteratively.

The evaluation is not independent of model training, which is the binding
limitation: the tracks are named training experiments, and the served model is
distilled from teachers trained across all folds. We did not identify a released
checkpoint holding this locus out, so these results speak to internal
consistency and not to generalisation.

The endpoint is positional. A differential reference-versus-alternate score
ranked against peak membership does not test which allele changes accessibility,
nor the direction or size of any effect, and we did not compare it against a
reference-accessibility baseline, so we cannot say what the differential scorer
adds over ranking by predicted accessibility alone.

The bootstrap resamples kidney-restricted peaks and holds everything else fixed,
including alleles in shared peaks, so its intervals are conditional on that
background rather than a general estimate of sampling uncertainty. They are
unadjusted for six comparisons. Peak resampling reduces but does not remove
dependence between neighbouring peaks. A non-significant interaction is not
evidence of no interaction, and we specify no margin below which a change would
be negligible.

"Kidney-restricted" is defined by a midpoint rule against six comparison
tissues and would change with a different panel or rule. The comparison locus is
approximately rather than exactly count-matched.

The prioritisation exercise that prompted this work did not yield a candidate
list we were prepared to report under the specified procedure of linkage and
prior-literature screening; those intermediate analyses are deposited rather
than presented, and nothing here should be read as establishing or excluding an
APOL1 regulatory modifier.

## Conclusion

Reference choice affected the observed positional enrichment and its
interpretation at the two evaluated loci. Consistent reference definitions and
reconciled partitions are necessary for interpretable comparisons. These
training-overlapping evaluations do not establish allele-specific regulatory
accuracy or an APOL1 modifier.

## Abbreviations

ATAC: assay for transposase-accessible chromatin; bp: base pairs; CI: confidence
interval; CKD: chronic kidney disease; DNase: deoxyribonuclease; ENCODE:
Encyclopedia of DNA Elements; eQTL: expression quantitative trait locus; GRCh38:
Genome Reference Consortium Human Build 38; OR: odds ratio; REF/ALT:
reference/alternate allele.

## Declarations

**Ethics approval and consent to participate.** Not applicable. No new
participant data were collected; the study analysed existing publicly available
ENCODE and gnomAD data together with model predictions.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All analysis code, the frozen allele
ledger, per-variant predictions, contingency counts for every endpoint, the
ENCODE peak manifest, partition labels, the ranking provenance record, the six
change-in-contrast estimates with intervals, and the bootstrap stability table
are available in the apol1-alphagenome repository [11]. Variant frequencies are
from gnomAD v4 [8, 9]; peak calls from the ENCODE portal [10]. Model predictions were
obtained under AlphaGenome's non-commercial research terms and are redistributed
subject to them; model weights were not accessed.

**Competing interests.** The author declares no competing interests.

**Funding.** No funding was received.

**Authors' contributions.** CL is the sole author and conceived the study, wrote
the code, performed the analyses and wrote the manuscript.

**Acknowledgements.** Not applicable.

## References

1. Genovese G, Friedman DJ, Ross MD, Lecordier L, Uzureau P, Freedman BI, et al. Association of trypanolytic ApoL1 variants with kidney disease in African Americans. Science. 2010;329:841-5.
2. Kopp JB, Smith MW, Nelson GW, Johnson RC, Freedman BI, Bowden DW, et al. MYH9 is a major-effect risk gene for focal segmental glomerulosclerosis. Nat Genet. 2008;40:1175-84.
3. Beckerman P, Bi-Karchin J, Park AS, Qiu C, Dummer PD, Soomro I, et al. Transgenic expression of human APOL1 risk variants in podocytes induces kidney disease in mice. Nat Med. 2017;23:429-38.
4. Ojo AO, Adu D, Bramham K, Freedman BI, Gbadegesin RA, Ilori TO, et al. APOL1 kidney disease: conclusions from a Kidney Disease: Improving Global Outcomes (KDIGO) Controversies Conference. Kidney Int. 2025;108:763-79.
5. Linder J, Srivastava D, Yuan H, Agarwal V, Kelley DR. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. Nat Genet. 2025;57:949-61.
6. Avsec Ž, Latysheva N, Cheng J, Novati G, Taylor KR, Ward T, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature. 2026;649:1206-18.
7. Kathail P, Shuai RW, Chung R, Ye CJ, Loeb GB, Ioannidis NM. Current genomic deep learning models display decreased performance in cell type-specific accessible regions. Genome Biol. 2024;25:202.
8. Chen S, Francioli LC, Goodrich JK, Collins RL, Kanai M, Wang Q, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625:92-100.
9. Genome Aggregation Database. gnomAD v4, GraphQL API dataset identifier `gnomad_r4`, GRCh38 region query, `afr` genetic-ancestry group, genome callset with exome fallback. https://gnomad.broadinstitute.org/api (accessed 13 September 2026). The raw API response is deposited with the analysis archive.
10. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature. 2020;583:699-710. Individual experiment and file accessions are listed in the deposited peak manifest.
11. Lawrence C. apol1-alphagenome: analysis code and frozen results, release v1.0-bmcrn. https://github.com/PatientThread/apol1-alphagenome/releases/tag/v1.0-bmcrn (accessed 21 September 2026).

## Figure legends

**Figure 1.** Positional enrichment of the top 200 alleles by the tissue output
used to rank, across the six locus and reference combinations of Table 1.
Panels: (A) APOL1-MYH9, whole kidney; (B) APOL1-MYH9, podocyte-inclusive union;
(C) union, kidney-restricted; (D) union, shared; (E) beta-globin, whole kidney;
(F) beta-globin, podocyte-inclusive union. Conditional maximum-likelihood odds
ratios with exact 95% intervals against the remaining common alleles at the same
locus, all at depth 200; eligible n = 1,930 at APOL1-MYH9 and 1,929 at
beta-globin. Common log scale with OR = 1 marked; tissue order identical in
every panel; podocyte and whole kidney labelled separately and drawn in black,
non-renal rankings in grey. Only the podocyte ranking was evaluated at the
comparison locus, so panels E and F show that estimate alone. Intervals assume
independence between alleles.

**Figure 2.** Change in tissue contrast on partitioning, for each of six
comparators. The quantity plotted is the podocyte-minus-comparator contrast in
kidney-restricted regions minus the same contrast in the podocyte-inclusive
union, in natural-log odds ratio, so a positive value means the podocyte ranking
gains relative advantage in kidney-restricted regions. Estimated on identical
resampled draws of kidney-restricted peaks with all other alleles held fixed;
500 replicates, percentile intervals, unadjusted for six comparisons. Numerical
interval limits are printed beside each point. No interval excludes zero.

## Tables

**Table 1.** Positional enrichment at depth 200 for each locus and measured
reference set. The podocyte ranking is held fixed across reference definitions;
the comparator column names the strongest non-renal ranking evaluated at the
same locus, reference and depth. Odds ratios are conditional maximum likelihood
with exact 95% intervals. Full contingency cells for every endpoint and depth
are deposited.

| Locus | Reference set | Peaks | bp (% locus) | n | Podocyte OR (95% CI) | Strongest non-renal OR (95% CI) |
|---|---|---|---|---|---|---|
| APOL1-MYH9 | Whole kidney | 67 | 30,358 (12.3) | 1,930 | 5.98 (4.18-8.52) | lung 5.46 (3.81-7.79) |
| APOL1-MYH9 | Podocyte-inclusive union | 110 | 44,576 (18.0) | 1,930 | 8.76 (6.30-12.20) | lung 5.61 (4.02-7.81) |
| APOL1-MYH9 | Union, kidney-restricted | 39 | 14,686 (5.9) | 1,930 | 4.30 (2.53-7.15) | lung 2.49 (1.36-4.36) |
| APOL1-MYH9 | Union, shared | 71 | 29,890 (12.1) | 1,930 | 7.94 (5.55-11.33) | lung 6.02 (4.18-8.62) |
| beta-globin | Whole kidney | 10 | 2,598 (1.0) | 1,929 | 1.44 (0.16-6.56) | not evaluated |
| beta-globin | Podocyte-inclusive union | 18 | 5,851 (2.4) | 1,929 | 7.54 (3.66-15.31) | not evaluated |

Percentages are of the 247,688 bp evaluated at each locus. Interpretations we
previously drew from the whole-kidney rows have been withdrawn; the recomputed
endpoint-specific estimates in those rows remain valid for the endpoint they
describe.
