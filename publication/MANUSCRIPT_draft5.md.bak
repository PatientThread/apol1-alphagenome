# Evaluating tissue-labelled AlphaGenome rankings at APOL1-MYH9: how the choice of measured reference changes the conclusion

**Christopher Lawrence**
Consultant Nephrologist, London, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net

**Keywords:** APOL1; MYH9; chromatin accessibility; benchmarking; machine
learning

---

## Abstract

**Objective.** Sequence-to-function models are increasingly used to prioritise
non-coding variants at kidney disease loci, and their tissue-labelled outputs are
read as tissue-specific. We evaluated AlphaGenome rankings at the APOL1-MYH9
locus and report how three successive conclusions were determined not by the
model but by which measured accessibility dataset we compared them against.

**Results.** Ranking 1,930 common variants by predicted podocyte accessibility
and evaluating against measured whole-kidney peaks alone, non-renal rankings
matched the podocyte ranking, suggesting the tissue label carried little
information. Against a kidney-plus-podocyte union the podocyte ranking led
clearly (odds ratio 8.76, 95% CI 6.30-12.20). Partitioning that union into
kidney-restricted and shared regions changed the podocyte-minus-comparator
contrast for none of six comparators. A beta-globin control locus appeared null
against whole-kidney peaks (odds ratio 1.44) but showed comparable enrichment to
the test locus against the matched union (7.54, 3.66-15.31). Each apparent
finding was an artefact of comparing a cell-type prediction against a
non-matching measurement. The evaluated region and tracks overlap model training
data, so none of these comparisons tests generalisation.

---

## Introduction

Two coding haplotypes at APOL1, G1 and G2, account for much of the excess risk of
non-diabetic kidney disease in people of recent African genetic ancestry [1]. The
association was first mapped to the adjacent gene MYH9 [2] and resolved to APOL1
two years later [1]. Penetrance is incomplete, and transgenic expression of the
risk variants in podocytes causes kidney disease in mice in a dose-dependent
manner [3], so regulatory modifiers of APOL1 expression are a plausible
explanation for who develops disease [4]. The locus is overwhelmingly non-coding
and that variation is largely uncharacterised.

Sequence-to-function models predict regulatory activity from DNA and are proposed
for prioritising such variants [5, 6]. Earlier models of this family perform less
well in regions that distinguish cell types, in an evaluation that already
stratified by cell-type specificity [7].

We set out to prioritise regulatory candidates at this locus. That attempt did
not succeed, and in the course of it we reached three successive conclusions
about tissue specificity, each of which we subsequently found to be an artefact
of the same design choice: the measured dataset used for evaluation did not match
the prediction being evaluated. We report this because each intermediate result
was internally coherent, none was implausible, and the error is easy to repeat.

## Main text

### Methods

The interval spanning APOL3 to MYH9, chr22:36,140,330-36,388,018 (GRCh38,
247,688 bp), fits within one model context window. Variants catalogued by gnomAD
v4 [8] were retrieved with alternate allele frequencies for the African
genetic-ancestry group. Retaining single-nucleotide substitutions at frequency at
least 1% gave **1,930 distinct variants at 1,921 positions**; nine positions are
multi-allelic, carrying two alternate alleles under one identifier, and are
counted separately. A deduplication ledger is deposited.

Predictions came from the AlphaGenome public API [6] using the **ALL_FOLDS**
model version, the default when none is specified, with position-local DNase
scorers. These are differential reference-versus-alternate scorers; the endpoint
below therefore tests whether predictions concentrate at accessible positions, not
which allele changes accessibility.

**This evaluation is not independent of model training.** AlphaGenome is trained
on ENCODE DNase-seq coverage and every kidney biosample used here corresponds to a
named training experiment (ENCSR206OJJ, ENCSR785BDQ, ENCSR000EPW, ENCSR000EOL,
ENCSR000EOK, ENCSR175IWT). The ALL_FOLDS model is trained on all eight genomic
sections with no held-out region; applying the published fold definitions this
locus lies in data fold 6, and no released checkpoint designates it as its test
region. These are consistency checks between predictions and the data underlying
them.

ENCODE DNase-seq peak files (GRCh38, released) [9] were merged over each locus.
Three measured reference sets were evaluated (Table 1): whole kidney alone
(`term_name=kidney`); whole kidney plus podocyte (adding
`glomerular visceral epithelial cell`); and that union partitioned into
**kidney-restricted** and **shared**, a kidney peak counting as shared if its
midpoint lies inside a peak from any of six comparison tissues (liver, lung,
heart left ventricle, stomach, spleen, thyroid gland), four files per tissue.
Partitions are mutually exclusive and sum to the union in peaks, base pairs and
variants, asserted in code. "Kidney-restricted" means no peak detected in the
selected panel, not closure. Every row of Table 1 uses the same variant universe,
the same predicted podocyte score and the same interval-merging convention, so
the rows differ only in the measured reference set.

Enrichment of the top-ranked variants was tested by Fisher exact test against the
other common variants at the same locus, with conditional maximum-likelihood odds
ratios and exact intervals. All four contingency cells and the eligible
denominator are deposited for every endpoint. Depth 200 is the primary analysis;
25, 50 and 100 are nested sensitivity analyses.

Because variants within one peak share the element and have correlated scores,
comparisons between rankings resample **peaks** with replacement, variants outside
peaks held fixed, 500 replicates, percentile intervals, unadjusted across six
comparisons. Uncertainty is therefore conditional on that fixed background. The
change in contrast between pooled and kidney-restricted endpoints was estimated on
the same draws.

The control locus was the beta-globin cluster, chr11:5,150,000-5,397,688, matched
on length, paralogue-cluster structure and, by frequency-stratified subsampling,
variant count.

### Results

**Against whole-kidney peaks alone** (67 peaks, 30,358 bp, 12.3% of the
locus), ranking by predicted **podocyte** accessibility gave an odds ratio of 5.98
(95% CI 4.18-8.52) at the primary depth, and non-renal rankings were of the same
order: lung 5.46 at that depth, and at depth 100 hepatocyte exceeded podocyte
(6.72 against 5.25). Read alone, this suggests the tissue label carries little
information. The endpoint evaluates a cell-type prediction against a
tissue-level measurement.

**Against the kidney-plus-podocyte union** (110 peaks, 44,576 bp, 18.0% of the
locus) the podocyte ranking led at the primary depth: odds ratio 8.76
(95% CI 6.30-12.20; 100 of the top 200 in peaks against 177 of 1,730). The
strongest non-renal ranking, lung, gave 5.61, and the whole-kidney ranking 2.00
(Figure 1).

**Partitioning that union** gave 39 kidney-restricted peaks (14,686 bp) and 71
shared (29,890 bp), containing 84 and 193 variants of the 277 in the union. The
podocyte ranking led in both partitions, 4.30 and 7.94. We had expected
partitioning to reveal discrimination that pooling hid. It did not: the change in
the podocyte-minus-comparator contrast ranged from -0.07 to +0.42 in log odds
ratio and no interval excluded zero for any of six comparators (Figure 2). Under
this resampling scheme we detected no change in relative enrichment on
partitioning.

**The control locus reverses with the same substitution.** Against whole-kidney
peaks the beta-globin cluster appeared null (10 peaks, 2,598 bp; odds ratio 1.44,
95% CI 0.16-6.56, one-sided p = 0.43), which we had read as evidence that the pipeline does not
produce enrichment indiscriminately. Against the matched kidney-plus-podocyte
union (18 peaks, 5,851 bp, 2.4% of that locus) it showed enrichment comparable to
the test locus at every depth: 18.36 (5.62-52.06) at depth 25 falling to 7.54
(3.66-15.31) at depth 200, against 8.76 at the test locus. The intervals overlap
substantially. A confidence interval from one locus against a point estimate from
another is not a test of their difference, and we make no such comparison; the
descriptive point is that the control is not null.

Table 1 sets out each measured reference set, its denominator and the conclusion
it supported.

### Discussion

Three conclusions, each internally coherent, were determined by which measured
dataset the predictions were compared against rather than by the predictions. The
common defect is evaluating a podocyte-level prediction against a measurement made
at a different level of cell-type resolution, or against a differently constituted
reference set.

Two of the three were detectable by arithmetic once suspected. Partitions that do
not sum to their parent set in peaks, base pairs and variants indicate that the
comparison is not of the same data, and asserting that sum in code would have
caught it. The first and third were not detectable that way: both parent sets were
internally consistent, and only matching the reference to the prediction exposes
them.

What survives is limited. The podocyte-labelled ranking is enriched in measured
podocyte-containing accessible regions, at a kidney disease locus and at an
erythroid control locus to a similar degree. That is the behaviour expected of a
model trained to predict accessibility, evaluated on data it was trained on, at
regions within its training folds. It says nothing about generalisation,
allele-specific accuracy, or APOL1.

For investigators evaluating these models at disease loci, three practical points
follow. Match the measured reference to the predicted track at the same level of
cell-type resolution, and report both explicitly. Assert in code that any
partition sums to its parent before interpreting differences between partitions.
Use the same measured reference at test and control loci, or the control tests
nothing.

## Limitations

This is a single-locus, exploratory analysis with no pre-registered plan, and the
comparisons reported were arrived at iteratively.

The evaluation is not independent of model training, which is the binding
limitation: the tracks are training experiments and the locus lies in a training
fold of the served checkpoint. No released checkpoint designates this region as
its test partition, so a genuine generalisation test was not available.

The endpoint is positional. A differential reference-versus-alternate score
ranked against peak membership does not test which allele changes accessibility,
nor the direction or size of any effect, and we did not compare it against a
simple reference-accessibility ranking to establish what the differential scorer
adds.

The bootstrap holds variants outside peaks fixed, so its intervals are conditional
on that background and are not a general estimate of sampling uncertainty. They
are unadjusted for six comparisons. A non-significant interaction is not evidence
of no interaction; we specify no margin below which a change would be negligible.

"Kidney-restricted" is defined relative to six comparison tissues and would change
with a different panel. Candidate prioritisation is not reported here: no
candidate survived linkage and prior-literature scrutiny as a novel independent
variant, and those analyses are deposited rather than presented.

## Declarations

**Availability of data and materials.** All analysis code, the frozen variant
ledger, per-variant predictions, contingency counts for every endpoint, ENCODE
file accessions, partition labels and the interval-to-fold mapping are available
at https://github.com/PatientThread/apol1-alphagenome, with a tagged release for
review. Variant frequencies are from gnomAD v4; peak calls from the ENCODE
portal. Model predictions were obtained under AlphaGenome's non-commercial
research terms and are redistributed subject to them; model weights were not
accessed.

**Ethics approval and consent to participate.** Not applicable. No new
participant data were collected. The study analysed existing publicly available
ENCODE and gnomAD data together with model predictions.

**Consent for publication.** Not applicable.

**Competing interests.** The author declares no competing interests.

**Funding.** No funding was received.

**Authors' contributions.** CL is the sole author and conceived the study, wrote
the code, performed the analyses and wrote the manuscript.

## References

1. Genovese G, Friedman DJ, Ross MD, Lecordier L, Uzureau P, Freedman BI, et al. Association of trypanolytic ApoL1 variants with kidney disease in African Americans. Science. 2010;329:841-5.
2. Kopp JB, Smith MW, Nelson GW, Johnson RC, Freedman BI, Bowden DW, et al. MYH9 is a major-effect risk gene for focal segmental glomerulosclerosis. Nat Genet. 2008;40:1175-84.
3. Beckerman P, Bi-Karchin J, Park AS, Qiu C, Dummer PD, Soomro I, et al. Transgenic expression of human APOL1 risk variants in podocytes induces kidney disease in mice. Nat Med. 2017;23:429-38.
4. Ojo AO, Adu D, Bramham K, Freedman BI, Gbadegesin RA, Ilori TO, et al. APOL1 kidney disease: conclusions from a Kidney Disease: Improving Global Outcomes (KDIGO) Controversies Conference. Kidney Int. 2025;108:763-79.
5. Linder J, Srivastava D, Yuan H, Agarwal V, Kelley DR. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. Nat Genet. 2025;57:949-61.
6. Avsec Ž, Latysheva N, Cheng J, Novati G, Taylor KR, Ward T, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature. 2026;649:1206-18.
7. Kathail P, Shuai RW, Chung R, Ye CJ, Loeb GB, Ioannidis NM. Current genomic deep learning models display decreased performance in cell type-specific accessible regions. Genome Biol. 2024;25:202.
8. Chen S, Francioli LC, Goodrich JK, Collins RL, Kanai M, Wang Q, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625:92-100.
9. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature. 2020;583:699-710.

## Figure legends

**Figure 1.** Enrichment of the top 200 variants in measured accessible regions,
by the tissue output used to rank, across the five measured reference sets of
Table 1. Conditional odds ratios with exact 95% intervals against the other
common variants at the same locus; eligible n = 1,930 at the test locus and
1,929 at the control. Panels, left to right: whole kidney only; kidney plus
podocyte; that union's kidney-restricted and shared partitions; and the control
locus against the same union definition. Kidney rankings are drawn in black
circles, non-renal in grey squares. The leftmost panel is the comparison that
suggested the tissue label carried little information and the rightmost is the
control that was read as null; both change with the reference set and neither
conclusion survives. Intervals assume variant-level independence.

**Figure 2.** Change in tissue contrast on partitioning. For each comparator, the
podocyte-minus-comparator contrast in kidney-restricted regions minus the same
contrast in the pooled union, in natural-log odds ratio, estimated on identical
resampled draws of peaks with variants outside peaks held fixed. Points are
observed values; intervals are percentile bootstrap, 500 replicates, unadjusted
for six comparisons. No interval excludes zero.

## Tables

**Table 1.** Each measured reference set evaluated, and the conclusion it
supported. All rankings use the same predicted podocyte accessibility score and
the same variant sets. Odds ratios are conditional maximum likelihood at ranking
depth 200 with exact 95% intervals. Full contingency counts for every endpoint
and depth are deposited.

| Measured reference set | Locus | Peaks | bp (% of locus) | Eligible n | Podocyte OR (95% CI) | Best non-renal OR | Conclusion it supported |
|---|---|---|---|---|---|---|---|
| Whole kidney only | APOL1-MYH9 | 67 | 30,358 (12.3) | 1,930 | 5.98 (4.18-8.52) | 6.72 (hepatocyte, depth 100) | label carries little information |
| Kidney + podocyte | APOL1-MYH9 | 110 | 44,576 (18.0) | 1,930 | 8.76 (6.30-12.20) | 5.61 (lung) | podocyte ranking leads clearly |
| Kidney + podocyte, kidney-restricted | APOL1-MYH9 | 39 | 14,686 (5.9) | 1,930 | 4.30 (2.53-7.15) | 2.49 | podocyte ranking leads |
| Kidney + podocyte, shared | APOL1-MYH9 | 71 | 29,890 (12.1) | 1,930 | 7.94 (5.55-11.33) | 6.02 | podocyte ranking leads |
| Whole kidney only | beta-globin | 10 | 2,598 (1.0) | 1,929 | 1.44 (0.16-6.56) | not evaluated | pipeline is locus-specific |
| Kidney + podocyte | beta-globin | 18 | 5,851 (2.4) | 1,929 | 7.54 (3.66-15.31) | not evaluated | pipeline is not locus-specific |

The first and fifth rows are superseded and are shown to document how the
conclusion changed with the reference set.
