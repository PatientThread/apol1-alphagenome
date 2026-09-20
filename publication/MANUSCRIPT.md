<!-- DRAFT 4. Rewritten after the second editorial review, which found that the
     pooled and partitioned kidney sets were different data, so the draft-3
     headline was unsupported. Reconciled onto one frozen universe (script 16),
     the podocyte ranking leads in every view and partitioning changes the tissue
     contrast for none of six comparators. TWO SUCCESSIVE HEADLINES WERE THE SAME
     ARTEFACT: mismatching the measured validation set to the predicted track.
     That is now the paper. Claims are deliberately small. -->

# Interpreting tissue-labelled sequence-model rankings at APOL1-MYH9: a benchmark-design case study

**Christopher Lawrence**
Consultant Nephrologist, London, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net

**Running head:** Benchmark design at APOL1-MYH9

**Keywords:** APOL1; chronic kidney disease; gene expression regulation;
machine learning; MYH9

---

## Abstract

**Background and hypothesis.** Sequence-to-function models are increasingly
applied to non-coding variation at kidney disease loci, and their tissue-labelled
outputs are read as tissue-specific. We examined what such rankings can and
cannot establish at the APOL1-MYH9 locus, and report two analytical errors we
made before arriving at an interpretable result.

**Methods.** We scored 1,921 variants common in the African genetic-ancestry
group across chr22:36,140,330-36,388,018 using AlphaGenome, and tested rankings
for enrichment in measured ENCODE accessible regions against the other common
variants at the same locus. Measured kidney regions were partitioned into those
with no peak detected in six comparison tissues and those shared with at least
one. Contrasts between rankings were estimated pairwise with peaks as the
resampling unit, and the change in contrast between pooled and partitioned
endpoints was estimated on the same resampled draws. Linkage disequilibrium with
G1 and G2 was computed from phased haplotypes in seven African populations,
reporting D' and the maximum r-squared the allele frequencies permit.

**Results.** Using one frozen variant set and one measured kidney union, the
podocyte ranking led in every view: pooled odds ratio 8.94, kidney-restricted
4.28, shared 8.13, against 5.58, 2.48 and 5.98 for the strongest non-renal
ranking. The change in contrast on partitioning excluded zero for none of six
comparators. Two earlier conclusions, that the tissue label was uninformative and
that partitioning reversed this, were both artefacts of validating podocyte
predictions against measured sets that did not match them.

**Conclusions.** In this exploratory single-locus analysis the podocyte-labelled
ranking showed consistent enrichment in measured kidney accessible regions.
Because the evaluated region and tracks overlap model training data, this concerns
benchmark construction and does not establish generalisation, allele-specific
accuracy or APOL1 modifiers.

---

## Key learning points

**What was known**

- Two coding haplotypes at APOL1, G1 and G2, carry most of the excess kidney
  disease risk in people of recent African ancestry, and disease depends on
  expression dose.
- Sequence-to-function models are proposed for interpreting non-coding variation
  at such loci, and earlier models perform less well in regions that distinguish
  cell types.
- Regulatory variation at this locus is under active investigation, including
  variants associated with circulating APOL1 and methylation differences across
  the APOL1-APOL4-MYH9 region.

**This study adds**

- A worked demonstration that mismatching the measured validation set to the
  predicted track can reverse the apparent conclusion, in either direction; we
  made that error twice before detecting it by arithmetic.
- Evaluation of this model against ENCODE data at this locus is not independent:
  the tracks used are named training experiments and the region lies in a
  training fold of the served checkpoint.
- Low r-squared with G1 or G2 does not establish haplotypic independence when
  allele frequencies differ, and most candidates failed once D' was examined.

**Potential impact**

- Investigators should match the measured evaluation set to the predicted track
  and verify that partitions sum to their parent set before interpreting
  differences between them.
- Model version and the fold containing the locus should be reported as routinely
  as software versions, since the default served checkpoint holds nothing out.
- Nothing here identifies a modifier of APOL1 penetrance, a target gene or a
  clinically actionable variant, and none is claimed.

---

## 1. Introduction

Two coding haplotypes at APOL1, G1 and G2, account for much of the excess risk of
non-diabetic kidney disease in people of recent African genetic ancestry [1]. The
association was first mapped to the adjacent gene MYH9 [2] and resolved to APOL1
two years later, the MYH9 signal being attributed to linkage disequilibrium [1].

Penetrance is incomplete. Transgenic expression of the risk variants in podocytes
causes kidney disease in mice in a dose-dependent manner [3], so genetic modifiers
of APOL1 expression are a plausible explanation for who develops disease [4].
Regulatory variation at this locus is under active investigation: cis-regulatory
variants have been associated with circulating APOL1 concentration [5], and APOL1
risk genotypes with methylation differences at enhancers and promoters across the
APOL1-APOL4-MYH9 region [6]. Kidney expression resources have improved, including
maps from 659 microdissected kidney samples [7], though ancestry representation
and variant coverage remain limiting.

Sequence-to-function models predict regulatory activity from DNA and are proposed
for prioritising non-coding variants [8, 9]. Earlier models of this family capture
promoter determinants while largely ignoring distal enhancers [10] and perform
less well in regions that distinguish cell types, an evaluation that already
stratified by cell-type specificity [11]. This study is a worked example of that
principle applied to one model at one locus, not a discovery of it.

We set out to prioritise regulatory candidates. That attempt did not succeed. What
it produced instead is a case study in how a single-locus benchmark can be
misread, and we report our own errors because they were not obvious and each
produced a clean result.

## 2. Materials and methods

### 2.1 Locus, variants and candidate selection

The interval spanning APOL3 to MYH9, chr22:36,140,330-36,388,018 (GRCh38,
247,688 bp), fits within one model context window. Variants catalogued by gnomAD
v4 [12] were retrieved with alternate allele frequencies for the African
genetic-ancestry group. We retained single-nucleotide substitutions with
alternate allele frequency at least 1%. After removing duplicate identifiers this
gave **1,921 unique variants**, of which 1,865 were non-coding. The approximately
99,000 non-coding variants at the locus were not all evaluated: frequency and
variant-type filtering, not model prioritisation, removed most of them.

Candidates were the 15 highest-ranked non-coding variants by absolute predicted
podocyte accessibility effect after the distance stratification in section 2.3,
excluding variants within the APOL1 gene body. No annotation, peak overlap or
prior literature informed selection. The full ranked list and the selection flow
are deposited.

### 2.2 Model, version and provenance

Predictions were generated through the AlphaGenome public API [9] using the
ALL_FOLDS model version, the default when none is specified. Channels used were a
gene-masked RNA scorer restricted to the APOL1 transcript and position-local DNase
and ATAC scorers, all differential REF-versus-ALT scorers. Kidney tracks were
named explicitly rather than matched by substring.

**This evaluation is not independent of model training.** AlphaGenome is trained
on ENCODE DNase-seq and ATAC-seq coverage, and every kidney biosample used here
corresponds to a named training experiment in the published track manifest
(ENCSR206OJJ, ENCSR785BDQ, ENCSR000EPW, ENCSR000EOL, ENCSR000EOK, ENCSR175IWT,
ENCSR297VGU). The ALL_FOLDS model is trained on all eight genomic sections and has
no held-out genomic region. Applying the published fold definitions, this locus
lies in data fold 6, which is in the training partition of the served checkpoint
and the validation partition of one other; **no released checkpoint designates it
as its test region**. We therefore describe these analyses as a consistency check
between predictions and the experimental data underlying them. Whether ENCODE-based
evaluation of this model is independent in general depends on checkpoint, region,
experiments and endpoint; the statement here concerns this case.

### 2.3 Expression channel and its control

The gene-masked scorer assigns larger predicted effects nearer the gene. We
stratified by gene body, regressed log absolute effect on log distance to the
nearer APOL1 gene-body boundary, and ranked on the residual. We asked whether the
coding risk variants remained extreme after correction. **This is a weak control.**
Coding variants are not expression-null: a coding change can alter splicing,
transcript stability or abundance, so a high predicted expression effect for G1 or
G2 is not by itself evidence of model failure. We report it as a positional
confound of unknown size and use the position-local accessibility channel for the
remaining analyses.

### 2.4 Measured accessibility, partitions and the reconciliation

ENCODE DNase-seq peak files (GRCh38, released) [13] were merged over the locus
for kidney and podocyte biosamples to form **one kidney union**, and separately for
six comparison tissues: liver, lung, heart left ventricle, stomach, spleen and
thyroid gland. Kidney peaks were then labelled **kidney-restricted** (no peak
detected in any comparison tissue) or **shared**. These labels are mutually
exclusive and their union is the kidney set by construction; the script asserts
that the base pairs and the variant memberships sum, because an earlier version of
this work compared a pooled result built on one peak set against partitions built
on another and the discrepancy was detectable only by that arithmetic.

"Kidney-restricted" means no peak was detected in the selected comparison panel,
not that the region is closed. The comparison panel was chosen to span germ layers
and to be well represented in ENCODE at this locus; it does not match the scored
outputs one for one, and sensitivity to that choice is deposited.

Enrichment was tested by Fisher exact test against the other common variants at
the same locus, with conditional maximum-likelihood odds ratios and exact 95%
intervals. **All four contingency cells and the eligible denominator are reported
for every endpoint** (Supplementary Table 1). Ranking depth 200 is the primary
analysis; depths 25, 50 and 100 are nested sensitivity analyses, not independent
replications.

Variants within one peak share the element and have correlated scores, so
comparisons between rankings resample **peaks** with replacement, with variants
outside peaks carried unchanged, 500 replicates, percentile intervals, unadjusted
for multiplicity across six comparisons. The change in contrast between the pooled
and kidney-restricted endpoints was estimated on the same resampled draws, so the
two contrasts are paired.

### 2.5 Linkage disequilibrium with G1 and G2

An earlier version reported G2 as untestable. That was incorrect. The six-base
deletion is absent from the biallelic GRCh38 lift first used but present in the
original 1000 Genomes phase 3 release at chr22:36,662,041 (GRCh37) as AATAATT>A,
carrying no identifier, which is why identifier-based searches failed [14, 15].
Identification rests on normalised position and alleles; frequency is
corroboration only (19.5% Gambian, 18.2% Mende, 0% in European and East Asian
panels).

We computed r-squared with G1a, G1b and G2 separately, in seven African
populations analysed individually, and report the maximum across populations with
the responsible population and partner identified. **D' and the maximum r-squared
attainable at the observed allele frequencies are reported alongside**, because
low r-squared between variants of different frequency does not establish
haplotypic independence. Neither does low D' establish it; both are filtering
thresholds, and candidates passing them are described as below those thresholds
rather than as independent.

### 2.6 Control locus and prior-record search

The pipeline was repeated at the beta-globin cluster, chr11:5,150,000-5,397,688,
matched on length, paralogue-cluster structure and, by frequency-stratified
subsampling, variant count. Candidates were searched against Europe PMC full text,
the GWAS Catalog [16] and GTEx v8 [17] on 19 September 2026; a positive control
confirmed retrieval of a known variant at this locus. "No prior record identified"
describes a bounded search, not proof of novelty.

### 2.7 Data and code

Analysis code, the frozen variant ledger, per-variant predictions, all contingency
counts, peak accessions, partition labels and the candidate selection flow are
deposited. No new participant data were collected; existing publicly available
1000 Genomes phased haplotypes and gnomAD summary frequencies were analysed.

## 3. Results

### 3.1 The expression channel tracks position

Variants inside the APOL1 gene body scored roughly three times higher than those
outside (median absolute effect 0.0022, n = 134, against 0.00077, n = 1,796), and
distance explained 5.7% of the variance in log absolute effect. Within the gene
body the 13 coding variants scored higher still (0.0062, against 0.0020 for the
121 non-coding intragenic variants). After residualising, the three coding
sentinel variants remained at the 88th to 98th percentile.

We do not read this as showing the channel is unusable, for the reason given in
section 2.3. We report it as a positional confound and use the accessibility
channel thereafter, where distance explained 0.3% of variance for podocyte DNase
and none for kidney ATAC.

### 3.2 The reconciled comparison

The kidney union comprised 110 peaks over 44,576 bp, of which 39 peaks
(14,686 bp) were kidney-restricted and 71 (29,890 bp) shared. Of 1,921 variants,
277 fell in the union, 84 restricted and 193 shared; the partitions sum to the
parent in both peaks and variants.

At the primary depth of 200, from that single frozen object (Table 1, Figure 2):

The podocyte ranking led in all three views, with odds ratios of 8.94 pooled, 4.28
kidney-restricted and 8.13 shared. The strongest non-renal ranking, lung, gave
5.58, 2.48 and 5.98. The whole-kidney ranking performed worst of the renal outputs
throughout, at 1.99, 0.90 and 2.46.

**Partitioning did not change the tissue contrast.** The difference between the
podocyte-minus-comparator contrast in kidney-restricted regions and the same
contrast pooled ranged from -0.10 to +0.40 in log odds ratio and excluded zero for
none of six comparators (Figure 3). Shared regions account for two-thirds of the
kidney union, and about 12% of the locus by length.

### 3.3 Two earlier conclusions were artefacts of the same error

We report these because each looked clean and neither was detectable without
checking that the sets matched.

First, ranking by predicted **podocyte** accessibility and evaluating against
**whole-kidney** peaks alone, several non-renal rankings matched or exceeded the
podocyte ranking at intermediate depths, suggesting the tissue label carried
little information. That comparison evaluates a cell-type prediction against a
tissue-level measurement.

Second, when the measured set was extended to include podocyte peaks for the
partitioned analysis but not for the pooled one, the partitioned result appeared
to reverse the pooled result. The two were different parent sets. The discrepancy
was visible only in the arithmetic: the partitions summed to 44,576 bp against
30,358 bp for the supposed whole.

Both errors reduce to mismatching the measured evaluation set to the predicted
track. Neither produced an implausible number.

### 3.4 Candidate prioritisation yields no novel independent modifier

Of 45 candidate-by-risk-variant comparisons, the maximum r-squared against any
risk variant was 0.177, which alone suggests independence. Twenty-two of the 45
combine D' at or near 1.0 with r-squared below 0.2, several at 96-100% of the
maximum their allele frequencies permit, consistent with nesting on the risk
background. Applying a threshold of D' below 0.8 against all three risk variants
leaves four candidates (Table 2).

Of those four, rs132708 is a reported GTEx eQTL for APOL4 and rs713797 for
FOXRED2, with three plasma protein associations catalogued; the specific records
are deposited. rs4820232 is correlated with rs5750250 (r-squared 0.55, D' 1.0),
one of several strongly associated variants in the MYH9 intron 13-15 region
identified by dense mapping [18]; that is a historical regional association, not a
demonstrated regulatory effect or an assigned target gene. rs6000250 falls below
both filtering thresholds and no prior record was identified.

An eQTL association does not establish that the indexed variant is causal, and an
association with APOL4 or FOXRED2 does not exclude an additional effect on APOL1.
A variant carried on a risk haplotype could still modify penetrance among
carriers. The defensible statement is therefore narrow: **this analysis does not
establish a novel, independent APOL1 modifier.**

Most candidates lie within MYH9. The explanation is mundane: the ranking channel
is position-local and gene-agnostic, MYH9 occupies much of the interval, and
accessible chromatin is denser there. It does mean "APOL1 locus" misdescribes what
was prioritised.

### 3.5 The control locus is too imprecise to exclude enrichment

At the beta-globin cluster, 2 of the top 200 variants fell in measured kidney
chromatin against 12 of 1,718: conditional odds ratio 1.44, one-sided p = 0.44,
exact 95% CI 0.155-6.520. No significant enrichment was detected, and the interval
includes the test-locus estimate at this and every other depth. **The control
estimate is too imprecise to exclude enrichment of the observed magnitude.** The
control locus also carries far less kidney chromatin, 1.0% of its length against
12.3%.

## 4. Discussion

The substantive finding of this study is about its own conduct. We reached two
opposite conclusions about whether a tissue-labelled ranking carries tissue
information, and both were artefacts of evaluating a prediction against a measured
set that did not match it. The first evaluated a podocyte prediction against
whole-kidney measurement. The second compared analyses built on two different peak
sets. Each produced a coherent result with a plausible mechanism, and neither was
visible except by checking that the sets summed.

When the analysis is reconciled onto one frozen variant set and one measured
kidney union, the picture is simpler than either earlier account. The
podocyte-labelled ranking is enriched in measured kidney accessible regions in
every view, most strongly in shared regions and still clearly in kidney-restricted
ones, and partitioning does not measurably change its advantage over non-renal
rankings. We had predicted that partitioning would matter and it did not.

Three limitations bound what this can mean, and the first is decisive.

The evaluation is not independent of training. Every kidney track used is a named
training experiment and the locus lies in a training fold of the served
checkpoint. What is demonstrated is that predictions are consistent with the data
underlying them, in a way that discriminates tissue. A test of generalisation
requires a held-out region and a checkpoint that excludes it, and none of the
released checkpoints designates this region as its test partition. A
validation-fold sensitivity analysis would be a weaker but distinguishable
alternative and is not attempted here.

The endpoint is positional. A differential REF-versus-ALT score ranked against
peak membership tests whether predictions concentrate at accessible positions. It
does not test which allele changes accessibility, nor the direction or magnitude
of any effect, and we have not compared it against a simple reference-accessibility
ranking to establish what the differential scorer adds.

The candidate list does not support a biological claim. Most candidates were
nested on the risk background once D' was examined, three of the four survivors
have published associations pointing elsewhere, and none has functional evidence.

For investigators applying these models at kidney loci, three practical points
follow. Match the measured evaluation set to the predicted track, at the same
level of cell-type resolution. Verify that partitions sum to their parent set
before interpreting any difference between them; that single check would have
caught both of our errors immediately. And report the model version and the fold
containing the locus, because the default served checkpoint holds nothing out.

What this study does not establish should be explicit. No candidate is shown to
alter expression of any gene. No target gene is assigned. No allele-specific
activity is demonstrated. No association with kidney disease is tested, and
nothing here bears on prognosis, on which carriers develop disease, or on
treatment.

## Data availability

Analysis code, the frozen variant ledger, per-variant predictions, all contingency
counts, ENCODE file accessions, partition labels, the candidate selection flow and
the interval-to-fold mapping are available at
https://github.com/PatientThread/apol1-alphagenome, with a tagged release for
review and a persistent identifier on publication. Model predictions were obtained
under AlphaGenome's non-commercial research terms and are redistributed subject to
them; model weights were not accessed.

## Acknowledgements

This work used data generated by the ENCODE Consortium, gnomAD and the 1000
Genomes Project. I thank the participants and investigators of those projects.

## Author contributions

CL is the sole author and conceived the study, wrote the analysis code, performed
the analyses, prepared the figures and wrote the manuscript.

## Funding

No funding was received in support of this study.

## Conflict of interest statement

None declared.

## Ethical approval

Ethical approval was not required. No new participant data were collected. The
study analysed existing publicly available 1000 Genomes phased haplotypes and
gnomAD summary frequencies, together with ENCODE consortium data and model
predictions.

## References

1. Genovese G, Friedman DJ, Ross MD, et al. Association of trypanolytic ApoL1 variants with kidney disease in African Americans. Science 2010; 329: 841-845
2. Kopp JB, Smith MW, Nelson GW, et al. MYH9 is a major-effect risk gene for focal segmental glomerulosclerosis. Nat Genet 2008; 40: 1175-1184
3. Beckerman P, Bi-Karchin J, Park AS, et al. Transgenic expression of human APOL1 risk variants in podocytes induces kidney disease in mice. Nat Med 2017; 23: 429-438
4. Ojo AO, Adu D, Bramham K, et al. APOL1 kidney disease: conclusions from a Kidney Disease: Improving Global Outcomes (KDIGO) Controversies Conference. Kidney Int 2025; 108: 763-779
5. Adamson WE, Noyes H, Ogunsola J, Parekh RS, Cooper A, MacLeod A. Coding, modifier, and regulatory effects shape circulating APOL1 levels. Hum Mol Genet 2026; 35: ddag087
6. Li Y, Bozack AK, Schlosser P, et al. APOL1 risk genotypes influence DNA methylation across multiple genomic elements in APOL1-APOL4-MYH9 region in African Americans. Clin Epigenetics 2026; 18: 147
7. Sheng X, Guan Y, Ma Z, et al. Mapping the genetic architecture of human traits to cell types in the kidney identifies mechanisms of disease and potential treatments. Nat Genet 2021; 53: 1322-1333
8. Linder J, Srivastava D, Yuan H, et al. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. Nat Genet 2025; 57: 949-961
9. Avsec Ž, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature 2026; 649: 1206-1218
10. Karollus A, Mauermeier T, Gagneur J. Current sequence-based models capture gene expression determinants in promoters but mostly ignore distal enhancers. Genome Biol 2023; 24: 56
11. Kathail P, Shuai RW, Chung R, Ye CJ, Loeb GB, Ioannidis NM. Current genomic deep learning models display decreased performance in cell type-specific accessible regions. Genome Biol 2024; 25: 202
12. Chen S, Francioli LC, Goodrich JK, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature 2024; 625: 92-100
13. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 2020; 583: 699-710
14. Zhang Z, Hao K, Ross MJ, Murphy B, Menon MC. APOL1 G2 risk allele-clarifying nomenclature. Kidney Int 2017; 92: 518-519
15. 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature 2015; 526: 68-74
16. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Res 2023; 51: D977-D985
17. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020; 369: 1318-1330
18. Nelson GW, Freedman BI, Bowden DW, et al. Dense mapping of MYH9 localizes the strongest kidney disease associations to the region of introns 13 to 15. Hum Mol Genet 2010; 19: 1805-1815

## Figure legends

**Figure 1.** Locus map. APOL1, APOL4 and MYH9 with transcription orientation, the
positions of G1 and G2, the four candidates passing linkage filtering, and
measured accessible regions drawn on separate rows for kidney-restricted and
shared, with the comparison-tissue union below.

**Figure 2.** Enrichment of the top 200 variants in measured accessible regions,
by the tissue output used to rank, for the pooled kidney union and each partition.
Conditional odds ratios with exact 95% intervals against the other common variants
at the same locus; eligible n = 1,921. Intervals assume variant-level
independence and are therefore narrower than the peak-resampled contrasts in
Figure 3. Kidney rankings are shown in black.

**Figure 3.** Change in tissue contrast on partitioning. For each comparator, the
podocyte-minus-comparator contrast in kidney-restricted regions minus the same
contrast pooled, in natural-log odds ratio, estimated on identical resampled draws
of peaks. Points are observed values; intervals are percentile bootstrap, 500
replicates, unadjusted for six comparisons. No interval excludes zero.

**Figure 4.** The control locus. Conditional odds ratios with exact 95% intervals
at the beta-globin cluster against the test-locus estimate, same pipeline and peak
source. Zero-event groups are labelled rather than plotted as finite estimates.

## Tables

**Table 1.** Enrichment at the primary ranking depth of 200, by tissue output and
endpoint, from one frozen variant set (n = 1,921) and one measured kidney union.
Conditional odds ratios; all contingency cells are in Supplementary Table 1.
Kidney-restricted and shared sum to the pooled kidney union: 39 + 71 = 110 peaks,
14,686 + 29,890 = 44,576 bp, 84 + 193 = 277 variants.

| Ranked by | Pooled kidney | Kidney-restricted | Shared |
|---|---|---|---|
| Podocyte (kidney) | 8.94 | 4.28 | 8.13 |
| Lung | 5.58 | 2.48 | 5.98 |
| Hepatocyte | 5.02 | 1.94 | 5.80 |
| Stomach | 4.63 | 1.46 | 5.80 |
| Brain | 3.72 | 1.94 | 4.07 |
| Liver | 2.12 | 0.90 | 2.66 |
| Kidney (whole) | 1.99 | 0.90 | 2.46 |

**Table 2.** The four candidates passing linkage filtering. Alternate allele
frequency is for the African genetic-ancestry group of gnomAD v4. r-squared and D'
are maxima against G1a, G1b or G2 across seven African populations analysed
separately; maxima for the two measures may arise from different comparisons and
the responsible population and partner for each are deposited.

| Variant | Alt AF | max r² | max D' | Nearest gene | Prior record identified |
|---|---|---|---|---|---|
| rs713797 | 0.396 | 0.070 | 0.38 | MYH9 | reported FOXRED2 eQTL; three plasma protein associations |
| rs132708 | 0.839 | 0.062 | 0.66 | APOL4 | reported APOL4 eQTL, multiple tissues |
| rs6000250 | 0.111 | 0.177 | 0.69 | MYH9 | none identified by the searches described |
| rs4820232 | 0.301 | 0.106 | 0.79 | MYH9 | correlated with rs5750250, a historical regional association marker |

Eleven further candidates were excluded on D' at or near 1.0 with G1 or G2 despite
low r-squared, consistent with haplotypic nesting rather than independence.

**Supplementary Table 1.** All contingency counts for every endpoint, ranking and
depth: top-set hits and misses, background hits and misses, and the eligible
denominator. Deposited as `results/reconciled_counts.csv`.
