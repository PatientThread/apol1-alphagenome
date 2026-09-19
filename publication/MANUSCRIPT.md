<!-- DRAFT 3. Rewritten after editorial review. THE CENTRAL CLAIM HAS INVERTED.
     Draft 2 said the tissue label carries no usable information. Partitioning
     measured accessible regions into tissue-specific and shared shows that it
     does, and that the earlier endpoint could not see it. The paper is now about
     benchmark design. Two further corrections: validation is NOT independent,
     the model was trained on these exact tracks at a locus in its training fold;
     and the control locus excludes nothing. The candidate list is demoted to an
     illustration. No claim about prognosis, treatment or causality survives. -->

# Benchmark design determines whether a sequence model appears tissue-blind: a worked example at the APOL1-MYH9 locus

**Christopher Lawrence**
Consultant Nephrologist, London, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net

**Running head:** Benchmark design and apparent tissue-blindness

**Keywords:** APOL1; chronic kidney disease; gene expression regulation;
machine learning; MYH9

---

## Abstract

**Background and hypothesis.** Sequence-to-function models are increasingly used
to prioritise non-coding variants at kidney disease loci, and their tissue-labelled
outputs are read as tissue-specific. We asked whether that reading is justified at
the APOL1-MYH9 locus, and found that the answer depends almost entirely on how the
benchmark is constructed.

**Methods.** We scored 1,930 variants common in African genetic-ancestry
populations across chr22:36,140,330-36,388,018 with AlphaGenome, and tested
rankings for enrichment in measured ENCODE kidney open chromatin against the other
common variants at the same locus. We then repeated the test after partitioning
measured accessible regions into those open in kidney but closed in six comparison
tissues, those shared, and those specific to other tissues. Rankings were compared
pairwise with peaks, not variants, as the resampling unit. Linkage disequilibrium
with G1 and G2 was computed directly from phased haplotypes in seven African
populations, reporting D' and the maximum r-squared the allele frequencies permit
alongside r-squared.

**Results.** Pooling all measured kidney peaks, a hepatocyte ranking identified
them almost as well as a kidney ranking, which suggests the tissue label is
uninformative. Partitioning reverses this. In kidney-specific regions the podocyte
ranking reached an odds ratio of 4.28 (95% CI 2.52-7.11) against 2.48 for the best
non-kidney ranking, and exceeded five of six alternatives with paired intervals
excluding zero. In regions specific to other tissues it was the weakest of seven.
Shared regions, where all rankings perform similarly, dominate the locus and
therefore dominate the pooled result.

**Conclusions.** The model does encode tissue identity; the original benchmark
could not detect it. Benchmarks that pool tissue-specific with shared regulatory
regions will understate tissue discrimination. We also show that ENCODE-based
validation of this model is not independent, and that candidate prioritisation at
this locus survives neither linkage nor prior-literature scrutiny.

---

## Key learning points

**What was known**

- Two coding haplotypes at APOL1, G1 and G2, carry most of the excess kidney
  disease risk in people of recent African ancestry, but only a minority of
  two-allele carriers develop disease.
- Sequence-to-function models are proposed for interpreting the non-coding
  variation at such loci, and their tissue-labelled outputs are widely read as
  carrying tissue-specific information.
- Regulatory variation at this locus is under active investigation, including
  variants associated with circulating APOL1 and methylation differences across
  the APOL1-APOL4-MYH9 region.

**This study adds**

- Whether a model appears tissue-blind depends on benchmark construction: pooling
  tissue-specific with shared regulatory regions hid a discrimination that is
  clear once the two are separated.
- Validation of this model against ENCODE data is not independent; the tracks used
  here are named training experiments and the locus lies in a training fold.
- Candidate prioritisation at this locus did not survive scrutiny: low r-squared
  with G1 and G2 reflected allele-frequency differences rather than haplotypic
  independence for most candidates.

**Potential impact**

- Investigators benchmarking these models in kidney should partition measured
  regulatory regions by tissue restriction before concluding anything about
  tissue specificity.
- Model provenance, specifically the version used and the fold containing the
  locus, should be reported as routinely as software versions.
- Nothing here identifies a modifier of APOL1 penetrance, a target gene or a
  clinically actionable variant, and none is claimed.

---

## 1. Introduction

Two coding haplotypes at APOL1, G1 and G2, account for much of the excess risk of
non-diabetic kidney disease in people of recent African genetic ancestry [1]. The
association was first mapped to the adjacent gene MYH9 [2] and resolved to APOL1
two years later, with the MYH9 signal attributed to linkage disequilibrium [1].

Penetrance is incomplete, and because the mechanism is dose-dependent, genetic
modifiers of APOL1 expression are a plausible explanation [3]. Regulatory variation
at this locus is under active investigation: cis-regulatory variants have been
associated with circulating APOL1 concentration [4], and APOL1 risk genotypes are
associated with methylation differences at enhancers and promoters across the
APOL1-APOL4-MYH9 region [5]. Kidney expression resources have improved, including
maps from 659 microdissected kidney samples [6], though ancestry representation and
variant coverage remain limiting.

Sequence-to-function models predict regulatory activity from DNA and are proposed
for prioritising non-coding variants [7, 8]. Earlier models of this family perform
less well in regions that distinguish cell types [9, 10], and a companion analysis
found that kidney-labelled outputs predicted kidney expression effects no better
than outputs labelled for other tissues, leaving open whether that reflected
genuinely shared regulatory effects or limited use of tissue identity [11].

We set out to prioritise regulatory candidates at this locus. That attempt largely
failed, for reasons reported here. In the course of testing it, we found that the
apparent answer to the tissue-specificity question depends on how the benchmark is
built, and that is the finding we consider most useful.

## 2. Materials and methods

### 2.1 Locus, variants and selection

The interval spanning APOL3 to MYH9, chr22:36,140,330-36,388,018 (GRCh38,
247,688 bp), fits inside one model context window. All variants catalogued by
gnomAD v4 were retrieved with African genetic-ancestry group alternate allele
frequencies. We retained single-nucleotide substitutions with frequency at least
1%, giving 1,930 variants of which 1,865 were non-coding. The approximately 99,000
non-coding variants at the locus were therefore not all evaluated; frequency and
variant-type filtering, not model prioritisation, removed most of them.

### 2.2 Model, version and provenance

Predictions were generated through the AlphaGenome public API [8] using the
ALL_FOLDS model version, which is the default when no version is specified.
Two channels were used: a gene-masked RNA scorer restricted to the APOL1
transcript, and position-local DNase and ATAC scorers. Kidney tracks were named
explicitly rather than matched by substring.

**This is not independent validation, and we do not describe it as such.**
AlphaGenome is trained on ENCODE DNase-seq and ATAC-seq coverage. Every kidney
biosample used here corresponds to a named training experiment in the model's
published track manifest (ENCSR206OJJ, ENCSR785BDQ, ENCSR000EPW, ENCSR000EOL,
ENCSR000EOK, ENCSR175IWT, ENCSR297VGU). Moreover, the ALL_FOLDS model is trained
on all eight genomic sections and has no held-out genomic regions; the developers
state that such models lack a dedicated holdout and were not used for track
prediction evaluations. Applying the published fold definitions, this locus lies
in fold 6, which is in the training partition of every released checkpoint except
one, where it is the validation partition. It is the test partition of none. We
therefore describe the enrichment analyses as a consistency check between
predictions and the experimental data underlying them.

### 2.3 Expression channel and its control

The gene-masked scorer assigns larger predicted effects nearer the gene. We
stratified by gene body, regressed log absolute effect on log distance and ranked
on the residual. As a control we asked whether the coding risk variants remained
extreme after correction. **This control is weaker than we first treated it as.**
Coding variants are not expression-null: a coding change can alter splicing,
transcript stability or abundance, so a high predicted expression effect for G1 or
G2 is not by itself evidence of model failure. We report the raw and adjusted
results as sensitivity analyses rather than as a pass-or-fail test.

### 2.4 Measured accessibility and the partition

ENCODE DNase-seq peak files (GRCh38, released) were merged over the locus for
kidney and podocyte, and separately for six comparison tissues: liver, lung, heart
left ventricle, stomach, spleen and thyroid gland. Measured accessible regions were
then partitioned into **kidney-specific** (open in kidney, closed in all comparison
tissues), **shared** (open in kidney and at least one comparison tissue) and
**other-tissue-specific**.

Enrichment of top-ranked variants was tested by Fisher exact test against the other
common variants at the same locus, with conditional maximum-likelihood odds ratios
and exact 95% intervals. We report all 2x2 counts. The same-locus background
matches sequence composition and gene density only approximately, since these vary
within the interval.

Variants inside one peak share the element and have highly correlated scores, so
paired comparisons between rankings resample **peaks**, not variants, with 600
replicates. The four ranking depths are nested sensitivity analyses, not
independent replications; a depth of 200 is the primary analysis.

### 2.5 Linkage disequilibrium with G1 and G2

Earlier analysis reported that G2 was untestable. That was incorrect. The six-base
deletion is absent from the biallelic GRCh38 lift we first used, but present in the
original 1000 Genomes phase 3 release at chr22:36,662,041 (GRCh37) as AATAATT>A,
carrying no identifier, which is why identifier-based searches failed [12]. Its
frequency confirms identity: 19.5% in Gambian, 18.2% in Mende, 0% in European and
East Asian panels.

We computed r-squared with G1 and G2 separately, in seven African populations
analysed individually, from phased haplotypes. **We also report D' and the maximum
r-squared attainable given the two allele frequencies**, because low r-squared
between variants of different frequency does not establish haplotypic independence:
at frequencies of 0.026 and 0.20, r-squared cannot exceed about 0.107 even under
complete nesting.

### 2.6 Control locus and prior-record search

The pipeline was repeated at the beta-globin cluster, chr11:5,150,000-5,397,688,
matched on length, paralogue-cluster structure and, by frequency-stratified
subsampling, variant count. Candidates were searched against Europe PMC full text,
the GWAS Catalog [13] and GTEx v8 [14]; a positive control confirmed that full-text
search reaches variants of this class.

### 2.7 Data and code

Analysis code, the frozen variant set, all per-variant predictions, raw overlap
counts and the selection flow are deposited (below). This study used publicly
available, de-identified genotype and haplotype data from 1000 Genomes and gnomAD,
and no new individual-level data were generated or accessed.

## 3. Results

### 3.1 The expression channel tracks position

Variants inside the APOL1 gene body scored roughly three times higher than those
outside (median absolute effect 0.0022, n = 134, against 0.00077, n = 1,796), and
distance alone explained 5.7% of the variance in log absolute effect. Within the
gene body the 13 coding variants scored higher still (0.0062). After residualising
on distance the coding risk variants remained in the upper tail, at the 88th to
98th percentile.

We do not read this as demonstrating that the channel is unusable. Coding variants
can genuinely affect expression, and 5.7% of variance is a modest positional
component. We report it as a confound of unknown size and use the position-local
accessibility channel for the remaining analyses, which do not carry it: distance
explained 0.3% of variance for podocyte DNase and none for kidney ATAC.

### 3.2 Pooled, the ranking appears tissue-blind

Measured kidney peaks covered 30,358 bp, 12.3% of the locus. Ranking by predicted
podocyte accessibility, 69 of the top 200 variants (34.5%) fell inside them
against 7.9% of the remainder: odds ratio 6.13.

Repeating this with eight deliberately non-renal tissue outputs gave closely
similar results (Table 1). A hepatocyte ranking reached 6.69 at a depth of 100
against 5.23 for podocyte. On this endpoint the tissue label appears to carry
little information.

### 3.3 Partitioned, the same data show clear tissue discrimination

Measured accessible regions divided into 39 kidney-specific peaks (14,686 bp), 71
shared (29,890 bp) and 47 other-tissue-specific (16,354 bp). Shared regions are
the largest partition.

In **kidney-specific** regions the podocyte ranking reached an odds ratio of 4.28
(2.52-7.11), against 2.48 (1.35-4.34) for the best non-kidney ranking (Table 2,
Figure 2). Paired comparison on the same variants, resampling peaks, gave
differences in log odds ratio excluding zero against five of six alternatives:
liver +1.66 (0.75-2.51), whole kidney +1.58 (0.91-2.34), stomach +1.04
(0.41-1.66), hepatocyte +0.84 (0.53-1.15) and brain +0.79 (0.17-1.39). Only lung
was not significantly exceeded, +0.54 (-0.17 to +1.14).

The mirror control agrees. In **other-tissue-specific** regions the podocyte
ranking was the weakest of all seven (1.51, 0.73-2.86) while stomach, hepatocyte
and lung led (3.42, 3.19, 2.97).

In **shared** regions every ranking was strongly enriched and the differences
narrowed: podocyte 8.14, lung 5.98, hepatocyte 5.80.

Because shared regions dominate the locus, an endpoint that pools all measured
kidney peaks is largely measuring shared accessibility. That is why the pooled
analysis in section 3.2 suggested the tissue label was uninformative when the
partitioned analysis shows it is not.

### 3.4 Candidate prioritisation does not survive scrutiny

Fifteen candidates were taken forward from the accessibility ranking. Testing
against G1 and G2 directly, the maximum r-squared against any risk variant was
0.177, which taken alone suggests independence. It does not. Of 45 comparisons, 22
combine D' at or near 1.0 with r-squared below 0.2, and several sit at 96-100% of
the maximum their allele frequencies permit, indicating complete haplotypic nesting
on the risk background. Requiring D' below 0.8 against all three risk variants
leaves four candidates: rs4820232, rs6000250, rs132708 and rs713797 (Table 3).

Three of those four already have published regulatory records, and none points at
APOL1. rs132708 is a GTEx eQTL for APOL4; rs713797 is a catalogued plasma protein
quantitative trait locus and an eQTL for FOXRED2; rs4820232 is correlated with
rs5750250 (r-squared 0.55, D' 1.0), among the strongest associations in the MYH9
intron 13-15 region identified by dense mapping [15]. An eQTL association does not
establish that the indexed variant is causal, and an association with APOL4 or
FOXRED2 does not exclude an additional effect on APOL1.

One candidate, rs6000250, is both haplotypically independent of the risk variants
and without prior record. We regard that as an illustration of what the workflow
yields rather than as a finding, and no functional evidence supports it.

### 3.5 The control locus excludes nothing

At the beta-globin cluster, 2 of the top 200 variants fell in measured kidney
chromatin against 12 of 1,718: odds ratio 1.44, one-sided p = 0.44, exact 95% CI
0.155-6.520. No significant enrichment was detected. **That interval includes the
test-locus estimate of 6.13, as does the interval at every other depth**, so this
control does not exclude an effect of the size observed at APOL1-MYH9. It is an
absence of evidence and we report it as such. The control locus also carries much
less kidney chromatin (1.0% of its length against 12.3%), so a single interval of
this width cannot establish that the pipeline never produces spurious enrichment.

## 4. Discussion

The useful result here is about benchmark construction rather than about this
model. Pooling every measured kidney peak, a hepatocyte ranking recovered them
almost as well as a kidney ranking, and we initially concluded that the tissue
label was close to decorative. Separating regions that are open only in kidney
from those open in many tissues reverses that conclusion: the podocyte ranking
leads clearly where tissue identity is the whole signal, and is the weakest of
seven where it should be. Shared regions dominate, so they dominated the pooled
endpoint and hid a discrimination that was plainly present.

The general point is that a benchmark which does not separate tissue-specific from
shared regulatory regions will understate tissue discrimination, possibly to the
point of concluding that a model has none. We made that error before correcting it,
and it is an easy one to make because the pooled analysis looks reasonable and
returns a clean answer. This also bears on our companion analysis, which left open
whether similar cross-tissue performance reflected genuinely shared effects or
limited use of tissue identity [11]. The evidence here favours the former.

Two limitations bound everything above, and neither is incidental.

First, none of this is independent validation. AlphaGenome is trained on ENCODE
accessibility data; every kidney track used here is a named training experiment;
and the locus sits in a training fold of the served model. What we demonstrate is
that predictions are consistent with the data underlying them, in a way that
discriminates tissue. A genuine test of generalisation would require a held-out
locus and a checkpoint that excludes it, and no released checkpoint holds this
locus out. We report model version, fold assignment and track accessions so that
the limitation is auditable rather than implicit.

Second, the candidate list does not survive. Our first linkage analysis used
r-squared alone and called seven candidates independent of the risk haplotypes.
Reporting D' and the maximum r-squared the allele frequencies permit shows that
most of those were completely nested on the risk background, and that low
r-squared reflected a frequency difference. Four survive, three of those are
already described with non-APOL1 targets, and one remains. We therefore make no
claim to have identified regulatory modifiers of APOL1.

The awkward observation deserves stating plainly rather than as a footnote. Most
candidates lie within MYH9, the gene to which this association was first mapped
[2] and from which it was moved in 2010 [1]. We are not reopening that question,
and nothing here indicates that MYH9 rather than APOL1 explains the disease
association. The explanation is mundane: the ranking channel is position-local and
gene-agnostic, MYH9 occupies much of the interval, and accessible chromatin is
denser there. It does mean that "APOL1 locus" is the wrong description of what was
prioritised.

What this study does not establish should be equally explicit. No candidate is
shown to alter expression of any gene. No target gene is assigned. No
tissue-specific activity is demonstrated for any individual variant, as distinct
from the ranking as a whole. No association with kidney disease is tested, and
nothing here bears on prognosis, on which carriers develop disease, or on
treatment.

For investigators applying these models at kidney loci we would suggest three
things. Partition measured regulatory regions by tissue restriction before drawing
conclusions about tissue specificity. Report the model version and the fold
containing the locus, because the default served model holds nothing out. And test
linkage with D' and the attainable maximum, not r-squared alone, because allele
frequency differences will otherwise manufacture apparent independence.

## Data availability

Analysis code, the frozen variant set, per-variant predictions, all raw overlap
counts, the candidate selection flow and the fold assignment of the locus are
available at https://github.com/PatientThread/apol1-alphagenome and will be
deposited with a persistent identifier on publication. Variant frequencies are from
gnomAD v4; peak calls from the ENCODE portal [16] with accessions listed;
haplotypes from 1000 Genomes phase 3 [17]. Model predictions were obtained under
AlphaGenome's non-commercial research terms and are redistributed subject to them;
model weights were not accessed.

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

Ethical approval was not required. This study analysed publicly available,
de-identified, summary-level and genotype data from open consortia. No new
individual-level data were generated or accessed and no participants were
recruited.

## References

1. Genovese G, Friedman DJ, Ross MD, et al. Association of trypanolytic ApoL1 variants with kidney disease in African Americans. Science 2010; 329: 841-845
2. Kopp JB, Smith MW, Nelson GW, et al. MYH9 is a major-effect risk gene for focal segmental glomerulosclerosis. Nat Genet 2008; 40: 1175-1184
3. Ojo AO, Adu D, Bramham K, et al. APOL1 kidney disease: conclusions from a Kidney Disease: Improving Global Outcomes (KDIGO) Controversies Conference. Kidney Int 2025; 108: 763-779
4. Adamson WE, Noyes H, Ogunsola J, Parekh RS, Cooper A, MacLeod A. Coding, modifier, and regulatory effects shape circulating APOL1 levels. Hum Mol Genet 2026; 35: ddag087
5. Li Y, Bozack AK, Schlosser P, et al. APOL1 risk genotypes influence DNA methylation across multiple genomic elements in APOL1-APOL4-MYH9 region in African Americans. Clin Epigenetics 2026; 18: 147
6. Sheng X, Guan Y, Ma Z, et al. Mapping the genetic architecture of human traits to cell types in the kidney identifies mechanisms of disease and potential treatments. Nat Genet 2021; 53: 1322-1333
7. Linder J, Srivastava D, Yuan H, et al. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. Nat Genet 2025; 57: 949-961
8. Avsec Ž, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature 2026; 649: 1206-1218
9. Karollus A, Mauermeier T, Gagneur J. Current sequence-based models capture gene expression determinants in promoters but mostly ignore distal enhancers. Genome Biol 2023; 24: 56
10. Kathail P, Shuai RW, Chung R, et al. Current genomic deep learning models display decreased performance in cell type-specific accessible regions. Genome Biol 2024; 25: 202
11. Lawrence C. Kidney representation and uncertainty in AlphaGenome regulatory variant prediction. Submitted manuscript, available from the author
12. Zhang Z, Hao K, Ross MJ, Murphy B, Menon MC. APOL1 G2 risk allele-clarifying nomenclature. Kidney Int 2017; 92: 518-519
13. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Res 2023; 51: D977-D985
14. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020; 369: 1318-1330
15. Nelson GW, Freedman BI, Bowden DW, et al. Dense mapping of MYH9 localizes the strongest kidney disease associations to the region of introns 13 to 15. Hum Mol Genet 2010; 19: 1805-1815
16. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 2020; 583: 699-710
17. 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature 2015; 526: 68-74

## Figure legends

**Figure 1.** Enrichment of top-ranked variants in measured accessible regions, by
the tissue output used to rank and by partition. Points are conditional odds
ratios with exact 95% intervals; the background is the other common variants at
the same locus. Kidney rankings are shown in black. The podocyte ranking leads in
the kidney-specific partition and is the weakest in the other-tissue-specific
partition, while rankings differ much less in the shared partition.

**Figure 2.** Paired differences in log odds ratio between the podocyte ranking
and each alternative, within the kidney-specific partition, with peaks resampled
rather than variants. Intervals exclude zero for five of six comparisons.

**Figure 3.** The expression channel. Predicted absolute effect on APOL1
expression by position, with the group median marked, and the distribution of the
residual after regressing log effect on log distance with the three coding risk
variants marked.

**Figure 4.** The control locus. Conditional odds ratios with exact 95% intervals
at the beta-globin cluster against the test locus, same pipeline and peak source.
The control intervals include the test-locus estimate at every depth.

## Tables

**Table 1.** Enrichment in all measured kidney open chromatin, pooled, by the
tissue output used to rank the same 1,930 variants. Conditional odds ratios at
four nested ranking depths; depth 200 is the primary analysis and the others are
sensitivity analyses. Raw counts are deposited.

| Ranked by | Top 25 | Top 50 | Top 100 | Top 200 |
|---|---|---|---|---|
| Podocyte (kidney) | 4.88 | 4.19 | 5.23 | 6.14 |
| Kidney (whole) | 4.06 | 3.79 | 3.04 | 2.27 |
| Hepatocyte | 4.88 | 5.07 | 6.69 | 4.81 |
| Lung | 4.06 | 4.19 | 5.23 | 5.44 |
| Frontal cortex | 4.88 | 4.19 | 4.25 | 3.85 |
| Brain | 3.33 | 2.75 | 4.48 | 4.37 |
| Stomach | 4.06 | 3.42 | 4.48 | 4.24 |
| Heart, left ventricle | 3.33 | 3.79 | 3.82 | 3.15 |
| Liver | 4.06 | 2.75 | 3.22 | 2.36 |
| Lung (left) | 1.61 | 2.44 | 2.37 | 1.49 |

**Table 2.** Enrichment by partition, top 200 by each tissue output. Conditional
odds ratios with exact 95% intervals, and the number of the top 200 falling in the
partition. Partition sizes: kidney-specific 39 peaks, shared 71,
other-tissue-specific 47.

| Ranked by | Kidney-specific | Shared | Other-tissue-specific |
|---|---|---|---|
| Podocyte (kidney) | **4.28 (2.52-7.11)**, n=26 | 8.14 (5.69-11.61), n=75 | 1.51 (0.73-2.86), n=12 |
| Lung | 2.48 (1.35-4.34), n=18 | 5.98 (4.16-8.57), n=65 | 2.97 (1.66-5.12), n=20 |
| Hepatocyte | 1.94 (1.01-3.51), n=15 | 5.80 (4.03-8.31), n=64 | 3.19 (1.80-5.46), n=21 |
| Brain | 1.94 (1.01-3.51), n=15 | 4.07 (2.78-5.89), n=53 | 2.37 (1.27-4.19), n=17 |
| Stomach | 1.46 (0.71-2.78), n=12 | 5.80 (4.03-8.31), n=64 | **3.42 (1.95-5.82)**, n=22 |
| Kidney (whole) | 0.90 (0.37-1.91), n=8 | 2.46 (1.63-3.67), n=39 | 2.37 (1.27-4.19), n=17 |
| Liver | 0.90 (0.37-1.91), n=8 | 2.66 (1.77-3.94), n=41 | 2.37 (1.27-4.19), n=17 |

**Table 3.** The four candidates surviving linkage filtering, with prior published
record. Alternate allele frequency is in the African genetic-ancestry group of
gnomAD v4. r-squared and D' are the maximum against G1a, G1b or G2 across seven
African populations analysed separately.

| Variant | Alt AF | max r² | max D' | Nearest gene | Prior published record |
|---|---|---|---|---|---|
| rs713797 | 0.396 | 0.070 | 0.38 | MYH9 | FOXRED2 eQTL; three plasma protein associations |
| rs132708 | 0.839 | 0.062 | 0.66 | APOL4 | APOL4 eQTL across multiple tissues |
| rs6000250 | 0.111 | 0.177 | 0.69 | MYH9 | none identified |
| rs4820232 | 0.301 | 0.106 | 0.79 | MYH9 | correlated with rs5750250 (r² 0.55, D' 1.0) |

Eleven further candidates were excluded: D' at or near 1.0 with G1 or G2 despite
low r-squared, indicating haplotypic nesting rather than independence.
