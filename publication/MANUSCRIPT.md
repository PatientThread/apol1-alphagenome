<!-- DRAFT 1. Written for Nephrology Dialysis Transplantation, Original Article.
     Framing decided before drafting: the primary finding is METHODOLOGICAL. The
     candidate list is secondary and heavily hedged. The MYH9 detail is faced in
     the Results and again in the Discussion rather than buried. No claim is made
     about prognosis or treatment, because none is supported: no candidate has
     been shown to change expression and no disease association is tested. -->

# What a sequence-to-function model can and cannot resolve at the APOL1-MYH9 locus

**Christopher Lawrence** BSc, MBBS, LLM, MD(Res), FRCP
Consultant Nephrologist, 9 Harley Street, London W1G 9QY, United Kingdom
ORCID 0000-0002-8159-0879
Correspondence: Christopher.lawrence3@nhs.net

**Keywords:** APOL1; chronic kidney disease; gene expression regulation;
machine learning; MYH9

---

## Abstract

**Background and hypothesis.** Kidney disease risk at the APOL1 locus is
attributed to two coding haplotypes, yet 91.8% of catalogued variation there is
non-coding and only a minority of two-allele carriers develop disease. Regulatory
variants altering APOL1 expression dose are a plausible modifier.
Sequence-to-function models predict regulatory activity from DNA alone and are
proposed for this task. We asked what such a model can resolve here, testing
every answer against data it never saw.

**Methods.** We scored 1,930 variants common in African populations across
chr22:36,140,330-36,388,018 with AlphaGenome, using a gene-masked expression
scorer and position-local accessibility scorers. Rankings were tested for
enrichment in independently measured ENCODE kidney chromatin against a
background of the other common variants at the same locus. Three controls were
pre-specified: linkage disequilibrium against all common APOL1 gene-body
variation using phased haplotypes in seven African populations kept separate; the
same ranking performed with eight deliberately non-renal tissue outputs; and the
identical pipeline at a length-matched beta-globin control locus.

**Results.** The expression scorer rediscovered proximity to the gene, failed its
pre-specified negative control and was discarded. The accessibility channel was
unconfounded and its top-ranked variants were enriched in measured kidney
chromatin (odds ratio 8.14 for the top 25 by kidney ATAC, p = 2.4x10-6). All
three controls then constrained that result. Seven of fifteen candidates were
correlated with gene-body haplotype structure. Enrichment was no stronger in
podocytes than in other kidney cell types, and no stronger for kidney than for
hepatocyte, lung or cortex outputs. At the control locus the ranking found
nothing.

**Conclusions.** The model detects regulatory positions, and detects them where
they exist, but the tissue label on the request carries no usable information.
Predictions labelled kidney should not be treated as kidney-specific without a
control of this kind.

---

## Key learning points

**What was known**

- Two coding haplotypes at APOL1, G1 and G2, carry most of the excess kidney
  disease risk in people of recent African ancestry, but only a minority of
  two-allele carriers develop disease.
- Disease depends on expression dose, so genetic modifiers of APOL1 expression
  are a plausible explanation for incomplete penetrance.
- Sequence-to-function models are increasingly proposed for prioritising
  non-coding variants, and tissue-labelled outputs are read as tissue-specific.

**This study adds**

- A tissue-labelled prediction is not a tissue-specific prediction: ranking by a
  hepatocyte output found measured kidney chromatin as well as ranking by kidney.
- The obvious readout, a gene-masked expression scorer, largely recovers
  proximity to the gene and fails a pre-specified negative control.
- Over half of an apparently novel candidate set was correlated with the
  haplotype structure that already explains the locus.

**Potential impact**

- Investigators applying these models to kidney disease loci should run a
  wrong-tissue control, a linkage control and a control locus before reporting
  candidates.
- No variant reported here is offered as prognostic or actionable; none has been
  shown to alter expression and no clinical association was tested.
- The approach narrows roughly 99,000 non-coding variants to a short list for
  functional work, which is its realistic value.

---

## 1. Introduction

Two coding haplotypes at APOL1, G1 and G2, account for much of the excess risk of
non-diabetic kidney disease in people of recent African ancestry [1]. The
association was first mapped to the adjacent gene MYH9 [2] and resolved to APOL1
two years later, with the MYH9 signal attributed to linkage disequilibrium [1].
That history is relevant to what follows.

Penetrance is incomplete: only a minority of two-allele carriers develop kidney
disease, and the mechanism is dose-dependent, so genetic modifiers of APOL1
expression are an obvious candidate explanation [3]. The locus is overwhelmingly
non-coding, with 98,875 of 107,737 catalogued variants (91.8%) outside coding
sequence, and that variation is largely unexamined. Kidney expression
quantitative trait loci are poorly powered to examine it, because kidney sample
sizes are the smallest in the major resources [4].

Sequence-to-function models predict regulatory activity directly from DNA and are
proposed as a route around that limitation [5, 6]. Earlier models of this family
perform less well in regions that distinguish one cell type from another [7, 8],
and a companion analysis found that kidney-labelled outputs of the model used
here predicted kidney expression effects no better than outputs labelled for
other tissues [9].

We therefore asked two questions together. Can such a model prioritise regulatory
candidates at this locus? And does the tissue label on its output carry the
information that its use implicitly assumes? Every positive answer was tested
against measured data the model never saw, and three controls were specified
before the results were examined.

## 2. Materials and methods

### 2.1 Locus and variants

The functional unit spanning APOL3 to MYH9, chr22:36,140,330-36,388,018 (GRCh38,
247,688 bp), fits inside a single model context window. All variants catalogued
by gnomAD v4 in that interval were retrieved with African allele frequencies. We
restricted to single-nucleotide substitutions with African allele frequency of at
least 1%, giving 1,930 variants, of which 1,865 were non-coding.

### 2.2 Scoring

Variants were scored with AlphaGenome [6] using a 1,048,576 bp context centred on
each variant, under the developers' non-commercial research terms. Two channels
were used: a gene-masked RNA scorer restricted to the APOL1 transcript, and
position-local DNase and ATAC scorers. Kidney tracks were named explicitly rather
than matched by substring, since matching on "renal" retrieves adrenal gland.

### 2.3 Proximity correction and negative control

The gene-masked scorer assigns larger predicted effects to variants nearer the
gene. We stratified by gene body, regressed log absolute effect on log distance,
and ranked on the residual. The negative control was specified in advance: the
coding risk variants should cease to be outliers once distance is controlled. If
they did not, the ranking was not to be used.

### 2.4 Independent validation

Eight ENCODE DNase-seq peak files from human kidney (GRCh38, released) were
merged over the locus. Enrichment of top-ranked variants in measured peaks was
tested by Fisher exact test against the other common variants at the same locus,
not against the genome, so that sequence composition, gene density and
mappability are matched by construction. Validation was repeated separately
against peaks from podocyte, proximal tubule, kidney epithelial, renal cortical
epithelial and kidney tubule biosamples.

### 2.5 The three controls

**Linkage disequilibrium.** Each candidate was tested against every common
variant in the APOL1 gene body, which contains both G1 and G2, using phased
1000 Genomes haplotypes in seven African populations analysed separately;
pooling inflates linkage through population structure. G2 could not be tested
directly because it is a six-base deletion absent from that call set. Intervals
were obtained by bootstrapping haplotypes within population for the strongest
partner identified in the original data.

**Wrong tissue.** The same variants were rescored and predicted accessibility
retained for kidney and for eight non-renal biosamples from a single call, so
that only the column selected differs. Each ranking was then tested for
enrichment in measured kidney chromatin.

**Control locus.** The pipeline was repeated at the beta-globin cluster,
chr11:5,150,000-5,397,688, matched on physical length, paralogue-cluster
structure and, by frequency-stratified subsampling, on variant count.

### 2.6 Code and data

All code, the frozen variant set and every result table are deposited (section
below). Thirteen numbered scripts reproduce every reported number.

## 3. Results

### 3.1 The expression channel is confounded and was discarded

Non-coding variants inside the gene body scored higher than coding ones (median
absolute effect 0.0022, n = 134, versus 0.00077, n = 1,796), so the ordering
reflects position rather than constraint. Distance explained 5.7% of the variance
in log absolute effect. After residualising, the coding risk variants remained at
the 97.6th and 96.9th percentile. The pre-specified control therefore failed, and
this channel is reported as unusable rather than omitted.

### 3.2 The accessibility channel is not confounded and its ranking is enriched

Distance explained 0.3% of the variance for podocyte DNase and none for kidney
ATAC. Measured kidney peaks covered 30,358 bp, 12.3% of the locus. Top-ranked
variants fell inside them far more often than the other common variants at the
same locus: odds ratio 8.14 for the top 25 by kidney ATAC (p = 2.4x10-6) and
6.13 for the top 200 by podocyte DNase (p = 1.5x10-22). This survived a stricter
peak definition.

### 3.3 Linkage disequilibrium removes half the candidates, and the rest cannot be called clean

Tested against 168 common gene-body variants, 7 of 15 candidates had point
estimates below r-squared 0.2, 7 were intermediate and 1 was a perfect proxy.
Bootstrapping changes how firmly this can be stated. No candidate is a proxy for
gene-body haplotype structure: 11 of 15 have an upper bound below 0.8, and that
exclusion is secure. But below 5% African allele frequency the median interval
width is 0.91, resting on 12 to 43 minor-allele copies, so those candidates
cannot be classified at all and are reported as untestable. Among the common
candidates the median interval width is 0.27, and only one, rs136204, has an
interval lying wholly below 0.2 (Table 1).

Notably, **all seven candidates are nearer a gene other than APOL1**: five lie
within MYH9 and two within or beside APOL4. This is a property of the method
rather than a surprise, since the accessibility channel is position-local and
gene-agnostic and never targeted APOL1. It does mean the candidates cannot be
described as APOL1 regulatory variants, and it returns the locus to the gene from
which the association was moved in 2010.

### 3.4 The ranking is not cell-type specific

The ranking predicts podocyte accessibility. Tested against each kidney cell type
separately, enrichment was strong in all (odds ratios 8 to 16) and podocyte led
at one of four thresholds. The original validation, against whole kidney only,
could not have detected this in either direction.

### 3.5 The tissue label carries no usable information

Ranking the same variants by predicted accessibility in eight non-renal tissues
and asking how enriched each ranking is in measured **kidney** chromatin gave
Table 2. Hepatocyte exceeded podocyte at two of four thresholds; kidney led at
one. Only the column selected from the model's output differs between rows.

### 3.6 The pipeline does not manufacture enrichment

At the beta-globin control locus the same ranking found 2 of the top 200 variants
in measured kidney chromatin against 0.7% expected (odds ratio 1.44, p = 0.44),
against 6.12 at the test locus. The control locus carries less kidney chromatin
(1.0% of its length against 12.3%), so power is lower; at the observed background
an odds ratio of 5 would have been detected at p = 0.002 and 8 at p = 5x10-6. An
effect of the size seen at APOL1-MYH9 is excluded; a small one is not.

## 4. Discussion

The headline is negative and worth stating plainly. A prediction labelled for a
tissue is not a prediction about that tissue. Asking this model about hepatocytes
identified measured kidney chromatin as well as asking it about kidney, and
asking about podocytes was no better than asking about proximal tubule. What the
accessibility channel does well is find regulatory positions, and the control
locus shows it finds them where they exist rather than everywhere. That is a
useful capability, but it is a general one, and the tissue label is close to
decorative for this task.

This is the second demonstration of that pattern in this model by our group, on a
different task with a different readout [9]. A companion analysis swapped tissue
outputs for predicting expression effects across 54 tissues; this study swaps them
for predicting chromatin accessibility at one locus. Both find that the label does
less than its use implies. It is consistent with reports on earlier models of this
family, which capture promoter determinants while largely ignoring distal
enhancers [7] and perform less well in precisely the regions that distinguish cell
types [8].

The awkward detail deserves a direct answer rather than a footnote. Five of the
seven surviving candidates lie inside MYH9, the gene to which this association was
first mapped [2] and from which it was moved in 2010 [1]. We are not reopening
that question, and nothing here is evidence that MYH9 rather than APOL1 explains
the disease association. The explanation is mundane: our ranking channel is
position-local and gene-agnostic, MYH9 occupies a large share of the interval, and
open chromatin is denser there. It nevertheless disciplines the interpretation.
Any claim that these variants regulate APOL1 would require evidence of target-gene
assignment that we do not have, and the honest description is that they are
candidate regulatory positions at the APOL1-MYH9 locus with no assigned target.

Three limitations bound the candidate list further. None of these variants has
been shown to alter expression of anything: enrichment in open chromatin is a
statement about position, not function. No clinical association was tested, so
nothing here bears on prognosis or on which carriers develop disease. And the
linkage analysis could not evaluate G2 directly, because the six-base deletion is
absent from the reference panel; we tested against all common gene-body variation
instead, which contains G2, but a reader should know the substitution was made.

We would add one methodological recommendation. Investigators applying these
models at kidney disease loci should run three controls before reporting
candidates: a wrong-tissue ranking, because tissue labels do not behave as
assumed; an explicit linkage analysis against local haplotype structure, because
it removed more than half of what initially looked novel here; and a matched
control locus, because it is the only cheap way to show the pipeline is not
producing enrichment everywhere. None requires new data.

What remains is modest and real. A method that narrows roughly 99,000 non-coding
variants at a locus of major clinical interest to a handful worth putting in front
of a bench, with honest error bars and no claim that the model knows which tissue
it is talking about.

## Data availability

Analysis code, the frozen variant set, per-variant predictions and all result
tables are available at https://github.com/PatientThread/apol1-alphagenome and
will be archived with a persistent identifier on publication. Variant frequencies
are from gnomAD v4; peak calls from the ENCODE portal; haplotypes from 1000
Genomes phase 3 (GRCh38). Model predictions were obtained under AlphaGenome's
non-commercial research terms and are redistributed subject to them; model weights
are not redistributed and were not accessed.

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

Ethical approval was not required. This study analyses publicly available,
de-identified, summary-level data and model predictions. No individual-level data
were accessed and no participants were recruited.

## References

1. Genovese G, Friedman DJ, Ross MD, et al. Association of trypanolytic ApoL1 variants with kidney disease in African Americans. Science 2010; 329: 841-845
2. Kopp JB, Smith MW, Nelson GW, et al. MYH9 is a major-effect risk gene for focal segmental glomerulosclerosis. Nat Genet 2008; 40: 1175-1184
3. Ojo AO, Adu D, Bramham K, et al. APOL1 kidney disease: a KDIGO Controversies Conference report. Kidney Int 2025; 108: 763-779
4. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. Nat Genet 2021; 53: 1290-1299
5. Linder J, Srivastava D, Yuan H, et al. Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation. Nat Genet 2025; 57: 949-961
6. Avsec Ž, Latysheva N, Cheng J, et al. Advancing regulatory variant effect prediction with AlphaGenome. Nature 2026; 649: 1206-1218
7. Karollus A, Mauermeier T, Gagneur J. Current sequence-based models capture gene expression determinants in promoters but mostly ignore distal enhancers. Genome Biol 2023; 24: 56
8. Kathail P, Shuai RW, Chung R, et al. Current genomic deep learning models display decreased performance in cell type-specific accessible regions. Genome Biol 2024; 25: 202
9. Lawrence C. Kidney representation and uncertainty in AlphaGenome regulatory variant prediction. Submitted
10. ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 2020; 583: 699-710
11. 1000 Genomes Project Consortium. A global reference for human genetic variation. Nature 2015; 526: 68-74
12. Muto Y, Wilson PC, Ledru N, et al. Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney. Nat Commun 2021; 12: 2190
13. Gisch DL, Brennan M, Lake BB, et al. The chromatin landscape of healthy and injured cell types in the human kidney. Nat Commun 2024; 15: 433
14. Sheng X, Guan Y, Ma Z, et al. Mapping the genetic architecture of human traits to cell types in the kidney identifies mechanisms of disease and potential treatments. Nat Genet 2021; 53: 1322-1333
15. Loeb GB, Kathail P, Shuai RW, et al. Variants in tubule epithelial regulatory elements mediate most heritable differences in human kidney function. Nat Genet 2024; 56: 2078-2092
16. Sasse A, Ng B, Spiro AE, et al. Benchmarking of deep neural networks for predicting personal gene expression from DNA sequence highlights shortcomings. Nat Genet 2023; 55: 2060-2064
17. Huang C, Shuai RW, Baokar P, et al. Personal transcriptome variation is poorly explained by current genomic deep learning models. Nat Genet 2023; 55: 2056-2059

## Figure legends

**Figure 1.** Predicted effect against distance to APOL1 for the gene-masked
expression scorer, with the coding risk variants marked. Their position after
residualising on distance is shown inset. The pre-specified negative control
failed, and this channel was discarded.

**Figure 2.** Enrichment of top-ranked variants in measured kidney chromatin at
the test locus, by ranking depth, for the accessibility channel. Background is
the other common variants at the same locus.

**Figure 3.** The tissue-label control. Odds ratio for enrichment in measured
kidney chromatin, by the tissue output used to rank, at four thresholds. Only the
column selected from the model's output differs between series.

**Figure 4.** The control locus. Enrichment at the beta-globin cluster against
the test locus, same pipeline and same peak source.

## Tables

**Table 1.** The seven candidates with point estimates below r-squared 0.2
against APOL1 gene-body haplotype structure, with bootstrap intervals, nearest
protein-coding gene, annotated regulatory element and the number of kidney
biosample types in which the position is measurably open.

| Variant | African AF | r² vs gene body (95% CI) | Minor allele copies | Nearest gene | Annotated element | Open in |
|---|---|---|---|---|---|---|
| rs132708 | 0.839 | 0.124 (0.02-0.29) | 1118 | APOL4, 0 bp | none | 1 of 6 |
| rs5750234 | 0.759 | 0.112 (0.02-0.26) | 1094 | APOL4, 4.3 kb | enhancer | 0 of 6 |
| rs136204 | 0.401 | 0.094 (0.03-0.17) | 443 | MYH9, 0 bp | enhancer | 6 of 6 |
| rs713797 | 0.396 | 0.120 (0.05-0.21) | 560 | MYH9, 0 bp | none | 0 of 6 |
| rs4820232 | 0.301 | 0.180 (0.07-0.32) | 318 | MYH9, 0 bp | none | 0 of 6 |
| rs6000250 | 0.111 | 0.173 (0.02-0.42) | 198 | MYH9, 0 bp | none | 5 of 6 |
| rs183925240 | 0.026 | 0.115 (0.00-0.66) | 43 | MYH9, 0 bp | enhancer | 6 of 6 |

Only rs136204 has an interval lying wholly below 0.2. rs183925240 is reported as
untestable on 43 minor-allele copies.

**Table 2.** Enrichment in measured kidney chromatin by the tissue output used to
rank the same 1,930 variants. Odds ratios against the other common variants at
the same locus.

| Ranked by | Top 25 | Top 50 | Top 100 | Top 200 |
|---|---|---|---|---|
| Podocyte (kidney) | 4.88 | 4.19 | 5.23 | 6.14 |
| Kidney (whole) | 4.06 | 3.79 | 3.04 | 2.27 |
| Hepatocyte | 4.88 | 5.07 | 6.69 | 4.81 |
| Lung | 4.06 | 4.19 | 5.23 | 5.44 |
| Frontal cortex | 4.88 | 4.19 | 4.25 | 3.85 |
| Stomach | 4.06 | 3.42 | 4.48 | 4.24 |
| Heart, left ventricle | 3.33 | 3.79 | 3.82 | 3.15 |
| Liver | 4.06 | 2.75 | 3.22 | 2.36 |
