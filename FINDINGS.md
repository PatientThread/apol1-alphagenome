# What this paper can claim

Written 18 September 2026, from the committed outputs of scripts 01 to 12. Every
number here is read from `results/`, not transcribed from memory. This exists so
the claims can be argued over before anything is drafted.

---

## The one-sentence version

A sequence-to-function model finds regulatory positions at the APOL1-MYH9 locus
competently, and finds them where they exist rather than everywhere, but it does
not use the tissue label to do it, and half the candidates it nominates turn out
to be correlated with the haplotype structure that already explains the locus.

**The primary finding is methodological. The candidate list is secondary.**

---

## What changed twice, and why the framing matters

**This is not an APOL1 paper.** All seven independent candidates are nearer a
gene other than APOL1: five sit inside **MYH9**, two in or beside APOL4. That is
defensible rather than embarrassing, because the accessibility channel used for
the ranking is position-local and gene-agnostic and never targeted APOL1. Only
the expression channel did, and that one was discarded. But it must be described
as the **APOL1-MYH9 locus** throughout, and handled carefully: MYH9 was reported
as the kidney disease risk gene in 2008 and displaced by APOL1 in 2010. Landing
back in MYH9 will draw scrutiny and should be met head on.

**The paper is now mainly a second test of the tissue-label claim**, the first
being the sister paper submitted to the European Journal of Human Genetics.

---

## The six findings, in the order they should be reported

### 1. The obvious readout is confounded and is not used

Ranking by the gene-masked expression scorer rediscovers proximity to the gene.
Non-coding variants inside the gene body outscore coding ones, so this is
position, not constraint. Distance explains 5.7% of the variance in log effect
(slope -0.104, p = 1.6e-26).

The negative control was pre-registered: if the coding risk variants remained
extreme after correction, the ranking must not be used. They did, at the 97.6th
and 96.9th percentile on the residual. **That channel is reported as unusable.**

### 2. The accessibility channel does not carry that confound

Distance explains 0.3% of the variance for podocyte DNase and 0.0% for kidney
ATAC. None of the top 15 candidates sits inside the gene body, against 6%
expected by chance.

### 3. Its top-ranked variants are enriched in independently measured kidney chromatin

Against the other common variants at the same locus, not a genome-wide
background, so the test isolates the ranking rather than rediscovering that the
locus is active.

| Ranking | Odds ratio | p |
|---|---|---|
| Top 25, kidney ATAC | 8.14 | 2.4e-06 |
| Top 200, podocyte DNase | 6.13 | 1.5e-22 |

Survives a stricter peak definition.

### 4. Linkage disequilibrium removes over half the candidates

The first analysis tested G1 and found all 15 independent. **That was falsely
reassuring**: G2 is a six-base deletion, Ensembl excludes indels from its linkage
service, and the deletion is absent from 1000 Genomes phase 3 altogether, so only
G1 was ever tested.

Testing instead against every common variant in the APOL1 gene body, which
contains both G1 and G2, using phased haplotypes across seven African populations
kept separate rather than pooled:

| Verdict | Count |
|---|---|
| Independent, r² < 0.2 | 7 |
| Partially correlated, 0.2 to 0.8 | 7 |
| Perfect proxy, r² = 1.00 | 1 |

Six of the seven independent candidates are common enough for the estimate to be
stable. **Only four of the seven sit in measured open chromatin in any kidney
cell type.** The strongest pair, rs136204 and rs183925240, are adjacent positions
inside an annotated MYH9 enhancer, open in all six sample types.

### 5. There is no cell-type specificity

The ranking predicts **podocyte** accessibility. Validated against each kidney
cell type separately, enrichment is strong everywhere (odds ratios 8 to 16) and
**podocyte leads at one of four thresholds**. The original validation used whole
kidney only and could not have detected this either way.

### 6. There is no tissue specificity either

Ranking the same variants by predicted accessibility in eight non-renal tissues,
then asking how enriched each ranking is in measured **kidney** chromatin:

| Ranked by | top 25 | top 50 | top 100 | top 200 |
|---|---|---|---|---|
| Podocyte (kidney) | 4.88 | 4.19 | 5.23 | 6.14 |
| Hepatocyte | 4.88 | 5.07 | 6.69 | 4.81 |

**Hepatocyte beats podocyte at two of four thresholds. Kidney leads at one.**
Only the tissue label changes between rows; same variants, same peaks, same
background, same test.

Internal consistency check: podocyte at top 25 gives 4.879 here against 4.877
reported by the earlier script, through a separately written code path.

### 7. The pipeline does not manufacture enrichment

Identical pipeline at the beta-globin cluster, matched on physical length,
paralogue-cluster structure and a frequency-stratified variant count:

| | Control locus | APOL1-MYH9 |
|---|---|---|
| top 100 | 1.40 | 5.22 |
| top 200 | 1.44 | 6.12 |

Two of the top 200 in measured kidney chromatin against 0.7% expected, p = 0.44.

**Honest power caveat, which must be stated.** The control locus carries far less
kidney chromatin, 1.0% of the locus against 12.3%. At the observed background a
true odds ratio of 5 would have been detected at p = 0.002 and 8 at p = 5e-6, so
an APOL1-sized effect is excluded. An odds ratio of 2 to 3 would not have been.
Do not claim a clean null.

---

## What this paper must NOT claim

- That any candidate regulates APOL1. None has been functionally validated, and
  all seven are nearer another gene.
- That the ranking is kidney-specific, podocyte-specific, or tissue-specific.
  Findings 5 and 6 exclude all three.
- That the accessibility enrichment is an APOL1 finding. It is a demonstration
  that the model detects open regulatory positions.
- That the control locus null is clean. It excludes a large effect only.

---

## Open before drafting

1. **Novelty check.** Are the six surviving candidates already described in the
   APOL1 or MYH9 regulatory literature? This is the main risk to the candidate
   list and has not been done.
2. **Rare-variant stability.** Several partially correlated candidates sit at 1
   to 3% African allele frequency, where r² between two rare variants is
   unstable. Some of those middling values may be noise.

## Target journal

Nephrology Dialysis Transplantation was proposed. Verified 14 September 2026:
free non-open-access route, no page or colour charges for submissions after
1 November 2024, impact factor 8.3, MEDLINE indexed. **Original Article ceiling
is 3,500 words including the 300-word abstract and excluding references, tables
and figures, leaving roughly 2,900 words of body text.** Whether figure legends
and the Key Learning Points box count is unverified; budget as though they do.

That is a tight envelope for a paper carrying three controls, and is an argument
for leading with the methodological finding and keeping the candidate table
short.
