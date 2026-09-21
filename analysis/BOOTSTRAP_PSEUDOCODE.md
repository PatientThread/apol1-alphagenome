# Bootstrap procedure, as implemented

Implemented in `analysis/20_shared_draw_interaction.py`. This file is the
plain-language statement of that code, written because an earlier description
conflated two different things: which allele records exist in a replicate, and
which of those records land in the top 200 after reranking.

## Definitions

- `R` = the 39 kidney-restricted merged peaks.
- `U` = the 27 members of `R` that contain at least one eligible allele.
  **These, not all 39, are the resampling units**; the other 12 contain no
  allele and can contribute nothing.
- `F` = every eligible allele record NOT inside a kidney-restricted peak.
  This includes alleles inside *shared* kidney peaks. `F` is fixed in its
  records and in their scores. It is **not** fixed in contingency membership.
- `B` = 500 replicates, seed 20260921, numpy PCG64.

## Procedure

```
observed = deltas(all records)

draws = [ sample_with_replacement(U, size=|U|) for b in 1..B ]   # ONE list,
                                                                 # generated once

for b in 1..B:
    pool = concat( records(peak) for peak in draws[b] )   # a peak drawn k times
         + records(F)                                    # contributes k copies
    for each comparator c in the six non-podocyte rankings:
        delta_b[c] = contrast(pool, podocyte, c, KIDNEY_RESTRICTED)
                   - contrast(pool, podocyte, c, UNION)

interval[c] = percentiles(delta_b[c], 2.5, 97.5)
```

`draws` is generated **once** and reused by every comparator and by both
reference definitions, so the six estimates are paired on identical draws. An
earlier version drew fresh peaks per comparator, which left the six intervals
non-comparable with each other.

```
contrast(pool, A, B, flag):
    for col in (A, B):
        rank the WHOLE pool by |score[col]| descending
        top = first 200 records of THAT ranking
        a = |top and flag|      b = |top and not flag|
        c = |rest and flag|     d = |rest and not flag|
        logor[col] = log( (a+0.5)(d+0.5) / ((b+0.5)(c+0.5)) )
    return logor[A] - logor[B]
```

## Consequences worth being explicit about

1. **Top-set membership is an output of reranking, not an input.** A record in
   the observed top 200 need not be in a replicate's top 200, and often is not.
   Copies first enter the eligible pool; membership follows.
2. **The pool size varies between replicates**, because a drawn peak carries a
   variable number of alleles. The denominator is recomputed each time.
3. **Alternate alleles at one position always travel together**, because the
   resampling unit is keyed on the peak containing the position.
4. **The estimator is Haldane-corrected** (the +0.5 terms), so no replicate is
   infinite. The observed point estimate uses the *same* estimator. This is
   deliberately **not** the conditional maximum-likelihood odds ratio reported
   in Table 1, and the two are not interchangeable.
5. **Resampling whole peaks preserves within-peak clustering.** It does not
   account for dependence between separate peaks, and the intervals are
   conditional on `F` being fixed.

## Stability

`--stability` reruns the whole procedure at 500, 2,000 and 5,000 replicates
across three seeds: nine configurations, each producing six comparator-specific
intervals, so 54 intervals. These support the stability of zero inclusion. They
do not establish precise endpoint stability.
