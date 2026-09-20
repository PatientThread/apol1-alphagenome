#!/usr/bin/env python3
"""
17_final_reconciliation.py

APOL1 paper, step 17: three sections left on stale inputs, now rebuilt.

Third editorial review found that reconciling the main comparison in step 16 had
left three other parts of the paper quoting the superseded universe. All three
are mine and all three are checkable.

  1. THE CONTROL LOCUS used ENCODE peaks queried for biosample "kidney" alone,
     while the test locus had moved to a kidney-plus-podocyte union. Comparing
     them was the same mismatch the paper is about. Its denominator was also
     1,918 against the test locus's 1,921, despite the text claiming a matched
     variant count. Rebuilt here on the same peak source.

  2. THE INTERVAL-CONTAINMENT CLAIM is stale. The control interval, 0.155 to
     6.520, contained the OLD pooled test estimate of 6.13. The reconciled
     estimate is 8.94, which the interval does NOT contain. The sentence must go
     rather than be reversed: an interval for one locus against a point estimate
     from another is not a test of the difference between them, so this is
     recomputed and reported descriptively.

  3. THE EXPRESSION SECTION still reported 134 intragenic plus 1,796 extragenic,
     totalling 1,930, the pre-deduplication universe. Recomputed on the frozen
     1,921.

  4. THE DEDUPLICATION LEDGER. Nine records were removed on rsID alone. An rsID
     can carry several alleles, so deduplication is redone on a build-specific
     chromosome/position/REF/ALT key and the removed records are listed.

Output: results/final_reconciliation.json
        results/dedup_ledger.csv
        results/control_locus_matched.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import glob
import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA, RESULTS = ROOT / "data", ROOT / "results"
CPEAKS = DATA / "control_locus_peaks_matched"
CPEAKS.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "apol1-research/1.0", "Accept": "application/json"}
CTRL_CHROM, CTRL_S, CTRL_E = "chr11", 5_150_000, 5_397_688
TEST_BP, LOCUS_BP = 44_576, 247_688
KIDNEY_TERMS = ["kidney", "glomerular visceral epithelial cell"]
PODO = "dnase_glomerular_visceral_epithelial_cell"
DEPTHS = [25, 50, 100, 200]
GENE_S, GENE_E = 36_253_010, 36_267_530
CODING = {"missense_variant", "synonymous_variant", "stop_gained",
          "frameshift_variant", "inframe_deletion", "start_lost", "stop_lost",
          "splice_acceptor_variant", "splice_donor_variant"}


def fetch(term, chrom, s, e):
    tag = term.replace(" ", "_")
    idx = CPEAKS / f"{tag}.index.json"
    if idx.exists():
        graph = json.loads(idx.read_text())
    else:
        url = ("https://www.encodeproject.org/search/?type=File&file_format=bed"
               "&output_type=peaks&assay_title=DNase-seq"
               f"&biosample_ontology.term_name={urllib.parse.quote(term)}"
               "&assembly=GRCh38&status=released"
               "&limit=4&format=json&field=accession&field=href")
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                    timeout=180) as fh:
            graph = json.load(fh)["@graph"]
        idx.write_text(json.dumps(graph)); time.sleep(0.4)
    iv, accs = [], []
    for f in graph[:4]:
        accs.append(f["accession"])
        p = CPEAKS / f"{f['accession']}.bed.gz"
        if not p.exists():
            with urllib.request.urlopen(
                    urllib.request.Request("https://www.encodeproject.org" + f["href"],
                                           headers=UA), timeout=300) as fh, open(p, "wb") as o:
                o.write(fh.read())
            time.sleep(0.4)
        with gzip.open(p, "rt") as fh:
            for line in fh:
                c = line.split("\t")
                if c and c[0] == chrom:
                    a, b = int(c[1]), int(c[2])
                    if b > s and a < e:
                        iv.append((a, b))
    return iv, accs


def merge(iv):
    if not iv:
        return []
    iv = sorted(iv); out = [list(iv[0])]
    for a, b in iv[1:]:
        if a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


def hit(p, iv):
    lo, hi = 0, len(iv) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        if p < iv[m][0]:
            hi = m - 1
        elif p >= iv[m][1]:
            lo = m + 1
        else:
            return True
    return False


def main() -> None:
    out = {}

    # ---------------------------------------------------- 4. dedup ledger
    raw = pd.concat([pd.read_csv(f, sep="\t") for f in
                     sorted(glob.glob(str(RESULTS / "scored_wrong_tissue" / "*.tsv.gz")))],
                    ignore_index=True)
    src = pd.concat([pd.read_csv(f, sep="\t") for f in
                     sorted(glob.glob(str(RESULTS / "scored" / "*.tsv.gz")))],
                    ignore_index=True)[["rsid", "pos", "ref", "alt", "consequence"]]
    raw = raw.merge(src.drop_duplicates(["pos", "ref", "alt"]),
                    on=["rsid", "pos"], how="left", suffixes=("", "_s"))
    raw["key"] = (raw["pos"].astype(str) + "_" + raw["ref"].astype(str)
                  + "_" + raw["alt"].astype(str))
    dup_key = raw[raw.duplicated("key", keep=False)].sort_values("key")
    dup_rs = raw[raw.duplicated("rsid", keep=False)].sort_values("rsid")
    ledger = dup_rs[["rsid", "pos", "ref", "alt", "key"]].copy()
    ledger["removed_by_rsid_rule"] = ledger.duplicated("rsid", keep="first")
    ledger["same_allele_key"] = ledger["key"].isin(set(dup_key["key"]))
    ledger.to_csv(RESULTS / "dedup_ledger.csv", index=False)
    by_key = raw.drop_duplicates("key")
    by_rs = raw.drop_duplicates("rsid")
    print("DEDUPLICATION")
    print(f"  raw rows                    {len(raw)}")
    print(f"  unique by rsID              {len(by_rs)}")
    print(f"  unique by pos_ref_alt key   {len(by_key)}")
    print(f"  rows sharing an rsID        {len(dup_rs)}")
    print(f"  of those, same allele too   {int(ledger['same_allele_key'].sum())}")
    if len(by_key) != len(by_rs):
        print("  NOTE: the two rules disagree; the allele key is authoritative.")
    out["dedup"] = {"raw_rows": int(len(raw)), "unique_rsid": int(len(by_rs)),
                    "unique_pos_ref_alt": int(len(by_key)),
                    "rows_sharing_rsid": int(len(dup_rs))}

    d = by_key.reset_index(drop=True)
    N = len(d)

    # ------------------------------------------- 3. expression on frozen set
    e = pd.concat([pd.read_csv(f, sep="\t") for f in
                   sorted(glob.glob(str(RESULTS / "scored" / "*.tsv.gz")))],
                  ignore_index=True)
    e["key"] = (e["pos"].astype(str) + "_" + e["ref"].astype(str)
                + "_" + e["alt"].astype(str))
    e = e.drop_duplicates("key")
    e = e.dropna(subset=["rna_cortex_of_kidney"]).copy()
    e["abs"] = e["rna_cortex_of_kidney"].abs()
    e["in_gene"] = e["pos"].between(GENE_S, GENE_E)
    e["coding"] = e["consequence"].isin(CODING)
    ig, og = e[e.in_gene], e[~e.in_gene]
    print(f"\nEXPRESSION CHANNEL, recomputed on {len(e)} deduplicated variants")
    print(f"  inside gene body   n={len(ig):>5}  median |effect| {ig['abs'].median():.5f}")
    print(f"    coding           n={int(ig.coding.sum()):>5}  "
          f"median {ig[ig.coding]['abs'].median():.5f}")
    print(f"    non-coding       n={int((~ig.coding).sum()):>5}  "
          f"median {ig[~ig.coding]['abs'].median():.5f}")
    print(f"  outside gene body  n={len(og):>5}  median |effect| {og['abs'].median():.5f}")
    print(f"  sum check          {len(ig)+len(og)} = {len(e)}  "
          f"{'OK' if len(ig)+len(og)==len(e) else 'MISMATCH'}")
    dist = np.where(e["in_gene"], 1.0,
                    np.minimum(np.abs(e["pos"] - GENE_S), np.abs(e["pos"] - GENE_E)))
    lg, ly = np.log10(np.maximum(dist, 1)), np.log10(e["abs"].values)
    ok = np.isfinite(lg) & np.isfinite(ly)
    r = np.corrcoef(lg[ok], ly[ok])[0, 1]
    print(f"  distance explains  {100*r*r:.1f}% of variance in log |effect|")
    out["expression"] = {
        "n": int(len(e)), "in_gene_n": int(len(ig)),
        "in_gene_median": float(ig["abs"].median()),
        "in_gene_coding_n": int(ig.coding.sum()),
        "in_gene_coding_median": float(ig[ig.coding]["abs"].median()),
        "in_gene_noncoding_median": float(ig[~ig.coding]["abs"].median()),
        "outside_n": int(len(og)), "outside_median": float(og["abs"].median()),
        "variance_explained_pct": round(100 * r * r, 1)}

    # ------------------------- 1 & 2. control locus on the MATCHED peak source
    iv, accs = [], []
    for t in KIDNEY_TERMS:
        a, b = fetch(t, CTRL_CHROM, CTRL_S, CTRL_E)
        iv += a; accs += b
    ctrl_peaks = merge(iv)
    bp = sum(b - a for a, b in ctrl_peaks)
    print(f"\nCONTROL LOCUS, matched peak source (kidney + podocyte)")
    print(f"  {len(ctrl_peaks)} peaks, {bp:,} bp, "
          f"{100*bp/(CTRL_E-CTRL_S):.1f}% of the control locus")
    print(f"  test locus for comparison: {TEST_BP:,} bp, "
          f"{100*TEST_BP/LOCUS_BP:.1f}% of the test locus")

    c = pd.concat([pd.read_csv(f, sep="\t") for f in
                   sorted(glob.glob(str(RESULTS / "scored_control_locus" / "*.tsv.gz")))],
                  ignore_index=True)
    c = c.drop_duplicates(["pos", "af_afr"]).dropna(subset=[PODO]).copy()
    c["k"] = c[PODO].abs()
    c = c.sort_values("k", ascending=False).reset_index(drop=True)
    c["hit"] = c["pos"].map(lambda p: hit(int(p), ctrl_peaks))
    rows = []
    print(f"\n  {len(c)} scored control variants")
    print(f"  {'depth':>6}{'a':>5}{'b':>6}{'c':>5}{'d':>7}{'OR':>8}{'p':>8}   95% CI")
    for n in DEPTHS:
        h = c["hit"].values
        a, b2 = int(h[:n].sum()), int(n - h[:n].sum())
        c2, d2 = int(h[n:].sum()), int(len(h) - n - h[n:].sum())
        t = stats.contingency.odds_ratio([[a, b2], [c2, d2]], kind="conditional")
        lo, hi = t.confidence_interval(confidence_level=0.95)
        _, p1 = stats.fisher_exact([[a, b2], [c2, d2]], alternative="greater")
        hs = "inf" if np.isinf(hi) else f"{hi:.2f}"
        print(f"  {n:>6}{a:>5}{b2:>6}{c2:>5}{d2:>7}{t.statistic:>8.2f}{p1:>8.2f}   "
              f"[{lo:.2f}, {hs}]")
        rows.append({"depth": n, "a_top_in": a, "b_top_out": b2,
                     "c_bg_in": c2, "d_bg_out": d2, "eligible_n": len(c),
                     "odds_ratio": round(float(t.statistic), 3),
                     "ci_low": round(float(lo), 3),
                     "ci_high": None if np.isinf(hi) else round(float(hi), 3),
                     "p_one_sided": round(float(p1), 4)})
    pd.DataFrame(rows).to_csv(RESULTS / "control_locus_matched.csv", index=False)
    out["control"] = {"peaks": len(ctrl_peaks), "bp": bp, "accessions": accs,
                      "eligible_n": int(len(c)), "rows": rows,
                      "test_locus_bp": TEST_BP,
                      "test_locus_pct": round(100 * TEST_BP / LOCUS_BP, 1),
                      "control_pct": round(100 * bp / (CTRL_E - CTRL_S), 1)}

    print("\n  READING: the control is reported descriptively. An interval for one")
    print("  locus against a point estimate from another is not a test of the")
    print("  difference between them, so no containment claim is made either way.")

    (RESULTS / "final_reconciliation.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote final_reconciliation.json, dedup_ledger.csv, "
          f"control_locus_matched.csv")


if __name__ == "__main__":
    main()
