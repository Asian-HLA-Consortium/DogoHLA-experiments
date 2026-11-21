# Public-only panel (2026-11-21)

The panel used in the main tests holds HPRC Asian haplotypes plus Arab and Saudi haplotypes; for one fold,
134 of its 241 added alleles occur only in the Arab and Saudi genomes, which are not public. Here the same
frozen T1K arm and frozen builder were rerun on the 66 real Asian genomes (re-pulled reads) with the panel
built by the same code from HPRC donors only (110-130 added alleles per fold, the tested donor's family left
out). Only the panel changed; a check on two donors with the original panel reproduced the original outputs
byte for byte. Descriptive: these donors had been looked at before.

| measure | T1K + IPD | T1K + full panel | T1K + public panel | builder on IPD | builder + full panel | builder + public panel |
|---|---|---|---|---|---|---|
| IPD-absent copies, whole gene rebuilt (of 134) | 0 | 31 | 24 | 18 | 45 | 40 |
| IPD-absent copies, exact final call (of 134) | 0 | 28 | 23 | 6 | 38 | 34 |
| IPD-present copies, four-field name right (of 914) | 868 | 852 | 864 | 879 | 865 | 876 |

Builder + public panel minus builder on IPD, exact final calls: +20.9 percentage points (family-clustered 95%
interval +14.0 to +27.9). Builder + public panel minus builder + full panel: -3.0 (-6.6 to 0.0). The public
panel names known alleles slightly better (876 vs 865), because the bigger panel gives wrong names to some
known alleles, mostly DPA1 and B.

Rows: `results/score_rows_pp_hprc.tsv` (arms B and D4R4n here are the public-panel arms). The full-panel and
IPD-only rows are in `../reads-repull-20261114/results/score_rows_test3.tsv` (readset T-real3).
