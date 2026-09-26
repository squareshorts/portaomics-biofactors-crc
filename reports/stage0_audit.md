# Stage-0 data-integrity and preprocessing audit

Protocol lock: 1.0 (2026-09-22)

Overall inclusion-count status: **PASS**

Genes measurable after harmonization in all nine transcriptomic cohorts: **10,998**

No disease-association statistic, differential-expression result, AUC, classifier, PAFS score, or external predictive performance was computed in this stage.

## Sample-count audit

| accession | role | n_series_samples | tumor_all | normal_all | ambiguous | excluded | primary_tumor | primary_normal | expected_tumor_all | expected_normal_all | expected_primary_tumor | expected_primary_normal |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | development | 390 | 186 | 54 | 0 | 150 | 186 | 54 | 186 | 54 | 186 | 54 |
| GSE39582 | development | 585 | 566 | 19 | 0 | 0 | 566 | 19 | 566 | 19 | 566 | 19 |
| GSE9348 | development | 82 | 70 | 12 | 0 | 0 | 70 | 12 | 70 | 12 | 70 | 12 |
| GSE23878 | development | 59 | 35 | 24 | 0 | 0 | 35 | 24 | 35 | 24 | 35 | 24 |
| GSE44861 | development | 111 | 56 | 55 | 0 | 0 | 56 | 55 | 56 | 55 | 56 | 55 |
| GSE103512 | development | 280 | 57 | 12 | 0 | 211 | 57 | 12 | 57 | 12 | 57 | 12 |
| GSE156451 | external_locked | 144 | 72 | 72 | 0 | 0 | 72 | 72 | 72 | 72 | 72 | 72 |
| GSE106582 | external_locked | 194 | 77 | 117 | 0 | 0 | 68 | 68 | 77 | 117 | 68 | 68 |
| GSE44076 | external_locked | 246 | 98 | 148 | 0 | 0 | 98 | 98 | 98 | 148 | 98 | 98 |

## Pairing audit

| accession | complete_pairs | primary_included_samples |
| --- | --- | --- |
| GSE156451 | 72 | 144 |
| GSE106582 | 68 | 136 |
| GSE44076 | 98 | 196 |

## Expression dimensions

| accession | genes_before_intersection | primary_samples | genes_in_all_9_cohorts | ranked_missing_fraction |
| --- | --- | --- | --- | --- |
| GSE41258 | 12501 | 240 | 10998 | 0.0 |
| GSE39582 | 20843 | 585 | 10998 | 0.0 |
| GSE9348 | 20843 | 82 | 10998 | 0.0 |
| GSE23878 | 20843 | 59 | 10998 | 0.0 |
| GSE44861 | 12788 | 111 | 10998 | 0.0 |
| GSE103512 | 21113 | 69 | 10998 | 0.0 |
| GSE156451 | 26363 | 144 | 10998 | 0.0 |
| GSE106582 | 20742 | 136 | 10998 | 0.0 |
| GSE44076 | 19009 | 196 | 10998 | 0.0 |

## Probe-to-gene mapping

| accession | platform | n_input_probes | n_annotated_unambiguous_probes | n_gene_symbols |
| --- | --- | --- | --- | --- |
| GSE41258 | GPL96 | 22283.0 | 19933.0 | 12501 |
| GSE39582 | GPL570 | 54675.0 | 42901.0 | 20843 |
| GSE9348 | GPL570 | 54675.0 | 42901.0 | 20843 |
| GSE23878 | GPL570 | 54675.0 | 42901.0 | 20843 |
| GSE44861 | GPL3921 | 22277.0 | 20704.0 | 12788 |
| GSE103512 | GPL13158 | 54715.0 | 43191.0 | 21113 |
| GSE156451 | GPL24676 |  |  | 26363 |
| GSE106582 | GPL10558 | 47290.0 | 31234.0 | 20742 |
| GSE44076 | GPL13667 | 49386.0 | 47216.0 | 19009 |

## Value-range / missingness audit

| accession | min | q01 | median | q99 | max | missing_fraction |
| --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | 0.351 | 4.3 | 56.75 | 1240.0 | 527000.0 | 0.0 |
| GSE39582 | 1.345960137 | 2.50291081522 | 5.1268154505000005 | 10.974026619999998 | 16.11154945 | 0.0 |
| GSE9348 | 0.135043 | 3.5281125 | 162.004 | 11333.61875 | 144757.5 | 0.0 |
| GSE23878 | -3.4446948 | -0.0149364345599999 | 4.8012357 | 10.582977359999996 | 17.314766 | 0.0 |
| GSE44861 | 4.37639 | 5.0239201 | 6.9868 | 13.1753 | 15.3355 | 0.0 |
| GSE103512 | 2.249208561 | 2.70608675448 | 4.8985301205 | 9.250925573000009 | 13.36269314 | 0.0 |
| GSE156451 | 0.0 | 0.0 | 0.917415692366241 | 219.33008251100884 | 26507.37768675032 | 0.0 |
| GSE106582 | 6.092692032 | 6.326693187 | 6.906635473 | 12.18666669 | 14.34822719 | 0.0 |
| GSE44076 | 1.4418 | 1.94525 | 3.6743 | 10.81143700000001 | 13.57575 | 0.0 |

## Frozen-analysis gate

`classifier_analysis_authorized` remains **false**. Stage 1 must not be run until this audit has been reviewed and any technical discrepancies have been resolved without reference to downstream predictive performance.
