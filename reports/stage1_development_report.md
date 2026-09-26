# Stage-1 development-only nested LOCO report

Protocol lock: 1.1 (2026-09-22)

External-expression loading status: **BLOCKED / NOT LOADED**

This stage used only the six frozen development cohorts. GSE156451, GSE106582, and GSE44076 were not loaded by the modelling process.

## Tuned outer-fold models

| outer_holdout | method | m | C | converged | calibration_intercept_train_crossfit | calibration_slope_train_crossfit |
| --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | PAFS | 100 | 3.1623 | True | 0.0415 | 0.9460 |
| GSE41258 | M_ONLY | 5 | 10.0000 | True | -0.6370 | 1.0695 |
| GSE41258 | COHORT_ADJ | 5 | 10.0000 | True | 0.0375 | 1.1760 |
| GSE39582 | PAFS | 100 | 316.2278 | True | 1.9817 | 0.6315 |
| GSE39582 | M_ONLY | 20 | 1000.0000 | True | -0.4702 | 0.5263 |
| GSE39582 | COHORT_ADJ | 5 | 0.0316 | True | -0.1926 | 52.7269 |
| GSE9348 | PAFS | 100 | 10000.0000 | True | -0.2376 | 0.2119 |
| GSE9348 | M_ONLY | 20 | 1000.0000 | True | -0.6670 | 0.4298 |
| GSE9348 | COHORT_ADJ | 100 | 31.6228 | True | -0.4733 | 0.6535 |
| GSE23878 | PAFS | 100 | 0.0001 | True | 1.8632 | 1087.1268 |
| GSE23878 | M_ONLY | 100 | 0.0316 | True | -0.3191 | 4.1761 |
| GSE23878 | COHORT_ADJ | 100 | 10.0000 | True | -0.2794 | 0.7247 |
| GSE44861 | PAFS | 100 | 31.6228 | True | -1.4242 | 0.4923 |
| GSE44861 | M_ONLY | 100 | 3.1623 | True | -2.0373 | 1.2781 |
| GSE44861 | COHORT_ADJ | 20 | 31.6228 | True | -1.2090 | 1.4393 |
| GSE103512 | PAFS | 20 | 316.2278 | True | -0.0632 | 1.0068 |
| GSE103512 | M_ONLY | 20 | 10000.0000 | True | -0.9302 | 0.6410 |
| GSE103512 | COHORT_ADJ | 20 | 1000.0000 | True | -1.4273 | 0.8865 |

## Held-out development-cohort performance

| outer_holdout | method | AUROC_calibrated | AUPRC_calibrated | balanced_Brier | calibration_slope_descriptive | sensitivity_at_0.5 | specificity_at_0.5 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | PAFS | 0.9941 | 0.9980 | 0.0167 | 1.6422 | 0.9892 | 0.9815 |
| GSE41258 | M_ONLY | 0.9894 | 0.9971 | 0.0836 | 2.3381 | 0.7957 | 0.9815 |
| GSE41258 | COHORT_ADJ | 0.9921 | 0.9975 | 0.0213 | 1.4640 | 0.9785 | 0.9815 |
| GSE39582 | PAFS | 0.9993 | 1.0000 | 0.0674 | 3.0694 | 0.9982 | 0.8947 |
| GSE39582 | M_ONLY | 0.9988 | 1.0000 | 0.0073 | 2.8466 | 0.9894 | 1.0000 |
| GSE39582 | COHORT_ADJ | 0.9984 | 0.9999 | 0.0149 | 3.3078 | 0.9841 | 1.0000 |
| GSE9348 | PAFS | 1.0000 | 1.0000 | 0.0580 | 675.9684 | 1.0000 | 0.8333 |
| GSE9348 | M_ONLY | 1.0000 | 1.0000 | 0.0105 | 12.7851 | 1.0000 | 1.0000 |
| GSE9348 | COHORT_ADJ | 1.0000 | 1.0000 | 0.0020 | 7.7921 | 1.0000 | 1.0000 |
| GSE23878 | PAFS | 1.0000 | 1.0000 | 0.0227 | 34.1250 | 0.9714 | 1.0000 |
| GSE23878 | M_ONLY | 0.9988 | 0.9992 | 0.0240 | 3.2394 | 0.9714 | 1.0000 |
| GSE23878 | COHORT_ADJ | 0.9976 | 0.9985 | 0.0208 | 2.8303 | 0.9714 | 1.0000 |
| GSE44861 | PAFS | 0.8747 | 0.8810 | 0.2710 | 1.0703 | 0.2679 | 0.9818 |
| GSE44861 | M_ONLY | 0.8841 | 0.8925 | 0.2564 | 0.6034 | 0.3929 | 0.9818 |
| GSE44861 | COHORT_ADJ | 0.8653 | 0.8759 | 0.1967 | 0.3581 | 0.6071 | 0.9455 |
| GSE103512 | PAFS | 0.8041 | 0.9449 | 0.4542 | 0.7643 | 1.0000 | 0.0000 |
| GSE103512 | M_ONLY | 0.8363 | 0.9568 | 0.1895 | 0.7880 | 0.9474 | 0.5833 |
| GSE103512 | COHORT_ADJ | 0.8699 | 0.9659 | 0.3024 | 1.1174 | 1.0000 | 0.2500 |

## Macro-average across the six outer held-out cohorts

| method | AUROC_raw | AUROC_calibrated | AUPRC_calibrated | balanced_Brier | sensitivity_at_0.5 | specificity_at_0.5 |
| --- | --- | --- | --- | --- | --- | --- |
| PAFS | 0.9454 | 0.9454 | 0.9706 | 0.1483 | 0.8711 | 0.7819 |
| M_ONLY | 0.9512 | 0.9512 | 0.9743 | 0.0952 | 0.8495 | 0.9244 |
| COHORT_ADJ | 0.9539 | 0.9539 | 0.9730 | 0.0930 | 0.9235 | 0.8628 |

## Feature-ranking stability across outer training partitions

| method | m | mean_pairwise_Jaccard | median_pairwise_Jaccard | mean_pairwise_Kuncheva | median_pairwise_Kuncheva | n_outer_fold_pairs |
| --- | --- | --- | --- | --- | --- | --- |
| PAFS | 5 | 0.0796 | 0.0000 | 0.1329 | -0.0005 | 15 |
| PAFS | 10 | 0.1837 | 0.1765 | 0.2994 | 0.2994 | 15 |
| PAFS | 20 | 0.1629 | 0.1429 | 0.2687 | 0.2486 | 15 |
| PAFS | 50 | 0.2852 | 0.2500 | 0.4334 | 0.3973 | 15 |
| PAFS | 100 | 0.3280 | 0.3333 | 0.4840 | 0.4954 | 15 |
| M_ONLY | 5 | 0.2561 | 0.2500 | 0.3597 | 0.3997 | 15 |
| M_ONLY | 10 | 0.3417 | 0.2500 | 0.4929 | 0.3995 | 15 |
| M_ONLY | 20 | 0.4180 | 0.4286 | 0.5759 | 0.5993 | 15 |
| M_ONLY | 50 | 0.4642 | 0.4493 | 0.6263 | 0.6183 | 15 |
| M_ONLY | 100 | 0.4597 | 0.4184 | 0.6192 | 0.5862 | 15 |
| COHORT_ADJ | 5 | 0.4603 | 0.4286 | 0.6132 | 0.5998 | 15 |
| COHORT_ADJ | 10 | 0.4245 | 0.4286 | 0.5863 | 0.5996 | 15 |
| COHORT_ADJ | 20 | 0.4557 | 0.4286 | 0.6193 | 0.5993 | 15 |
| COHORT_ADJ | 50 | 0.5642 | 0.5385 | 0.7147 | 0.6986 | 15 |
| COHORT_ADJ | 100 | 0.5541 | 0.5748 | 0.7060 | 0.7275 | 15 |

## Primary comparison: PAFS minus disease-magnitude-only AUROC

| outer_holdout | delta_AUROC_PAFS_minus_M_ONLY | bootstrap_ci95_low | bootstrap_ci95_high | bootstrap_SE |
| --- | --- | --- | --- | --- |
| GSE41258 | 0.0047 | -0.0022 | 0.0161 | 0.0053 |
| GSE39582 | 0.0005 | -0.0016 | 0.0034 | 0.0013 |
| GSE9348 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| GSE23878 | 0.0012 | 0.0000 | 0.0063 | 0.0019 |
| GSE44861 | -0.0094 | -0.0292 | 0.0099 | 0.0098 |
| GSE103512 | -0.0322 | -0.1645 | 0.0857 | 0.0603 |

Random-effects pooling (DerSimonian-Laird; descriptive with six cohorts):

| pooled_delta_AUROC_DL | ci95_low | ci95_high | tau2 | I2_percent | Q | k |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0000 | -0.0000 | 0.0000 | 0.0000 | 0.0000 | 2.5087 | 6 |

## Convergence audit

| method | sum | count | all_outer_final_fits_converged |
| --- | --- | --- | --- |
| COHORT_ADJ | 6 | 6 | True |
| M_ONLY | 6 | 6 | True |
| PAFS | 6 | 6 | True |

Inner-grid optimizer fits not reporting convergence: **72**.

## Gate

`external_expression_loading_authorized` remains **false**.

`external_performance_evaluation_authorized` remains **false**.

Do not run locked external validation until this development-only report and the selected pipeline are reviewed and frozen.

Runtime: 1.4 minutes.