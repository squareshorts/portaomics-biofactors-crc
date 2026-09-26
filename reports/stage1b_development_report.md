# Stage-1B adaptive portability development report

Protocol lock: 1.2 (2026-09-22)

External-expression loading status: **BLOCKED / NOT LOADED**

Adaptive lambda grid: [0, 0.25, 0.5, 1, 2, 3]; calibration ridge alpha: 0.1.

This stage used only the six frozen development cohorts. GSE156451, GSE106582, and GSE44076 were not loaded. PAFS denotes adaptive PAFS with lambda selected strictly inside inner LOCO.

## Tuned outer-fold models

| outer_holdout | method | lambda_portability | m | C | converged | calibration_intercept_train_crossfit | calibration_slope_train_crossfit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | PAFS | 0.0000 | 5 | 10.0000 | True | -0.3159 | 1.0192 |
| GSE41258 | M_ONLY |  | 5 | 10.0000 | True | -0.3159 | 1.0192 |
| GSE41258 | COHORT_ADJ |  | 5 | 10.0000 | True | 0.0132 | 1.1171 |
| GSE39582 | PAFS | 1.0000 | 500 | 10000.0000 | True | -0.0027 | 0.3222 |
| GSE39582 | M_ONLY |  | 20 | 1000.0000 | True | -0.2097 | 0.5507 |
| GSE39582 | COHORT_ADJ |  | 5 | 0.0316 | True | -0.0029 | 1.2041 |
| GSE9348 | PAFS | 0.0000 | 20 | 1000.0000 | True | -0.3200 | 0.4290 |
| GSE9348 | M_ONLY |  | 20 | 1000.0000 | True | -0.3200 | 0.4290 |
| GSE9348 | COHORT_ADJ |  | 200 | 100.0000 | True | -0.1760 | 0.5986 |
| GSE23878 | PAFS | 0.0000 | 100 | 0.0316 | True | -0.0847 | 1.8586 |
| GSE23878 | M_ONLY |  | 100 | 0.0316 | True | -0.0847 | 1.8586 |
| GSE23878 | COHORT_ADJ |  | 500 | 3.1623 | True | -0.2373 | 0.7058 |
| GSE44861 | PAFS | 2.0000 | 200 | 31.6228 | True | -0.5965 | 0.4155 |
| GSE44861 | M_ONLY |  | 200 | 1.0000 | True | -0.4653 | 0.9765 |
| GSE44861 | COHORT_ADJ |  | 20 | 31.6228 | True | -0.2097 | 1.1465 |
| GSE103512 | PAFS | 0.0000 | 20 | 10000.0000 | True | -0.2634 | 0.7027 |
| GSE103512 | M_ONLY |  | 20 | 10000.0000 | True | -0.2634 | 0.7027 |
| GSE103512 | COHORT_ADJ |  | 20 | 1000.0000 | True | -0.3707 | 0.8470 |

## Held-out development-cohort performance

| outer_holdout | method | AUROC_calibrated | AUPRC_calibrated | balanced_Brier | calibration_slope_descriptive | sensitivity_at_0.5 | specificity_at_0.5 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | PAFS | 0.9894 | 0.9971 | 0.0656 | 2.4536 | 0.9086 | 0.9815 |
| GSE41258 | M_ONLY | 0.9894 | 0.9971 | 0.0656 | 2.4536 | 0.9086 | 0.9815 |
| GSE41258 | COHORT_ADJ | 0.9921 | 0.9975 | 0.0221 | 1.5412 | 0.9785 | 0.9815 |
| GSE39582 | PAFS | 0.9978 | 0.9999 | 0.0108 | 1.9173 | 0.9929 | 1.0000 |
| GSE39582 | M_ONLY | 0.9988 | 1.0000 | 0.0077 | 2.7204 | 0.9912 | 1.0000 |
| GSE39582 | COHORT_ADJ | 0.9984 | 0.9999 | 0.2352 | 144.8493 | 0.9894 | 1.0000 |
| GSE9348 | PAFS | 1.0000 | 1.0000 | 0.0164 | 12.0072 | 1.0000 | 1.0000 |
| GSE9348 | M_ONLY | 1.0000 | 1.0000 | 0.0164 | 12.0072 | 1.0000 | 1.0000 |
| GSE9348 | COHORT_ADJ | 1.0000 | 1.0000 | 0.0041 | 8.4597 | 1.0000 | 1.0000 |
| GSE23878 | PAFS | 0.9988 | 0.9992 | 0.0781 | 7.2785 | 0.9714 | 1.0000 |
| GSE23878 | M_ONLY | 0.9988 | 0.9992 | 0.0781 | 7.2785 | 0.9714 | 1.0000 |
| GSE23878 | COHORT_ADJ | 0.9976 | 0.9985 | 0.0207 | 2.5498 | 0.9714 | 1.0000 |
| GSE44861 | PAFS | 0.8932 | 0.8856 | 0.1851 | 1.1488 | 0.5536 | 0.9636 |
| GSE44861 | M_ONLY | 0.8695 | 0.8805 | 0.1863 | 0.8503 | 0.6071 | 0.9455 |
| GSE44861 | COHORT_ADJ | 0.8653 | 0.8759 | 0.1779 | 0.4496 | 0.6786 | 0.9455 |
| GSE103512 | PAFS | 0.8363 | 0.9568 | 0.2197 | 0.7188 | 0.9649 | 0.5000 |
| GSE103512 | M_ONLY | 0.8363 | 0.9568 | 0.2197 | 0.7188 | 0.9649 | 0.5000 |
| GSE103512 | COHORT_ADJ | 0.8699 | 0.9659 | 0.3766 | 1.1695 | 1.0000 | 0.0000 |

## Macro-average across the six outer held-out cohorts

| method | AUROC_raw | AUROC_calibrated | AUPRC_calibrated | balanced_Brier | sensitivity_at_0.5 | specificity_at_0.5 |
| --- | --- | --- | --- | --- | --- | --- |
| PAFS | 0.9526 | 0.9526 | 0.9731 | 0.0959 | 0.8986 | 0.9075 |
| M_ONLY | 0.9488 | 0.9488 | 0.9723 | 0.0956 | 0.9072 | 0.9045 |
| COHORT_ADJ | 0.9539 | 0.9539 | 0.9730 | 0.1394 | 0.9363 | 0.8212 |

## Feature-ranking stability across outer training partitions

| method | m | mean_pairwise_Jaccard | median_pairwise_Jaccard | mean_pairwise_Kuncheva | median_pairwise_Kuncheva | n_outer_fold_pairs |
| --- | --- | --- | --- | --- | --- | --- |
| PAFS | 5 | 0.1127 | 0.0000 | 0.1730 | -0.0005 | 15 |
| PAFS | 10 | 0.2149 | 0.1765 | 0.3394 | 0.2994 | 15 |
| PAFS | 20 | 0.2241 | 0.1429 | 0.3455 | 0.2486 | 15 |
| PAFS | 50 | 0.3035 | 0.2658 | 0.4468 | 0.4174 | 15 |
| PAFS | 100 | 0.3089 | 0.2903 | 0.4557 | 0.4450 | 15 |
| PAFS | 200 | 0.3664 | 0.3115 | 0.5176 | 0.4653 | 15 |
| PAFS | 500 | 0.4448 | 0.4144 | 0.5888 | 0.5663 | 15 |
| M_ONLY | 5 | 0.2561 | 0.2500 | 0.3597 | 0.3997 | 15 |
| M_ONLY | 10 | 0.3417 | 0.2500 | 0.4929 | 0.3995 | 15 |
| M_ONLY | 20 | 0.4180 | 0.4286 | 0.5759 | 0.5993 | 15 |
| M_ONLY | 50 | 0.4642 | 0.4493 | 0.6263 | 0.6183 | 15 |
| M_ONLY | 100 | 0.4597 | 0.4184 | 0.6192 | 0.5862 | 15 |
| M_ONLY | 200 | 0.5262 | 0.4760 | 0.6788 | 0.6384 | 15 |
| M_ONLY | 500 | 0.6000 | 0.5723 | 0.7345 | 0.7150 | 15 |
| COHORT_ADJ | 5 | 0.4603 | 0.4286 | 0.6132 | 0.5998 | 15 |
| COHORT_ADJ | 10 | 0.4245 | 0.4286 | 0.5863 | 0.5996 | 15 |
| COHORT_ADJ | 20 | 0.4557 | 0.4286 | 0.6193 | 0.5993 | 15 |
| COHORT_ADJ | 50 | 0.5642 | 0.5385 | 0.7147 | 0.6986 | 15 |
| COHORT_ADJ | 100 | 0.5541 | 0.5748 | 0.7060 | 0.7275 | 15 |
| COHORT_ADJ | 200 | 0.6139 | 0.6667 | 0.7515 | 0.7963 | 15 |
| COHORT_ADJ | 500 | 0.6939 | 0.7094 | 0.8092 | 0.8219 | 15 |

## Primary development comparison: adaptive PAFS minus disease-magnitude-only AUROC

| outer_holdout | delta_AUROC_PAFS_minus_M_ONLY | bootstrap_ci95_low | bootstrap_ci95_high | bootstrap_SE |
| --- | --- | --- | --- | --- |
| GSE41258 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| GSE39582 | -0.0010 | -0.0035 | 0.0004 | 0.0010 |
| GSE9348 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| GSE23878 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| GSE44861 | 0.0237 | -0.0152 | 0.0588 | 0.0188 |
| GSE103512 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

Macro mean delta AUROC (adaptive PAFS - M_ONLY): 0.0038.

Legacy inverse-variance pooling (descriptive only; not used for method selection because a near-zero cohort SE can dominate):

| pooled_delta_AUROC_DL | ci95_low | ci95_high | tau2 | I2_percent | Q | k |
| --- | --- | --- | --- | --- | --- | --- |
| -0.0000 | -0.0000 | 0.0000 | 0.0000 | 0.0000 | 2.5718 | 6 |

## Convergence audit

| method | sum | count | all_outer_final_fits_converged |
| --- | --- | --- | --- |
| COHORT_ADJ | 6 | 6 | True |
| M_ONLY | 6 | 6 | True |
| PAFS | 6 | 6 | True |

Inner-grid optimizer fits not reporting convergence: **0**.

## Gate

`external_expression_loading_authorized` remains **false**.

`external_performance_evaluation_authorized` remains **false**.

Do not run locked external validation until this development-only report and the selected pipeline are reviewed and frozen.

Runtime: 4.6 minutes.