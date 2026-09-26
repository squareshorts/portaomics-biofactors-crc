# Stage-1C full-universe ElasticNet development report

Protocol lock: 1.3 (2026-09-22)

External-expression loading status: **BLOCKED / NOT LOADED**

Only the six frozen development cohorts were loaded. The ElasticNet used the complete frozen 10,998-gene universe.

## Targeted convergence repair

The original GSE44861 final outer fit reached the prespecified computational ceiling of 3,000 iterations without reporting convergence. The already-selected hyperparameters were retained unchanged (l1_ratio=0.1, C=0.1), and only that final outer model was refitted with max_iter=12000. No tuning grid was rerun and no hyperparameter was changed.

The repaired GSE44861 fit converged after 3388 iterations with 734 nonzero coefficients.

## Tuned outer-fold ElasticNet models

| outer_holdout | l1_ratio | C | converged | iterations | nonzero_coefficients | calibration_intercept_train_crossfit | calibration_slope_train_crossfit | calibration_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | 0.1000 | 0.1000 | True | 2584 | 993 | 0.5372 | 0.7966 | True |
| GSE39582 | 0.9000 | 1000.0000 | True | 1817 | 10979 | 0.2764 | 0.5030 | True |
| GSE9348 | 0.1000 | 0.0100 | True | 1917 | 325 | 0.0269 | 1.0081 | True |
| GSE23878 | 0.1000 | 0.0100 | True | 2175 | 315 | 0.0075 | 1.0361 | True |
| GSE44861 | 0.1000 | 0.1000 | True | 3388 | 734 | -0.0413 | 1.1626 | True |
| GSE103512 | 0.9000 | 0.0100 | True | 1469 | 16 | 0.0193 | 1.6067 | True |

## Held-out development-cohort performance

| outer_holdout | AUROC_raw | AUROC_calibrated | AUPRC_calibrated | balanced_Brier | calibration_intercept_descriptive | calibration_slope_descriptive | calibration_fit_ok | sensitivity_at_0.5 | specificity_at_0.5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE41258 | 0.9870 | 0.9870 | 0.9937 | 0.0163 | 0.9758 | 1.5094 | True | 0.9946 | 0.9815 |
| GSE39582 | 0.9980 | 0.9972 | 1.0000 | 0.0115 | 18.9782 | 3.1227 | True | 0.9753 | 1.0000 |
| GSE9348 | 1.0000 | 1.0000 | 1.0000 | 0.0025 | 12.5796 | 9.0441 | True | 0.9857 | 1.0000 |
| GSE23878 | 1.0000 | 1.0000 | 1.0000 | 0.0078 | 11.2604 | 13.5186 | True | 1.0000 | 1.0000 |
| GSE44861 | 0.9185 | 0.9185 | 0.9093 | 0.1752 | 2.2867 | 0.5579 | True | 0.6071 | 0.9636 |
| GSE103512 | 0.9591 | 0.9591 | 0.9911 | 0.1122 | -0.6475 | 4.0053 | True | 0.9825 | 0.7500 |

## Macro-average across the six outer held-out cohorts

| AUROC_raw | AUROC_calibrated | AUPRC_calibrated | balanced_Brier | sensitivity_at_0.5 | specificity_at_0.5 |
| --- | --- | --- | --- | --- | --- |
| 0.9771 | 0.9770 | 0.9824 | 0.0543 | 0.9242 | 0.9492 |

## Convergence audit

Inner-grid fits not reporting convergence: **74**.

Final outer fits converged: **6/6**.

## Gate

`external_expression_loading_authorized` remains **false**.

`external_performance_evaluation_authorized` remains **false**.

No external cohort was loaded or evaluated in Stage 1C.