# Final development freeze before external validation

Protocol lock: 1.4 (2026-09-23)

External-expression loading status: **BLOCKED / NOT LOADED**

External-performance evaluation status: **BLOCKED / NOT RUN**

All final hyperparameters were selected by leave-one-cohort-out validation across the six development cohorts only. No external expression matrix or external outcome was loaded.

## Frozen feature-selection models

| method | lambda_portability | m | C | macro_loco_auc | min_loco_auc | max_loco_auc |
| --- | --- | --- | --- | --- | --- | --- |
| PAFS | 0.0000 | 20 | 0.0316 | 0.9805 | 0.9029 | 1.0000 |
| M_ONLY |  | 20 | 0.0316 | 0.9805 | 0.9029 | 1.0000 |
| COHORT_ADJ |  | 100 | 100.0000 | 0.9776 | 0.9003 | 1.0000 |

## Frozen full-universe ElasticNet benchmark

| method | l1_ratio | C | macro_loco_auc | min_loco_auc | max_loco_auc |
| --- | --- | --- | --- | --- | --- |
| ELASTICNET | 0.1000 | 0.0100 | 0.9795 | 0.9065 | 1.0000 |

## Final all-development fitted models

| method | lambda_portability | m | C | intercept_standardized | fit_converged | fit_iterations | objective | calibration_intercept | calibration_slope | calibration_converged | selected_macro_loco_auc | selected_min_loco_auc | selected_max_loco_auc | l1_ratio | nonzero_coefficients |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PAFS | 0.0000 | 20.0000 | 0.0316 | 0.1360 | True | 12 | 0.6515 | -0.0104 | 1.5986 | True | 0.9805 | 0.9029 | 1.0000 |  |  |
| M_ONLY |  | 20.0000 | 0.0316 | 0.1360 | True | 12 | 0.6515 | -0.0104 | 1.5986 | True | 0.9805 | 0.9029 | 1.0000 |  |  |
| COHORT_ADJ |  | 100.0000 | 100.0000 | 4.6235 | True | 46 | 0.0770 | -0.3213 | 0.5577 | True | 0.9776 | 0.9003 | 1.0000 |  |  |
| ELASTICNET |  |  | 0.0100 | 2.3592 | True | 2056 |  | -0.0281 | 1.0002 | True | 0.9795 | 0.9065 | 1.0000 | 0.1000 | 313.0000 |

## Convergence audit

Ridge-grid fits not reporting convergence: **0**.

ElasticNet-grid fits not reporting convergence: **13**.

The selected final ElasticNet hyperparameter pair was required to converge in all six LOCO folds. Every final all-development model was also required to converge.

## Frozen external estimand

Primary comparison: adaptive PAFS minus M_ONLY AUROC within each external cohort, with 2,000-replicate participant-cluster bootstrap uncertainty. Report the three cohort-specific effects and their unweighted macro mean.

Absolute AUROC for all four methods is also reported separately by cohort and macro-averaged.

The GSE44076 healthy-control analysis is a secondary stress test performed only after the primary paired tumor/adjacent-normal evaluation.

Orthogonal methylation/protein validation uses the final adaptive-PAFS gene set and cannot feed back into transcriptomic model selection.

## External gate

`external_expression_loading_authorized` remains **false**.

`external_performance_evaluation_authorized` remains **false**.

Frozen model artifacts and SHA-256 hashes are recorded in `config/final_external_lock.json`. External validation must not be run until this report is reviewed.

Runtime: 141.4 minutes.