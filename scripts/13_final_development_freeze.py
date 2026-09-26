from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"

LOCK_PATH = CONFIG / "final_development_lock.json"
OUT = RESULTS / "final_development"

RIDGE_SCRIPT = ROOT / "scripts" / "11_stage1b_adaptive_loco.py"
ELASTIC_SCRIPT = ROOT / "scripts" / "12_stage1c_elasticnet.py"

FINAL_EXTERNAL_LOCK = CONFIG / "final_external_lock.json"
REPORT_PATH = REPORTS / "final_development_freeze.md"


def load_module(name, path):
    if not path.exists():
        raise SystemExit(f"Missing required script: {path}")

    spec = importlib.util.spec_from_file_location(name, path)

    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")

    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)

    return mod


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()

    with path.open("rb") as f:
        while True:
            b = f.read(chunk)

            if not b:
                break

            h.update(b)

    return h.hexdigest()


def fit_final_ridge(
    r,
    data,
    genes,
    dev,
    method,
    m,
    C,
    score_cache,
    alpha,
    lam=None
):
    r.ACTIVE_OUTER = "FINAL"

    if method == "PAFS":
        r.SELECTED_LAMBDA["FINAL"] = float(lam)

    cf_p, cf_y, cf_w = r.crossfit_for_calibration(
        data,
        genes,
        dev,
        method,
        int(m),
        float(C),
        score_cache=score_cache
    )

    r.CALIBRATION_RIDGE_ALPHA = float(alpha)

    cal_a, cal_b, cal_ok = r.calibration_fit_ridge(
        cf_p,
        cf_y,
        cf_w,
        nonnegative_slope=True
    )

    if not cal_ok:
        raise RuntimeError(
            f"{method}: calibration did not converge"
        )

    scores = r.feature_scores(
        data,
        dev,
        genes
    )

    rank = r.ranking_indices(
        scores,
        method
    )

    idx = rank[:int(m)]

    X, y, _, _, _, w = r.pooled_data(
        data,
        dev,
        idx
    )

    mu = X.mean(
        axis=0,
        dtype=np.float64
    )

    sd = X.std(
        axis=0,
        dtype=np.float64
    )

    sd = np.where(
        sd < 1e-8,
        1.0,
        sd
    )

    Xz = (
        X - mu
    ) / sd

    theta, ok, nit, obj = r.ridge_logistic_fit(
        Xz,
        y,
        w,
        float(C)
    )

    if not ok:
        raise RuntimeError(
            f"{method}: final all-six fit did not converge"
        )

    selected = scores.iloc[
        idx
    ].copy()

    selected.insert(
        0,
        "selected_rank",
        np.arange(
            1,
            len(selected) + 1
        )
    )

    selected.insert(
        1,
        "method",
        method
    )

    selected.insert(
        2,
        "lambda_portability",
        ""
        if lam is None
        else float(lam)
    )

    selected.insert(
        3,
        "m",
        int(m)
    )

    selected.insert(
        4,
        "C",
        float(C)
    )

    selected[
        "training_mean"
    ] = mu

    selected[
        "training_sd"
    ] = sd

    selected[
        "coefficient_standardized"
    ] = theta[1:]

    info = {
        "method":
            method,

        "lambda_portability":
            None
            if lam is None
            else float(lam),

        "m":
            int(m),

        "C":
            float(C),

        "intercept_standardized":
            float(theta[0]),

        "fit_converged":
            bool(ok),

        "fit_iterations":
            int(nit),

        "objective":
            float(obj),

        "calibration_intercept":
            float(cal_a),

        "calibration_slope":
            float(cal_b),

        "calibration_converged":
            bool(cal_ok)
    }

    return selected, info


def fit_final_elastic(
    e,
    data,
    genes,
    dev,
    l1,
    C,
    pred_cache,
    alpha
):
    pp = []
    yy = []
    cc = []

    for h in dev:

        pp.append(
            pred_cache[
                (
                    h,
                    float(l1),
                    float(C)
                )
            ]
        )

        yy.append(
            data[h].y
        )

        cc.extend(
            [h]
            * len(
                data[h].y
            )
        )

    pp = np.concatenate(pp)
    yy = np.concatenate(yy)

    cc = np.asarray(
        cc,
        dtype=object
    )

    cw = e.b.cohort_class_weights(
        yy,
        cc
    )

    cal_a, cal_b, cal_ok = e.calfit(
        pp,
        yy,
        cw,
        alpha
    )

    if not cal_ok:
        raise RuntimeError(
            "ElasticNet calibration did not converge"
        )

    X, y, _, _, _, w = e.b.pooled_data(
        data,
        dev
    )

    mu = X.mean(
        axis=0,
        dtype=np.float64
    )

    sd = X.std(
        axis=0,
        dtype=np.float64
    )

    sd = np.where(
        sd < 1e-8,
        1.0,
        sd
    )

    Xz = (
        X - mu
    ) / sd

    model, ok, nit, nz = e.fit_model(
        Xz,
        y,
        w,
        float(l1),
        float(C),
        20260922,
        12000,
        False,
        None
    )

    if not ok:
        raise RuntimeError(
            "ElasticNet final all-six fit did not converge"
        )

    model_df = pd.DataFrame(
        {
            "gene":
                genes.astype(str),

            "training_mean":
                mu,

            "training_sd":
                sd,

            "coefficient_standardized":
                model.coef_[0]
        }
    )

    info = {
        "method":
            "ELASTICNET",

        "l1_ratio":
            float(l1),

        "C":
            float(C),

        "intercept_standardized":
            float(
                model.intercept_[0]
            ),

        "fit_converged":
            bool(ok),

        "fit_iterations":
            int(nit),

        "nonzero_coefficients":
            int(nz),

        "calibration_intercept":
            float(cal_a),

        "calibration_slope":
            float(cal_b),

        "calibration_converged":
            bool(cal_ok)
    }

    return model_df, info


def main():

    t0 = time.time()

    lock = json.loads(
        LOCK_PATH.read_text(
            encoding="utf-8-sig"
        )
    )

    if (
        lock.get(
            "stage0_review_status"
        )
        != "PASS"
    ):
        raise SystemExit(
            "Final development lock does not record Stage-0 PASS"
        )

    if (
        lock.get(
            "external_expression_loading_authorized"
        )
        is not False
    ):
        raise SystemExit(
            "External loading must remain false"
        )

    if (
        lock.get(
            "external_performance_evaluation_authorized"
        )
        is not False
    ):
        raise SystemExit(
            "External evaluation must remain false"
        )

    r = load_module(
        "freeze_ridge",
        RIDGE_SCRIPT
    )

    e = load_module(
        "freeze_elastic",
        ELASTIC_SCRIPT
    )

    r.verify_stage0()

    r.CALIBRATION_RIDGE_ALPHA = float(
        lock.get(
            "calibration_ridge_alpha",
            0.1
        )
    )

    r_lock = r.load_lock()

    if (
        [
            float(x)
            for x in r_lock[
                "adaptive_lambda_grid"
            ]
        ]
        !=
        [
            float(x)
            for x in lock[
                "adaptive_lambda_grid"
            ]
        ]
    ):
        raise RuntimeError(
            "Stage-1B lambda grid differs from final lock"
        )

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    genes, data, _ = r.load_development_data(
        lock
    )

    dev = list(
        lock[
            "development_accessions"
        ]
    )

    external = set(
        lock[
            "external_accessions_locked"
        ]
    )

    if (
        set(dev)
        & external
    ):
        raise RuntimeError(
            "Development/external accession overlap"
        )

    m_grid = [
        int(x)
        for x in lock[
            "signature_sizes"
        ]
    ]

    ridge_C = [
        10.0 ** float(x)
        for x in lock[
            "ridge_C_grid_log10"
        ]
    ]

    l1_grid = [
        float(x)
        for x in lock[
            "elasticnet_l1_ratio_grid"
        ]
    ]

    elastic_C = [
        10.0 ** float(x)
        for x in lock[
            "elasticnet_C_grid_log10"
        ]
    ]

    alpha = float(
        lock.get(
            "calibration_ridge_alpha",
            0.1
        )
    )

    print(
        "\n=== FINAL DEVELOPMENT FREEZE ==="
    )

    print(
        "External cohorts remain BLOCKED / NOT LOADED."
    )

    print(
        "\n[1/2] Final six-cohort LOCO tuning: "
        "adaptive PAFS, M_ONLY, COHORT_ADJ"
    )

    _, ridge_tuning, ridge_diag, score_cache = (
        r.tune_one_outer(
            data,
            genes,
            dev,
            m_grid,
            ridge_C,
            "FINAL"
        )
    )

    # Require convergence in all six folds.
    rt = ridge_tuning.copy()
    rd = ridge_diag.copy()

    rt["lambda_key"] = pd.to_numeric(
        rt[
            "lambda_portability"
        ],
        errors="coerce"
    ).fillna(-1.0)

    rd["lambda_key"] = pd.to_numeric(
        rd[
            "lambda_portability"
        ],
        errors="coerce"
    ).fillna(-1.0)

    conv = (
        rd.groupby(
            [
                "method",
                "lambda_key",
                "m",
                "C"
            ],
            dropna=False
        )[
            "converged"
        ]
        .all()
        .reset_index(
            name="all_six_folds_converged"
        )
    )

    rt = rt.merge(
        conv,
        on=[
            "method",
            "lambda_key",
            "m",
            "C"
        ],
        how="left"
    )

    rt[
        "all_six_folds_converged"
    ] = (
        rt[
            "all_six_folds_converged"
        ]
        .fillna(False)
        .astype(bool)
    )

    ridge_best = {}
    selected_tuning_rows = {}

    for method in [
        "PAFS",
        "M_ONLY",
        "COHORT_ADJ"
    ]:

        cand = rt[
            (
                rt[
                    "method"
                ]
                == method
            )
            &
            (
                rt[
                    "all_six_folds_converged"
                ]
            )
        ].copy()

        if cand.empty:
            raise RuntimeError(
                f"No fully converged final candidate for {method}"
            )

        cand = cand.sort_values(
            [
                "macro_inner_auc",
                "m",
                "C",
                "lambda_key"
            ],
            ascending=[
                False,
                True,
                True,
                True
            ],
            kind="mergesort"
        )

        row = cand.iloc[0]

        ridge_best[
            method
        ] = (
            int(
                row["m"]
            ),
            float(
                row["C"]
            )
        )

        selected_tuning_rows[
            method
        ] = row

    lam = float(
        selected_tuning_rows[
            "PAFS"
        ][
            "lambda_key"
        ]
    )

    r.SELECTED_LAMBDA[
        "FINAL"
    ] = lam

    r.ACTIVE_OUTER = "FINAL"

    rt.drop(
        columns=[
            "lambda_key"
        ]
    ).to_csv(
        OUT
        / "ridge_final_loco_tuning.csv.gz",
        index=False,
        compression="gzip"
    )

    ridge_diag.to_csv(
        OUT
        / "ridge_final_loco_optimizer_diagnostics.csv.gz",
        index=False,
        compression="gzip"
    )

    ridge_rows = []
    final_rows = []
    artifact_paths = {}

    for method in [
        "PAFS",
        "M_ONLY",
        "COHORT_ADJ"
    ]:

        m, C = ridge_best[
            method
        ]

        tr = selected_tuning_rows[
            method
        ]

        ridge_rows.append(
            {
                "method":
                    method,

                "lambda_portability":
                    lam
                    if method == "PAFS"
                    else np.nan,

                "m":
                    int(m),

                "C":
                    float(C),

                "macro_loco_auc":
                    float(
                        tr[
                            "macro_inner_auc"
                        ]
                    ),

                "min_loco_auc":
                    float(
                        tr[
                            "inner_auc_min"
                        ]
                    ),

                "max_loco_auc":
                    float(
                        tr[
                            "inner_auc_max"
                        ]
                    )
            }
        )

        print(
            f"  {method}: "
            f"lambda="
            f"{lam if method == 'PAFS' else '-'}, "
            f"m={m}, "
            f"C={C:.6g}, "
            f"macro LOCO AUROC="
            f"{float(tr['macro_inner_auc']):.4f}"
        )

        selected, info = fit_final_ridge(
            r,
            data,
            genes,
            dev,
            method,
            m,
            C,
            score_cache,
            alpha,
            lam=(
                lam
                if method == "PAFS"
                else None
            )
        )

        info[
            "selected_macro_loco_auc"
        ] = float(
            tr[
                "macro_inner_auc"
            ]
        )

        info[
            "selected_min_loco_auc"
        ] = float(
            tr[
                "inner_auc_min"
            ]
        )

        info[
            "selected_max_loco_auc"
        ] = float(
            tr[
                "inner_auc_max"
            ]
        )

        csv_path = (
            OUT
            / f"{method.lower()}_final_model.csv.gz"
        )

        json_path = (
            OUT
            / f"{method.lower()}_final_model.json"
        )

        selected.to_csv(
            csv_path,
            index=False,
            compression="gzip"
        )

        json_path.write_text(
            json.dumps(
                info,
                indent=2
            ),
            encoding="utf-8"
        )

        artifact_paths[
            f"{method.lower()}_model_csv"
        ] = csv_path

        artifact_paths[
            f"{method.lower()}_model_json"
        ] = json_path

        final_rows.append(
            info
        )

    print(
        "\n[2/2] Final six-cohort LOCO tuning: "
        "full-universe ElasticNet"
    )

    (
        elastic_best,
        elastic_tuning,
        elastic_diag,
        pred_cache
    ) = e.tune_outer(
        data,
        dev,
        l1_grid,
        elastic_C,
        "FINAL"
    )

    elastic_tuning.to_csv(
        OUT
        / "elasticnet_final_loco_tuning.csv.gz",
        index=False,
        compression="gzip"
    )

    elastic_diag.to_csv(
        OUT
        / "elasticnet_final_loco_optimizer_diagnostics.csv.gz",
        index=False,
        compression="gzip"
    )

    l1, C = elastic_best

    erow = elastic_tuning[
        np.isclose(
            elastic_tuning[
                "l1_ratio"
            ].astype(float),
            float(l1)
        )
        &
        np.isclose(
            elastic_tuning[
                "C"
            ].astype(float),
            float(C)
        )
    ]

    if len(erow) != 1:
        raise RuntimeError(
            "Could not recover selected ElasticNet tuning row"
        )

    erow = erow.iloc[0]

    print(
        f"  ELASTICNET: "
        f"l1_ratio={l1:g}, "
        f"C={C:g}, "
        f"macro LOCO AUROC="
        f"{float(erow['macro_inner_auc']):.4f}"
    )

    elastic_model, elastic_info = fit_final_elastic(
        e,
        data,
        genes,
        dev,
        l1,
        C,
        pred_cache,
        alpha
    )

    elastic_info[
        "selected_macro_loco_auc"
    ] = float(
        erow[
            "macro_inner_auc"
        ]
    )

    elastic_info[
        "selected_min_loco_auc"
    ] = float(
        erow[
            "inner_auc_min"
        ]
    )

    elastic_info[
        "selected_max_loco_auc"
    ] = float(
        erow[
            "inner_auc_max"
        ]
    )

    elastic_csv = (
        OUT
        / "elasticnet_final_model.csv.gz"
    )

    elastic_json = (
        OUT
        / "elasticnet_final_model.json"
    )

    elastic_model.to_csv(
        elastic_csv,
        index=False,
        compression="gzip"
    )

    elastic_json.write_text(
        json.dumps(
            elastic_info,
            indent=2
        ),
        encoding="utf-8"
    )

    artifact_paths[
        "elasticnet_model_csv"
    ] = elastic_csv

    artifact_paths[
        "elasticnet_model_json"
    ] = elastic_json

    final_rows.append(
        elastic_info
    )

    ridge_view = pd.DataFrame(
        ridge_rows
    )

    elastic_view = pd.DataFrame(
        [
            {
                "method":
                    "ELASTICNET",

                "l1_ratio":
                    float(l1),

                "C":
                    float(C),

                "macro_loco_auc":
                    float(
                        erow[
                            "macro_inner_auc"
                        ]
                    ),

                "min_loco_auc":
                    float(
                        erow[
                            "inner_auc_min"
                        ]
                    ),

                "max_loco_auc":
                    float(
                        erow[
                            "inner_auc_max"
                        ]
                    )
            }
        ]
    )

    final_summary = pd.DataFrame(
        final_rows
    )

    summary_path = (
        OUT
        / "final_model_summary.csv"
    )

    final_summary.to_csv(
        summary_path,
        index=False
    )

    artifact_paths[
        "final_model_summary"
    ] = summary_path

    ridge_nonconv = (
        int(
            (
                ~ridge_diag[
                    "converged"
                ].astype(bool)
            ).sum()
        )
        if not ridge_diag.empty
        else 0
    )

    elastic_nonconv = (
        int(
            (
                ~elastic_diag[
                    "converged"
                ].astype(bool)
            ).sum()
        )
        if not elastic_diag.empty
        else 0
    )

    artifact_hashes = {
        k: {
            "path":
                str(
                    p.relative_to(
                        ROOT
                    )
                ),

            "sha256":
                sha256(p)
        }

        for k, p
        in artifact_paths.items()
    }

    safe_models = (
        final_summary
        .astype(object)
        .where(
            pd.notna(
                final_summary
            ),
            None
        )
        .to_dict(
            orient="records"
        )
    )

    external_lock = {
        "lock_version":
            "2.0-pre-external",

        "lock_date":
            lock.get(
                "lock_date"
            ),

        "status":
            "DEVELOPMENT_FROZEN_AWAITING_EXTERNAL_AUTHORIZATION",

        "development_accessions":
            dev,

        "external_accessions":
            list(
                lock[
                    "external_accessions_locked"
                ]
            ),

        "feature_universe_n":
            int(
                lock[
                    "feature_universe_n"
                ]
            ),

        "primary_method":
            "PAFS",

        "primary_comparator":
            "M_ONLY",

        "secondary_feature_comparator":
            "COHORT_ADJ",

        "predictive_benchmark":
            "ELASTICNET",

        "primary_external_estimand":
            "Within each external cohort, delta AUROC = adaptive PAFS minus M_ONLY",

        "primary_external_uncertainty":
            "2000-replicate participant-cluster bootstrap within each external cohort",

        "primary_external_summary":
            "Report each external-cohort delta AUROC and the unweighted macro-mean delta across the three external cohorts",

        "absolute_external_performance":
            "Report AUROC separately by method and external cohort plus unweighted macro-average",

        "secondary_external_endpoints": [
            "AUPRC",
            "balanced_Brier",
            "calibration_intercept",
            "calibration_slope",
            "sensitivity_at_0.5",
            "specificity_at_0.5"
        ],

        "classification_threshold":
            0.5,

        "external_primary_samples":
            "Use Stage-0 primary_include samples only",

        "GSE44076_secondary_stress_test":
            "After primary paired validation, evaluate the prespecified healthy-control stress test separately; it must not replace the primary adjacent-normal result",

        "orthogonal_multiomics_signature_source":
            "Final all-six-development adaptive PAFS selected genes only; methylation/protein results cannot feed back into transcriptomic selection",

        "external_expression_loading_authorized":
            False,

        "external_performance_evaluation_authorized":
            False,

        "frozen_models":
            safe_models,

        "artifact_hashes":
            artifact_hashes,

        "source_hashes": {
            "final_development_lock":
                sha256(
                    LOCK_PATH
                ),

            "stage1b_working_script":
                sha256(
                    RIDGE_SCRIPT
                ),

            "stage1c_working_script":
                sha256(
                    ELASTIC_SCRIPT
                )
        },

        "notes":
            "No external cohort was loaded or evaluated during the final development freeze. External gates remain false until review."
    }

    FINAL_EXTERNAL_LOCK.write_text(
        json.dumps(
            external_lock,
            indent=2
        ),
        encoding="utf-8"
    )

    report = [
        "# Final development freeze before external validation",
        "",
        (
            f"Protocol lock: "
            f"{lock.get('lock_version', '')} "
            f"({lock.get('lock_date', '')})"
        ),
        "",
        "External-expression loading status: **BLOCKED / NOT LOADED**",
        "",
        "External-performance evaluation status: **BLOCKED / NOT RUN**",
        "",
        (
            "All final hyperparameters were selected by "
            "leave-one-cohort-out validation across the six "
            "development cohorts only. No external expression "
            "matrix or external outcome was loaded."
        ),
        "",
        "## Frozen feature-selection models",
        "",
        r.write_markdown_table(
            ridge_view
        ),
        "",
        "## Frozen full-universe ElasticNet benchmark",
        "",
        r.write_markdown_table(
            elastic_view
        ),
        "",
        "## Final all-development fitted models",
        "",
        r.write_markdown_table(
            final_summary
        ),
        "",
        "## Convergence audit",
        "",
        (
            "Ridge-grid fits not reporting convergence: "
            f"**{ridge_nonconv}**."
        ),
        "",
        (
            "ElasticNet-grid fits not reporting convergence: "
            f"**{elastic_nonconv}**."
        ),
        "",
        (
            "The selected final ElasticNet hyperparameter pair "
            "was required to converge in all six LOCO folds. "
            "Every final all-development model was also required "
            "to converge."
        ),
        "",
        "## Frozen external estimand",
        "",
        (
            "Primary comparison: adaptive PAFS minus M_ONLY "
            "AUROC within each external cohort, with "
            "2,000-replicate participant-cluster bootstrap "
            "uncertainty. Report the three cohort-specific "
            "effects and their unweighted macro mean."
        ),
        "",
        (
            "Absolute AUROC for all four methods is also reported "
            "separately by cohort and macro-averaged."
        ),
        "",
        (
            "The GSE44076 healthy-control analysis is a secondary "
            "stress test performed only after the primary paired "
            "tumor/adjacent-normal evaluation."
        ),
        "",
        (
            "Orthogonal methylation/protein validation uses the "
            "final adaptive-PAFS gene set and cannot feed back "
            "into transcriptomic model selection."
        ),
        "",
        "## External gate",
        "",
        "`external_expression_loading_authorized` remains **false**.",
        "",
        "`external_performance_evaluation_authorized` remains **false**.",
        "",
        (
            "Frozen model artifacts and SHA-256 hashes are "
            "recorded in `config/final_external_lock.json`. "
            "External validation must not be run until this "
            "report is reviewed."
        ),
        "",
        (
            f"Runtime: "
            f"{(time.time() - t0) / 60.0:.1f} minutes."
        )
    ]

    REPORT_PATH.write_text(
        "\n".join(
            report
        ),
        encoding="utf-8"
    )

    print(
        "\n=== FINAL DEVELOPMENT FREEZE COMPLETE ==="
    )

    print(
        f"Report: {REPORT_PATH}"
    )

    print(
        f"External lock: {FINAL_EXTERNAL_LOCK}"
    )

    print(
        "External expression loaded: NO"
    )

    print(
        "External performance evaluated: NO"
    )

    print(
        "External gates remain FALSE pending review."
    )


if __name__ == "__main__":
    main()
