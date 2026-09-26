from __future__ import annotations
import importlib.util, json, sys, time, warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "config" / "stage1c_lock.json"
OUT = ROOT / "results" / "stage1c"
REPORTS = ROOT / "reports"
BASE = ROOT / "scripts" / "10_stage1_nested_loco.py"


def load_base():
    sp = importlib.util.spec_from_file_location("s1base", BASE)

    if sp is None or sp.loader is None:
        raise RuntimeError(f"Cannot load {BASE}")

    m = importlib.util.module_from_spec(sp)
    sys.modules[sp.name] = m
    sp.loader.exec_module(m)

    return m


b = load_base()


def load_lock():

    x = json.loads(
        LOCK.read_text(
            encoding="utf-8-sig"
        )
    )

    if (
        x.get("stage0_review_status") != "PASS"
        or not x.get(
            "classifier_analysis_authorized",
            False
        )
    ):
        raise SystemExit(
            "Stage-1C is not authorized"
        )

    if (
        x.get(
            "external_expression_loading_authorized",
            False
        )
        or x.get(
            "external_performance_evaluation_authorized",
            False
        )
    ):
        raise SystemExit(
            "External-data safety gate failed"
        )

    return x


def calfit(p, y, w, alpha):

    x = logit(
        np.clip(
            np.asarray(p, float),
            1e-6,
            1 - 1e-6
        )
    )

    y = np.asarray(y, float)
    w = np.asarray(w, float)

    sw = w.sum()

    def fg(t):

        eta = (
            t[0]
            + t[1] * x
        )

        pp = expit(eta)

        r = (
            w
            * (pp - y)
            / sw
        )

        loss = (
            np.sum(
                w
                * (
                    np.logaddexp(
                        0,
                        eta
                    )
                    - y * eta
                )
            )
            / sw
        )

        loss += (
            0.5
            * alpha
            * (
                t[0] ** 2
                + (t[1] - 1) ** 2
            )
        )

        grad = np.array(
            [
                r.sum()
                + alpha * t[0],

                np.sum(r * x)
                + alpha * (t[1] - 1)
            ]
        )

        return float(loss), grad

    r = minimize(
        lambda t: fg(t)[0],
        [0.0, 1.0],
        jac=lambda t: fg(t)[1],
        method="L-BFGS-B",
        bounds=[
            (None, None),
            (0, None)
        ],
        options={
            "maxiter": 200,
            "ftol": 1e-10,
            "gtol": 1e-7
        }
    )

    return (
        float(r.x[0]),
        float(r.x[1]),
        bool(r.success)
    )


def calapply(p, a, c):

    return expit(
        a
        + c
        * logit(
            np.clip(
                np.asarray(
                    p,
                    float
                ),
                1e-6,
                1 - 1e-6
            )
        )
    )


def fit_model(
    X,
    y,
    w,
    l1,
    C,
    seed,
    max_iter=2000,
    warm=False,
    model=None
):

    if model is None:

        model = LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            l1_ratio=float(l1),
            C=float(C),
            fit_intercept=True,
            max_iter=max_iter,
            tol=1e-4,
            warm_start=warm,
            random_state=seed,
            n_jobs=1
        )

    else:

        model.set_params(
            C=float(C),
            l1_ratio=float(l1),
            max_iter=max_iter
        )

    with warnings.catch_warnings(
        record=True
    ) as ww:

        warnings.simplefilter(
            "always",
            ConvergenceWarning
        )

        model.fit(
            X,
            y,
            sample_weight=w
        )

    ni = int(
        np.max(
            model.n_iter_
        )
    )

    warned = any(
        issubclass(
            z.category,
            ConvergenceWarning
        )
        for z in ww
    )

    ok = (
        (not warned)
        and ni < max_iter
    )

    nz = int(
        np.count_nonzero(
            np.abs(
                model.coef_[0]
            )
            > 1e-12
        )
    )

    return (
        model,
        ok,
        ni,
        nz
    )


def tune_outer(
    data,
    train,
    l1s,
    Cs,
    outer
):

    store = {
        (
            float(l),
            float(C)
        ): []
        for l in l1s
        for C in Cs
    }

    conv = {
        (
            float(l),
            float(C)
        ): []
        for l in l1s
        for C in Cs
    }

    cache = {}
    di = []

    for ii, h in enumerate(
        train
    ):

        tr = [
            c
            for c in train
            if c != h
        ]

        print(
            f"    inner holdout {h}: "
            f"fitting full-universe grid..."
        )

        (
            X,
            y,
            _,
            _,
            _,
            w
        ) = b.pooled_data(
            data,
            tr
        )

        Xt = data[h].X

        (
            Xz,
            Xtz,
            _,
            _
        ) = b.standardize_train_test(
            X,
            Xt
        )

        for li, l1 in enumerate(
            l1s
        ):

            model = None

            for C in sorted(Cs):

                (
                    model,
                    ok,
                    ni,
                    nz
                ) = fit_model(
                    Xz,
                    y,
                    w,
                    l1,
                    C,
                    20260922
                    + 1000 * ii
                    + li,
                    2000,
                    True,
                    model
                )

                p = model.predict_proba(
                    Xtz
                )[:, 1]

                key = (
                    float(l1),
                    float(C)
                )

                store[key].append(
                    b.auc_binary(
                        data[h].y,
                        p
                    )
                )

                conv[key].append(
                    ok
                )

                cache[
                    (
                        h,
                        float(l1),
                        float(C)
                    )
                ] = p

                di.append(
                    {
                        "outer_holdout":
                            outer,
                        "inner_holdout":
                            h,
                        "l1_ratio":
                            float(l1),
                        "C":
                            float(C),
                        "converged":
                            ok,
                        "n_iter":
                            ni,
                        "nonzero_coefficients":
                            nz
                    }
                )

    rows = []
    cand = []

    for l1 in l1s:

        for C in Cs:

            key = (
                float(l1),
                float(C)
            )

            v = store[key]

            eligible = bool(
                all(
                    conv[key]
                )
            )

            a = float(
                np.mean(v)
            )

            rows.append(
                {
                    "outer_holdout":
                        outer,
                    "l1_ratio":
                        float(l1),
                    "C":
                        float(C),
                    "macro_inner_auc":
                        a,
                    "inner_auc_min":
                        float(
                            np.min(v)
                        ),
                    "inner_auc_max":
                        float(
                            np.max(v)
                        ),
                    "n_inner_cohorts":
                        len(v),
                    "all_inner_fits_converged":
                        eligible
                }
            )

            if eligible:

                # Highest AUROC.
                # Tie: smaller C.
                # Then larger l1 ratio.
                cand.append(
                    (
                        -a,
                        float(C),
                        -float(l1),
                        float(l1)
                    )
                )

    if not cand:

        raise RuntimeError(
            f"{outer}: no ElasticNet "
            "hyperparameter combination "
            "converged in all inner folds"
        )

    cand.sort()

    (
        neg,
        Cbest,
        _,
        lbest
    ) = cand[0]

    print(
        "    tuned ElasticNet: "
        f"l1_ratio={lbest:g}, "
        f"C={Cbest:g}, "
        f"inner macro AUC={-neg:.4f}"
    )

    return (
        (
            lbest,
            Cbest
        ),
        pd.DataFrame(rows),
        pd.DataFrame(di),
        cache
    )


def fit_outer(
    data,
    train,
    outer,
    l1,
    C,
    cache,
    alpha
):

    (
        X,
        y,
        _,
        _,
        _,
        w
    ) = b.pooled_data(
        data,
        train
    )

    Xt = data[outer].X

    (
        Xz,
        Xtz,
        _,
        _
    ) = b.standardize_train_test(
        X,
        Xt
    )

    (
        model,
        ok,
        ni,
        nz
    ) = fit_model(
        Xz,
        y,
        w,
        l1,
        C,
        20260922,
        12000,
        False,
        None
    )

    raw = model.predict_proba(
        Xtz
    )[:, 1]

    pp = []
    yy = []
    cc = []

    for h in train:

        pp.append(
            cache[
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
        object
    )

    cw = b.cohort_class_weights(
        yy,
        cc
    )

    (
        a,
        s,
        calok
    ) = calfit(
        pp,
        yy,
        cw,
        alpha
    )

    p = calapply(
        raw,
        a,
        s
    )

    info = {
        "l1_ratio":
            float(l1),

        "C":
            float(C),

        "converged":
            ok,

        "iterations":
            ni,

        "nonzero_coefficients":
            nz,

        "calibration_intercept_train_crossfit":
            a,

        "calibration_slope_train_crossfit":
            s,

        "calibration_converged":
            calok
    }

    coef = (
        model.coef_[0]
        .copy()
    )

    return (
        raw,
        p,
        info,
        coef
    )


def mrow(
    outer,
    y,
    raw,
    p
):

    (
        se,
        sp
    ) = b.sensitivity_specificity(
        y,
        p,
        0.5
    )

    (
        ci,
        cs,
        ok
    ) = b.evaluation_calibration(
        y,
        p
    )

    return {
        "outer_holdout":
            outer,

        "AUROC_raw":
            b.auc_binary(
                y,
                raw
            ),

        "AUROC_calibrated":
            b.auc_binary(
                y,
                p
            ),

        "AUPRC_calibrated":
            b.average_precision(
                y,
                p
            ),

        "balanced_Brier":
            b.balanced_brier(
                y,
                p
            ),

        "calibration_intercept_descriptive":
            ci,

        "calibration_slope_descriptive":
            cs,

        "calibration_fit_ok":
            ok,

        "sensitivity_at_0.5":
            se,

        "specificity_at_0.5":
            sp
    }


def main():

    t = time.time()

    lock = load_lock()

    b.verify_stage0()

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    (
        genes,
        data,
        _
    ) = b.load_development_data(
        lock
    )

    dev = list(
        lock[
            "development_accessions"
        ]
    )

    l1s = [
        float(x)
        for x in lock[
            "elasticnet_l1_ratio_grid"
        ]
    ]

    Cs = [
        10 ** float(x)
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

    tun = []
    dia = []
    met = []
    hyp = []
    pred = []
    coefs = []

    for oi, outer in enumerate(
        dev,
        1
    ):

        print(
            f"\n=== OUTER "
            f"{oi}/{len(dev)}: "
            f"hold out {outer} ==="
        )

        train = [
            c
            for c in dev
            if c != outer
        ]

        (
            best,
            tt,
            dd,
            cache
        ) = tune_outer(
            data,
            train,
            l1s,
            Cs,
            outer
        )

        tun.append(tt)
        dia.append(dd)

        (
            l1,
            C
        ) = best

        print(
            "  fitting final "
            "ElasticNet: "
            f"l1_ratio={l1:g}, "
            f"C={C:g}"
        )

        (
            raw,
            p,
            info,
            coef
        ) = fit_outer(
            data,
            train,
            outer,
            l1,
            C,
            cache,
            alpha
        )

        hyp.append(
            {
                "outer_holdout":
                    outer,
                **info
            }
        )

        met.append(
            mrow(
                outer,
                data[outer].y,
                raw,
                p
            )
        )

        for j in np.flatnonzero(
            np.abs(coef)
            > 1e-12
        ):

            coefs.append(
                {
                    "outer_holdout":
                        outer,

                    "gene":
                        str(
                            genes[j]
                        ),

                    "coefficient":
                        float(
                            coef[j]
                        )
                }
            )

        for i in range(
            len(
                data[outer].y
            )
        ):

            pred.append(
                {
                    "outer_holdout":
                        outer,

                    "sample_id":
                        data[
                            outer
                        ].samples[i],

                    "participant_id":
                        data[
                            outer
                        ].participant[i],

                    "case":
                        int(
                            data[
                                outer
                            ].y[i]
                        ),

                    "raw_probability":
                        float(
                            raw[i]
                        ),

                    "calibrated_probability":
                        float(
                            p[i]
                        )
                }
            )

        pd.concat(
            tun,
            ignore_index=True
        ).to_csv(
            OUT
            / "inner_tuning.csv.gz",
            index=False,
            compression="gzip"
        )

        pd.concat(
            dia,
            ignore_index=True
        ).to_csv(
            OUT
            / "optimizer_diagnostics.csv.gz",
            index=False,
            compression="gzip"
        )

        pd.DataFrame(
            hyp
        ).to_csv(
            OUT
            / "outer_hyperparameters.csv",
            index=False
        )

        pd.DataFrame(
            met
        ).to_csv(
            OUT
            / "outer_metrics.csv",
            index=False
        )

        pd.DataFrame(
            pred
        ).to_csv(
            OUT
            / "outer_predictions.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip"
        )

        pd.DataFrame(
            coefs
        ).to_csv(
            OUT
            / "nonzero_coefficients.csv.gz",
            index=False,
            compression="gzip"
        )

        print(
            "  checkpoint written "
            f"({(time.time() - t) / 60:.1f} "
            "min elapsed)"
        )

    M = pd.DataFrame(met)
    H = pd.DataFrame(hyp)

    D = pd.concat(
        dia,
        ignore_index=True
    )

    macro = pd.DataFrame(
        [
            {
                "AUROC_raw":
                    M[
                        "AUROC_raw"
                    ].mean(),

                "AUROC_calibrated":
                    M[
                        "AUROC_calibrated"
                    ].mean(),

                "AUPRC_calibrated":
                    M[
                        "AUPRC_calibrated"
                    ].mean(),

                "balanced_Brier":
                    M[
                        "balanced_Brier"
                    ].mean(),

                "sensitivity_at_0.5":
                    M[
                        "sensitivity_at_0.5"
                    ].mean(),

                "specificity_at_0.5":
                    M[
                        "specificity_at_0.5"
                    ].mean()
            }
        ]
    )

    macro.to_csv(
        OUT
        / "macro_metrics.csv",
        index=False
    )

    nonconv = (
        int(
            (
                ~D[
                    "converged"
                ].astype(bool)
            ).sum()
        )
        if len(D)
        else 0
    )

    R = [
        "# Stage-1C full-universe ElasticNet development report",
        "",
        (
            f"Protocol lock: "
            f"{lock.get('lock_version', '')} "
            f"({lock.get('lock_date', '')})"
        ),
        "",
        "External-expression loading status: **BLOCKED / NOT LOADED**",
        "",
        (
            "Only the six frozen development cohorts were loaded. "
            "The ElasticNet used the complete frozen "
            "10,998-gene universe."
        ),
        "",
        "## Tuned outer-fold ElasticNet models",
        "",
        b.write_markdown_table(H),
        "",
        "## Held-out development-cohort performance",
        "",
        b.write_markdown_table(M),
        "",
        "## Macro-average across the six outer held-out cohorts",
        "",
        b.write_markdown_table(macro),
        "",
        "## Convergence audit",
        "",
        (
            "Inner-grid fits not reporting convergence: "
            f"**{nonconv}**."
        ),
        "",
        (
            "Final outer fits converged: "
            f"**{int(H['converged'].sum())}/{len(H)}**."
        ),
        "",
        "## Gate",
        "",
        "`external_expression_loading_authorized` remains **false**.",
        "",
        "`external_performance_evaluation_authorized` remains **false**.",
        "",
        "No external cohort was loaded or evaluated in Stage 1C.",
        "",
        (
            f"Runtime: "
            f"{(time.time() - t) / 60:.1f} minutes."
        )
    ]

    out = (
        REPORTS
        / "stage1c_elasticnet_report.md"
    )

    out.write_text(
        "\n".join(R),
        encoding="utf-8"
    )

    print(
        "\n=== STAGE 1C COMPLETE ==="
    )

    print(
        f"Report: {out}"
    )

    print(
        "External cohorts remain locked."
    )


if __name__ == "__main__":
    main()
