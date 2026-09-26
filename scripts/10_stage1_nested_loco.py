from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import __version__ as scipy_version
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
RESULTS = ROOT / "results"
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"
MANIFEST_PATH = RESULTS / "manifests" / "sample_manifest.tsv"
GENE_UNIVERSE_PATH = RESULTS / "qc" / "gene_intersection.txt"
STAGE0_REPORT = REPORTS / "stage0_audit.md"
STAGE1_LOCK = CONFIG / "stage1_lock.json"
OUTDIR = RESULTS / "stage1"

METHODS = ["PAFS", "M_ONLY", "COHORT_ADJ"]
EPS = 1e-12


@dataclass
class CohortData:
    accession: str
    genes: np.ndarray
    samples: np.ndarray
    X: np.ndarray  # samples x genes, float32
    y: np.ndarray  # 0/1
    participant: np.ndarray


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_lock() -> dict:
    if not STAGE1_LOCK.exists():
        raise SystemExit(
            f"Missing {STAGE1_LOCK}. Stage 1 is not authorized until stage1_lock.json exists."
        )
    lock = json.loads(STAGE1_LOCK.read_text(encoding="utf-8-sig"))
    if lock.get("stage0_review_status") != "PASS":
        raise SystemExit("stage1_lock.json does not record Stage-0 PASS.")
    if not bool(lock.get("classifier_analysis_authorized", False)):
        raise SystemExit("Stage-1 classifier analysis is not authorized in stage1_lock.json.")
    if bool(lock.get("external_expression_loading_authorized", False)):
        raise SystemExit(
            "Safety gate failed: external_expression_loading_authorized must remain false in Stage 1."
        )
    if bool(lock.get("external_performance_evaluation_authorized", False)):
        raise SystemExit(
            "Safety gate failed: external_performance_evaluation_authorized must remain false in Stage 1."
        )
    return lock


def verify_stage0() -> None:
    if not STAGE0_REPORT.exists():
        raise SystemExit(f"Missing Stage-0 report: {STAGE0_REPORT}")
    txt = STAGE0_REPORT.read_text(encoding="utf-8", errors="replace")
    required = [
        "Overall inclusion-count status: **PASS**",
        "10,998",
        "classifier_analysis_authorized` remains **false**",
    ]
    missing = [x for x in required if x not in txt]
    if missing:
        raise SystemExit(f"Stage-0 report does not match reviewed PASS record; missing={missing}")


def truthy_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes", "y"])


def load_gene_universe(expected_n: int) -> np.ndarray:
    genes = [x.strip() for x in GENE_UNIVERSE_PATH.read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(genes) != expected_n:
        raise SystemExit(
            f"Frozen gene universe mismatch: observed {len(genes)}, expected {expected_n}."
        )
    if genes != sorted(genes):
        raise SystemExit("Frozen gene universe is not sorted as expected from Stage 0.")
    if len(set(genes)) != len(genes):
        raise SystemExit("Frozen gene universe contains duplicates.")
    return np.asarray(genes, dtype=object)


def read_ranked_expression(acc: str) -> pd.DataFrame:
    # Stage-0 on this machine was patched to gzip-pickle to avoid pyarrow.
    pkl = PROCESSED / "rank_normal_scores" / f"{acc}.pkl.gz"
    parquet = PROCESSED / "rank_normal_scores" / f"{acc}.parquet"
    if pkl.exists():
        return pd.read_pickle(pkl, compression="gzip")
    if parquet.exists():
        # Kept only for compatibility with an unpatched Stage-0 environment.
        return pd.read_parquet(parquet)
    raise SystemExit(f"Missing rank-normal expression for development cohort {acc}: {pkl}")


def load_development_data(lock: dict) -> Tuple[np.ndarray, Dict[str, CohortData], pd.DataFrame]:
    development = list(lock["development_accessions"])
    external = set(lock.get("external_accessions_locked", []))
    if set(development) & external:
        raise SystemExit("Stage-1 lock error: development and external accession sets overlap.")

    manifest = pd.read_csv(MANIFEST_PATH, sep="\t", dtype=str).fillna("")
    if not set(development).issubset(set(manifest["accession"])):
        raise SystemExit("Sample manifest is missing one or more frozen development cohorts.")

    genes = load_gene_universe(int(lock["feature_universe_n"]))
    out: Dict[str, CohortData] = {}

    print("Loading DEVELOPMENT cohorts only:")
    for acc in development:
        if acc in external:
            raise AssertionError("External cohort reached development loader.")
        m = manifest[
            (manifest["accession"] == acc) & truthy_series(manifest["primary_include"])
        ].copy()
        if m.empty:
            raise SystemExit(f"No primary samples in manifest for {acc}")
        m["case"] = pd.to_numeric(m["case"], errors="raise").astype(int)
        if set(m["case"].unique()) != {0, 1}:
            raise SystemExit(f"{acc}: expected both cases and controls.")

        expr = read_ranked_expression(acc)
        if set(genes) != set(expr.index.astype(str)):
            missing = sorted(set(genes) - set(expr.index.astype(str)))[:10]
            extra = sorted(set(expr.index.astype(str)) - set(genes))[:10]
            raise SystemExit(f"{acc}: gene-universe mismatch; missing={missing}, extra={extra}")
        expr = expr.reindex(genes)

        sample_order = m["sample_id"].tolist()
        missing_samples = [s for s in sample_order if s not in expr.columns]
        if missing_samples:
            raise SystemExit(f"{acc}: missing expression samples: {missing_samples[:10]}")
        expr = expr[sample_order]
        X = expr.to_numpy(dtype=np.float32).T
        if not np.isfinite(X).all():
            raise SystemExit(f"{acc}: non-finite values in Stage-0 rank-normal matrix.")

        participant = m["participant_id"].astype(str).to_numpy()
        empty = participant == ""
        participant[empty] = m.loc[empty, "sample_id"].astype(str).to_numpy()
        out[acc] = CohortData(
            accession=acc,
            genes=genes,
            samples=np.asarray(sample_order, dtype=object),
            X=X,
            y=m["case"].to_numpy(dtype=np.int8),
            participant=participant,
        )
        print(
            f"  {acc}: n={X.shape[0]} (case={int(out[acc].y.sum())}, "
            f"control={int((out[acc].y == 0).sum())}), genes={X.shape[1]}"
        )

    # Safety check: this process never reads any external expression path.
    print("External cohorts remain locked and were not loaded:")
    for acc in lock.get("external_accessions_locked", []):
        print(f"  {acc}")
    return genes, out, manifest


def cohort_class_weights(y: np.ndarray, cohort: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=int)
    cohort = np.asarray(cohort, dtype=object)
    levels = list(dict.fromkeys(cohort.tolist()))
    K = len(levels)
    w = np.zeros(len(y), dtype=float)
    for k in levels:
        idxk = cohort == k
        for cls in (0, 1):
            idx = idxk & (y == cls)
            n = int(idx.sum())
            if n == 0:
                raise ValueError(f"Cohort {k} lacks class {cls}; cannot apply frozen weighting.")
            w[idx] = 1.0 / (2.0 * K * n)
    # Preserve relative frozen weights while making mean weight 1 for numerical stability.
    w *= len(w) / w.sum()
    return w


def pooled_data(data: Mapping[str, CohortData], cohorts: Sequence[str], idx: np.ndarray | None = None):
    Xs, ys, cs, ss, ps = [], [], [], [], []
    for c in cohorts:
        d = data[c]
        Xs.append(d.X if idx is None else d.X[:, idx])
        ys.append(d.y)
        cs.extend([c] * len(d.y))
        ss.extend(d.samples.tolist())
        ps.extend(d.participant.tolist())
    X = np.vstack(Xs)
    y = np.concatenate(ys).astype(np.int8)
    cohort = np.asarray(cs, dtype=object)
    samples = np.asarray(ss, dtype=object)
    participant = np.asarray(ps, dtype=object)
    w = cohort_class_weights(y, cohort)
    return X, y, cohort, samples, participant, w


def auc_binary(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    pos = y == 1
    n1 = int(pos.sum())
    n0 = int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(score, method="average")
    u = float(r[pos].sum() - n1 * (n1 + 1) / 2.0)
    return u / (n1 * n0)


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    npos = int(y.sum())
    if npos == 0:
        return float("nan")
    order = np.argsort(-score, kind="mergesort")
    yy = y[order]
    tp = np.cumsum(yy)
    precision = tp / np.arange(1, len(y) + 1)
    return float(precision[yy == 1].sum() / npos)


def balanced_brier(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    e = (p - y) ** 2
    return 0.5 * (float(e[y == 1].mean()) + float(e[y == 0].mean()))


def sensitivity_specificity(y: np.ndarray, p: np.ndarray, threshold: float = 0.5) -> Tuple[float, float]:
    pred = np.asarray(p) >= threshold
    y = np.asarray(y, dtype=int)
    sens = float(pred[y == 1].mean()) if np.any(y == 1) else float("nan")
    spec = float((~pred[y == 0]).mean()) if np.any(y == 0) else float("nan")
    return sens, spec


def per_gene_signed_auc_effect(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    # X = samples x genes. scipy rankdata supports axis and uses average ties.
    y = np.asarray(y, dtype=int)
    pos = y == 1
    n1 = int(pos.sum())
    n0 = int((~pos).sum())
    if n1 == 0 or n0 == 0:
        raise ValueError("Both classes are required for per-gene AUC.")
    ranks = rankdata(X, axis=0, method="average")
    u = ranks[pos].sum(axis=0) - n1 * (n1 + 1) / 2.0
    auc = u / (n1 * n0)
    return 2.0 * (auc - 0.5)


def design_matrices(y: np.ndarray, cohort: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    y = np.asarray(y, dtype=float)
    cohort = np.asarray(cohort, dtype=object)
    levels = list(dict.fromkeys(cohort.tolist()))
    reduced = np.column_stack([np.ones(len(y)), y])
    dummies = []
    for level in levels[1:]:
        dummies.append((cohort == level).astype(float))
    full = reduced if not dummies else np.column_stack([reduced] + dummies)
    return reduced, full


def weighted_fit_summary(Y: np.ndarray, D: np.ndarray, w: np.ndarray):
    # Solve weighted least squares for all genes simultaneously without allocating residual n x g matrices.
    WD = D * w[:, None]
    XtWX = D.T @ WD
    XtWY = D.T @ (w[:, None] * Y)
    inv = np.linalg.pinv(XtWX, rcond=1e-12)
    B = inv @ XtWY
    ywy = (w[:, None] * (Y ** 2)).sum(axis=0)
    sse = ywy - np.sum(B * XtWY, axis=0)
    sse = np.maximum(sse, EPS)
    return B, sse, inv


def percentile_rank(v: np.ndarray, higher_better: bool = True) -> np.ndarray:
    x = np.asarray(v, dtype=float)
    if not higher_better:
        x = -x
    r = rankdata(x, method="average")
    return r / len(r)


def feature_scores(data: Mapping[str, CohortData], cohorts: Sequence[str], genes: np.ndarray) -> pd.DataFrame:
    effects = []
    for c in cohorts:
        d = data[c]
        effects.append(per_gene_signed_auc_effect(d.X, d.y))
    E = np.vstack(effects)
    med = np.median(E, axis=0)
    M = np.abs(med)
    med_sign = np.sign(med)
    signs = np.sign(E)
    S = np.mean(signs == med_sign[None, :], axis=0).astype(float)
    S[med_sign == 0] = 0.5  # frozen neutral treatment of exact-zero median effects
    H = np.median(np.abs(E - med[None, :]), axis=0)

    X, y, cohort, _, _, w = pooled_data(data, cohorts)
    reduced, full = design_matrices(y, cohort)
    _, sse_red, _ = weighted_fit_summary(X, reduced, w)
    B_full, sse_full, inv_full = weighted_fit_summary(X, full, w)
    C = 1.0 - sse_full / sse_red
    C = np.clip(C, 0.0, 1.0)

    # Cohort-adjusted disease association statistic: |weighted t| for the disease term.
    dof = max(1, len(y) - full.shape[1])
    sigma2 = sse_full / dof
    se_disease = np.sqrt(np.maximum(inv_full[1, 1] * sigma2, EPS))
    disease_t = B_full[1, :] / se_disease
    cohort_adj = np.abs(disease_t)

    p_M = percentile_rank(M, True)
    p_S = percentile_rank(S, True)
    p_H = percentile_rank(H, False)
    p_C = percentile_rank(C, False)
    pafs = (p_M + p_S + p_H + p_C) / 4.0

    return pd.DataFrame(
        {
            "gene": genes,
            "median_signed_auc_effect": med,
            "M_disease_magnitude": M,
            "S_directional_concordance": S,
            "H_effect_MAD": H,
            "C_partial_R2_cohort_given_disease": C,
            "cohort_adjusted_abs_t": cohort_adj,
            "rank_M": p_M,
            "rank_S": p_S,
            "rank_low_H": p_H,
            "rank_low_C": p_C,
            "PAFS": pafs,
        }
    )


def ranking_indices(scores: pd.DataFrame, method: str) -> np.ndarray:
    if method == "PAFS":
        v = scores["PAFS"].to_numpy(float)
    elif method == "M_ONLY":
        v = scores["M_disease_magnitude"].to_numpy(float)
    elif method == "COHORT_ADJ":
        v = scores["cohort_adjusted_abs_t"].to_numpy(float)
    else:
        raise ValueError(method)
    # Stable sorting preserves alphabetical gene order on exact ties.
    return np.argsort(-v, kind="mergesort")


def standardize_train_test(Xtr: np.ndarray, Xte: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mu = Xtr.mean(axis=0, dtype=np.float64)
    sd = Xtr.std(axis=0, dtype=np.float64)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (Xtr - mu) / sd, (Xte - mu) / sd, mu, sd


def ridge_logistic_fit(
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    C: float,
    start: np.ndarray | None = None,
) -> Tuple[np.ndarray, bool, int, float]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    sw = float(w.sum())
    n, p = X.shape
    lam = 1.0 / float(C)

    if start is None or len(start) != p + 1:
        prev = np.clip(np.average(y, weights=w), 1e-5, 1 - 1e-5)
        theta0 = np.zeros(p + 1, dtype=float)
        theta0[0] = logit(prev)
    else:
        theta0 = np.asarray(start, dtype=float).copy()

    def fun_grad(theta):
        b0 = theta[0]
        b = theta[1:]
        eta = b0 + X @ b
        # Stable logistic negative log likelihood.
        loss = float(np.sum(w * (np.logaddexp(0.0, eta) - y * eta)) / sw)
        loss += 0.5 * lam * float(b @ b)
        p1 = expit(eta)
        r = w * (p1 - y) / sw
        g0 = float(r.sum())
        gb = X.T @ r + lam * b
        grad = np.concatenate([[g0], gb])
        return loss, grad

    res = minimize(
        lambda th: fun_grad(th)[0],
        theta0,
        jac=lambda th: fun_grad(th)[1],
        method="L-BFGS-B",
        options={"maxiter": 120, "ftol": 1e-9, "gtol": 1e-6, "maxls": 30},
    )
    return res.x, bool(res.success), int(res.nit), float(res.fun)


def predict_logistic(theta: np.ndarray, X: np.ndarray) -> np.ndarray:
    return expit(theta[0] + np.asarray(X, dtype=float) @ theta[1:])


def fit_ridge_grid(
    Xtr: np.ndarray,
    ytr: np.ndarray,
    wtr: np.ndarray,
    Xte: np.ndarray,
    C_grid: Sequence[float],
) -> Tuple[Dict[float, np.ndarray], List[dict]]:
    # Strong -> weak regularization, warm-started.
    preds: Dict[float, np.ndarray] = {}
    diagnostics = []
    theta = None
    for C in sorted(C_grid):
        theta, ok, nit, obj = ridge_logistic_fit(Xtr, ytr, wtr, C, start=theta)
        preds[float(C)] = predict_logistic(theta, Xte)
        diagnostics.append({"C": float(C), "converged": ok, "iterations": nit, "objective": obj})
    return preds, diagnostics


def parse_C_grid(lock: dict) -> List[float]:
    if "ridge_C_grid_log10" in lock:
        return [10.0 ** float(x) for x in lock["ridge_C_grid_log10"]]
    raise SystemExit("stage1_lock.json is missing ridge_C_grid_log10")


def tune_one_outer(
    data: Mapping[str, CohortData],
    genes: np.ndarray,
    train_cohorts: Sequence[str],
    m_grid: Sequence[int],
    C_grid: Sequence[float],
    outer_acc: str,
) -> Tuple[Dict[str, Tuple[int, float]], pd.DataFrame, pd.DataFrame, Dict[str, pd.DataFrame]]:
    # rows keyed by method,m,C and aggregate AUC over inner held-out cohorts.
    auc_store = {
        (method, int(m), float(C)): []
        for method in METHODS for m in m_grid for C in C_grid
    }
    diag_rows = []
    score_cache: Dict[str, pd.DataFrame] = {}

    for inner_holdout in train_cohorts:
        inner_train = [c for c in train_cohorts if c != inner_holdout]
        print(f"    inner holdout {inner_holdout}: scoring {len(inner_train)} training cohorts...")
        scores = feature_scores(data, inner_train, genes)
        score_cache[inner_holdout] = scores
        rankings = {method: ranking_indices(scores, method) for method in METHODS}
        dtest = data[inner_holdout]

        for method in METHODS:
            ranking = rankings[method]
            for m in m_grid:
                idx = ranking[: int(m)]
                Xtr, ytr, ctr, _, _, wtr = pooled_data(data, inner_train, idx)
                Xte = dtest.X[:, idx]
                Xtrz, Xtez, _, _ = standardize_train_test(Xtr, Xte)
                preds, diagnostics = fit_ridge_grid(Xtrz, ytr, wtr, Xtez, C_grid)
                for C, p in preds.items():
                    auc_store[(method, int(m), float(C))].append(auc_binary(dtest.y, p))
                for d in diagnostics:
                    diag_rows.append(
                        {
                            "outer_holdout": outer_acc,
                            "inner_holdout": inner_holdout,
                            "method": method,
                            "m": int(m),
                            **d,
                        }
                    )

    tuning_rows = []
    best = {}
    for method in METHODS:
        candidates = []
        for m in m_grid:
            for C in C_grid:
                vals = auc_store[(method, int(m), float(C))]
                macro = float(np.mean(vals))
                tuning_rows.append(
                    {
                        "outer_holdout": outer_acc,
                        "method": method,
                        "m": int(m),
                        "C": float(C),
                        "macro_inner_auc": macro,
                        "inner_auc_min": float(np.min(vals)),
                        "inner_auc_max": float(np.max(vals)),
                        "n_inner_cohorts": len(vals),
                    }
                )
                # max AUC, then smaller signature, then smaller C (stronger regularization)
                candidates.append((-macro, int(m), float(C)))
        candidates.sort()
        neg, m_best, C_best = candidates[0]
        best[method] = (m_best, C_best)
        print(f"    tuned {method}: m={m_best}, C={C_best:.6g}, inner macro AUC={-neg:.4f}")

    return best, pd.DataFrame(tuning_rows), pd.DataFrame(diag_rows), score_cache


def calibration_fit(raw_p: np.ndarray, y: np.ndarray, w: np.ndarray, nonnegative_slope: bool = False) -> Tuple[float, float, bool]:
    p = np.clip(np.asarray(raw_p, dtype=float), 1e-6, 1 - 1e-6)
    x = logit(p)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    sw = float(w.sum())

    def fg(theta):
        eta = theta[0] + theta[1] * x
        pp = expit(eta)
        loss = float(np.sum(w * (np.logaddexp(0.0, eta) - y * eta)) / sw)
        r = w * (pp - y) / sw
        grad = np.array([r.sum(), np.sum(r * x)], dtype=float)
        return loss, grad

    bounds = [(None, None), (0.0, None)] if nonnegative_slope else None
    res = minimize(
        lambda th: fg(th)[0],
        np.array([0.0, 1.0]),
        jac=lambda th: fg(th)[1],
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 100, "ftol": 1e-10, "gtol": 1e-7},
    )
    return float(res.x[0]), float(res.x[1]), bool(res.success)


def apply_calibration(p: np.ndarray, a: float, b: float) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    return expit(a + b * logit(p))


def crossfit_for_calibration(
    data: Mapping[str, CohortData],
    genes: np.ndarray,
    train_cohorts: Sequence[str],
    method: str,
    m: int,
    C: float,
    score_cache: Mapping[str, pd.DataFrame] | None = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    pred_all, y_all, cohort_all = [], [], []
    for holdout in train_cohorts:
        subtrain = [c for c in train_cohorts if c != holdout]
        scores = score_cache[holdout] if score_cache is not None and holdout in score_cache else feature_scores(data, subtrain, genes)
        rank = ranking_indices(scores, method)
        idx = rank[:m]
        Xtr, ytr, ctr, _, _, wtr = pooled_data(data, subtrain, idx)
        Xte = data[holdout].X[:, idx]
        Xtrz, Xtez, _, _ = standardize_train_test(Xtr, Xte)
        theta, _, _, _ = ridge_logistic_fit(Xtrz, ytr, wtr, C)
        pred_all.append(predict_logistic(theta, Xtez))
        y_all.append(data[holdout].y)
        cohort_all.extend([holdout] * len(data[holdout].y))
    p = np.concatenate(pred_all)
    y = np.concatenate(y_all)
    cohort = np.asarray(cohort_all, dtype=object)
    w = cohort_class_weights(y, cohort)
    return p, y, w


def fit_outer_model(
    data: Mapping[str, CohortData],
    genes: np.ndarray,
    train_cohorts: Sequence[str],
    outer_acc: str,
    method: str,
    m: int,
    C: float,
    outer_scores: pd.DataFrame,
    calibration_score_cache: Mapping[str, pd.DataFrame] | None = None,
):
    rank = ranking_indices(outer_scores, method)
    idx = rank[:m]
    Xtr, ytr, ctr, _, _, wtr = pooled_data(data, train_cohorts, idx)
    dtest = data[outer_acc]
    Xte = dtest.X[:, idx]
    Xtrz, Xtez, mu, sd = standardize_train_test(Xtr, Xte)
    theta, ok, nit, obj = ridge_logistic_fit(Xtrz, ytr, wtr, C)
    raw = predict_logistic(theta, Xtez)

    # Frozen cross-fitted recalibration using training cohorts only.
    cf_p, cf_y, cf_w = crossfit_for_calibration(
        data, genes, train_cohorts, method, m, C, score_cache=calibration_score_cache
    )
    cal_a, cal_b, cal_ok = calibration_fit(cf_p, cf_y, cf_w, nonnegative_slope=True)
    calibrated = apply_calibration(raw, cal_a, cal_b)

    selected = outer_scores.iloc[idx].copy()
    selected.insert(0, "outer_holdout", outer_acc)
    selected.insert(1, "method", method)
    selected.insert(2, "selected_rank", np.arange(1, len(selected) + 1))
    selected.insert(3, "m", int(m))
    selected.insert(4, "C", float(C))

    info = {
        "converged": ok,
        "iterations": nit,
        "objective": obj,
        "calibration_intercept_train_crossfit": cal_a,
        "calibration_slope_train_crossfit": cal_b,
        "calibration_converged": cal_ok,
    }
    return raw, calibrated, selected, info


def evaluation_calibration(y: np.ndarray, p: np.ndarray) -> Tuple[float, float, bool]:
    # Descriptive held-out calibration intercept/slope; never fed back into fitting.
    w = np.ones(len(y), dtype=float)
    return calibration_fit(p, y, w, nonnegative_slope=False)


def metrics_row(outer: str, method: str, y: np.ndarray, raw: np.ndarray, p: np.ndarray) -> dict:
    sens, spec = sensitivity_specificity(y, p, 0.5)
    ci, cs, ok = evaluation_calibration(y, p)
    return {
        "outer_holdout": outer,
        "method": method,
        "n": len(y),
        "cases": int(np.sum(y == 1)),
        "controls": int(np.sum(y == 0)),
        "AUROC_raw": auc_binary(y, raw),
        "AUROC_calibrated": auc_binary(y, p),
        "AUPRC_calibrated": average_precision(y, p),
        "balanced_Brier": balanced_brier(y, p),
        "calibration_intercept_descriptive": ci,
        "calibration_slope_descriptive": cs,
        "calibration_fit_ok": ok,
        "sensitivity_at_0.5": sens,
        "specificity_at_0.5": spec,
    }


def bootstrap_auc_delta(
    y: np.ndarray,
    p1: np.ndarray,
    p0: np.ndarray,
    participant: np.ndarray,
    n_boot: int = 2000,
    seed: int = 20260922,
) -> Tuple[float, float, float, float]:
    # Cluster bootstrap by participant ID when available; singleton samples are their own clusters.
    rng = np.random.default_rng(seed)
    participant = np.asarray(participant, dtype=object)
    clusters = np.unique(participant)
    index_by_cluster = {c: np.flatnonzero(participant == c) for c in clusters}
    obs = auc_binary(y, p1) - auc_binary(y, p0)
    vals = []
    attempts = 0
    while len(vals) < n_boot and attempts < n_boot * 10:
        attempts += 1
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        idx = np.concatenate([index_by_cluster[c] for c in sampled])
        yy = y[idx]
        if len(np.unique(yy)) < 2:
            continue
        vals.append(auc_binary(yy, p1[idx]) - auc_binary(yy, p0[idx]))
    if len(vals) < max(100, n_boot // 2):
        raise RuntimeError("Too few valid bootstrap replicates for AUC difference.")
    arr = np.asarray(vals, dtype=float)
    return obs, float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975)), float(arr.std(ddof=1))


def dersimonian_laird(effects: np.ndarray, ses: np.ndarray) -> dict:
    effects = np.asarray(effects, dtype=float)
    variances = np.maximum(np.asarray(ses, dtype=float) ** 2, 1e-12)
    w = 1.0 / variances
    mu_fixed = float(np.sum(w * effects) / np.sum(w))
    Q = float(np.sum(w * (effects - mu_fixed) ** 2))
    k = len(effects)
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    wr = 1.0 / (variances + tau2)
    mu = float(np.sum(wr * effects) / np.sum(wr))
    se = math.sqrt(1.0 / float(np.sum(wr)))
    i2 = max(0.0, (Q - (k - 1)) / Q) * 100.0 if Q > 0 else 0.0
    return {
        "pooled_delta_AUROC_DL": mu,
        "ci95_low": mu - 1.96 * se,
        "ci95_high": mu + 1.96 * se,
        "tau2": tau2,
        "I2_percent": i2,
        "Q": Q,
        "k": k,
    }


def write_markdown_table(df: pd.DataFrame, digits: int = 4) -> str:
    if df.empty:
        return "(none)"
    cols = [str(c) for c in df.columns]
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return ""
            return f"{float(v):.{digits}f}"
        return str(v).replace("|", "\\|").replace("\n", " ")
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for row in df.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(fmt(v) for v in row) + " |")
    return "\n".join(lines)


def main() -> None:
    t0 = time.time()
    lock = load_lock()
    verify_stage0()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "feature_scores").mkdir(parents=True, exist_ok=True)

    genes, data, manifest = load_development_data(lock)
    development = list(lock["development_accessions"])
    m_grid = [int(x) for x in lock["signature_sizes"]]
    C_grid = parse_C_grid(lock)

    # Input provenance. No external expression file is hashed or opened.
    provenance = {
        "stage1_lock_sha256": sha256(STAGE1_LOCK),
        "stage0_report_sha256": sha256(STAGE0_REPORT),
        "sample_manifest_sha256": sha256(MANIFEST_PATH),
        "gene_universe_sha256": sha256(GENE_UNIVERSE_PATH),
        "development_expression": {},
        "external_expression_loaded": False,
        "methods": METHODS,
        "software": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy_version},
    }
    for acc in development:
        pkl = PROCESSED / "rank_normal_scores" / f"{acc}.pkl.gz"
        parquet = PROCESSED / "rank_normal_scores" / f"{acc}.parquet"
        p = pkl if pkl.exists() else parquet
        provenance["development_expression"][acc] = {"path": str(p.relative_to(ROOT)), "sha256": sha256(p)}
    (OUTDIR / "run_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    all_tuning = []
    all_diag = []
    all_metrics = []
    all_predictions = []
    all_selected = []
    all_hyper = []
    outer_rankings: Dict[str, Dict[str, np.ndarray]] = {}

    for oi, outer in enumerate(development, start=1):
        print(f"\n=== OUTER {oi}/{len(development)}: hold out {outer} ===")
        train_cohorts = [c for c in development if c != outer]
        best, tuning, diag, inner_score_cache = tune_one_outer(data, genes, train_cohorts, m_grid, C_grid, outer)
        all_tuning.append(tuning)
        all_diag.append(diag)

        print("  scoring full outer-training partition...")
        outer_scores = feature_scores(data, train_cohorts, genes)
        outer_rankings[outer] = {
            method: ranking_indices(outer_scores, method) for method in METHODS
        }
        outer_scores.insert(0, "outer_holdout", outer)
        outer_scores.to_csv(OUTDIR / "feature_scores" / f"{outer}.csv.gz", index=False, compression="gzip")

        dtest = data[outer]
        for method in METHODS:
            m, C = best[method]
            print(f"  fitting final {method}: m={m}, C={C:.6g}")
            raw, p, selected, info = fit_outer_model(
                data, genes, train_cohorts, outer, method, m, C,
                outer_scores.drop(columns=["outer_holdout"]),
                calibration_score_cache=inner_score_cache,
            )
            all_selected.append(selected)
            all_metrics.append(metrics_row(outer, method, dtest.y, raw, p))
            all_hyper.append({"outer_holdout": outer, "method": method, "m": m, "C": C, **info})
            for i in range(len(dtest.y)):
                all_predictions.append(
                    {
                        "outer_holdout": outer,
                        "method": method,
                        "sample_id": dtest.samples[i],
                        "participant_id": dtest.participant[i],
                        "case": int(dtest.y[i]),
                        "raw_probability": float(raw[i]),
                        "calibrated_probability": float(p[i]),
                    }
                )

        # Checkpoint after every outer fold.
        pd.concat(all_tuning, ignore_index=True).to_csv(OUTDIR / "inner_tuning.csv.gz", index=False, compression="gzip")
        pd.concat(all_diag, ignore_index=True).to_csv(OUTDIR / "optimizer_diagnostics.csv.gz", index=False, compression="gzip")
        pd.DataFrame(all_metrics).to_csv(OUTDIR / "outer_metrics.csv", index=False)
        pd.DataFrame(all_hyper).to_csv(OUTDIR / "outer_hyperparameters.csv", index=False)
        pd.DataFrame(all_predictions).to_csv(OUTDIR / "outer_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
        pd.concat(all_selected, ignore_index=True).to_csv(OUTDIR / "selected_features.csv.gz", index=False, compression="gzip")
        print(f"  checkpoint written ({time.time() - t0:.1f} s elapsed)")

    metrics = pd.DataFrame(all_metrics)
    predictions = pd.DataFrame(all_predictions)
    hyper = pd.DataFrame(all_hyper)

    # Primary comparison: PAFS versus M_ONLY.
    delta_rows = []
    for j, outer in enumerate(development):
        dtest = data[outer]
        sub = predictions[predictions["outer_holdout"] == outer]
        p_pafs = sub[sub["method"] == "PAFS"].set_index("sample_id").reindex(dtest.samples)["calibrated_probability"].to_numpy(float)
        p_m = sub[sub["method"] == "M_ONLY"].set_index("sample_id").reindex(dtest.samples)["calibrated_probability"].to_numpy(float)
        obs, lo, hi, se = bootstrap_auc_delta(
            dtest.y, p_pafs, p_m, dtest.participant,
            n_boot=2000, seed=20260922 + j,
        )
        delta_rows.append(
            {
                "outer_holdout": outer,
                "delta_AUROC_PAFS_minus_M_ONLY": obs,
                "bootstrap_ci95_low": lo,
                "bootstrap_ci95_high": hi,
                "bootstrap_SE": se,
            }
        )
    delta = pd.DataFrame(delta_rows)
    delta.to_csv(OUTDIR / "primary_auc_delta_bootstrap.csv", index=False)
    meta = dersimonian_laird(
        delta["delta_AUROC_PAFS_minus_M_ONLY"].to_numpy(float),
        delta["bootstrap_SE"].to_numpy(float),
    )
    pd.DataFrame([meta]).to_csv(OUTDIR / "primary_auc_delta_meta.csv", index=False)

    macro = (
        metrics.groupby("method", sort=False)[
            ["AUROC_raw", "AUROC_calibrated", "AUPRC_calibrated", "balanced_Brier", "sensitivity_at_0.5", "specificity_at_0.5"]
        ]
        .mean()
        .reset_index()
    )
    macro.to_csv(OUTDIR / "macro_metrics.csv", index=False)

    stability_rows = []
    N = len(genes)
    for method in METHODS:
        for m in m_grid:
            jac = []
            kun = []
            for i in range(len(development)):
                a = set(outer_rankings[development[i]][method][:m].tolist())
                for j in range(i + 1, len(development)):
                    b = set(outer_rankings[development[j]][method][:m].tolist())
                    r = len(a & b)
                    jac.append(r / len(a | b))
                    denom = m * (N - m)
                    kun.append((r * N - m * m) / denom if denom > 0 else float("nan"))
            stability_rows.append({
                "method": method,
                "m": int(m),
                "mean_pairwise_Jaccard": float(np.mean(jac)),
                "median_pairwise_Jaccard": float(np.median(jac)),
                "mean_pairwise_Kuncheva": float(np.mean(kun)),
                "median_pairwise_Kuncheva": float(np.median(kun)),
                "n_outer_fold_pairs": len(jac),
            })
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(OUTDIR / "feature_stability.csv", index=False)

    convergence = hyper.groupby("method")["converged"].agg(["sum", "count"]).reset_index()
    convergence["all_outer_final_fits_converged"] = convergence["sum"] == convergence["count"]
    optimizer_diag = pd.concat(all_diag, ignore_index=True)
    inner_nonconverged = int((~optimizer_diag["converged"].astype(bool)).sum()) if not optimizer_diag.empty else 0

    report = []
    report.append("# Stage-1 development-only nested LOCO report")
    report.append("")
    report.append(f"Protocol lock: {lock.get('lock_version', '')} ({lock.get('lock_date', '')})")
    report.append("")
    report.append("External-expression loading status: **BLOCKED / NOT LOADED**")
    report.append("")
    report.append(
        "This stage used only the six frozen development cohorts. GSE156451, GSE106582, and GSE44076 were not loaded by the modelling process."
    )
    report.append("")
    report.append("## Tuned outer-fold models")
    report.append("")
    report.append(write_markdown_table(hyper[["outer_holdout", "method", "m", "C", "converged", "calibration_intercept_train_crossfit", "calibration_slope_train_crossfit"]]))
    report.append("")
    report.append("## Held-out development-cohort performance")
    report.append("")
    report.append(write_markdown_table(metrics[["outer_holdout", "method", "AUROC_calibrated", "AUPRC_calibrated", "balanced_Brier", "calibration_slope_descriptive", "sensitivity_at_0.5", "specificity_at_0.5"]]))
    report.append("")
    report.append("## Macro-average across the six outer held-out cohorts")
    report.append("")
    report.append(write_markdown_table(macro))
    report.append("")
    report.append("## Feature-ranking stability across outer training partitions")
    report.append("")
    report.append(write_markdown_table(stability))
    report.append("")
    report.append("## Primary comparison: PAFS minus disease-magnitude-only AUROC")
    report.append("")
    report.append(write_markdown_table(delta))
    report.append("")
    report.append("Random-effects pooling (DerSimonian-Laird; descriptive with six cohorts):")
    report.append("")
    report.append(write_markdown_table(pd.DataFrame([meta])))
    report.append("")
    report.append("## Convergence audit")
    report.append("")
    report.append(write_markdown_table(convergence))
    report.append("")
    report.append(f"Inner-grid optimizer fits not reporting convergence: **{inner_nonconverged}**.")
    report.append("")
    report.append("## Gate")
    report.append("")
    report.append("`external_expression_loading_authorized` remains **false**.")
    report.append("")
    report.append("`external_performance_evaluation_authorized` remains **false**.")
    report.append("")
    report.append("Do not run locked external validation until this development-only report and the selected pipeline are reviewed and frozen.")
    report.append("")
    report.append(f"Runtime: {(time.time() - t0) / 60.0:.1f} minutes.")
    (REPORTS / "stage1_development_report.md").write_text("\n".join(report), encoding="utf-8")

    print("\n=== STAGE 1 COMPLETE ===")
    print(f"Report: {REPORTS / 'stage1_development_report.md'}")
    print("External cohorts remain locked.")


if __name__ == "__main__":
    main()
