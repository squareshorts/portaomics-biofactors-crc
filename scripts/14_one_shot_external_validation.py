import sys
import json
import argparse
import hashlib
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, confusion_matrix
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"
OUT = RESULTS / "external_validation"
MANIFEST_PATH = RESULTS / "manifests" / "sample_manifest.tsv"

LOCK_PATH = CONFIG / "final_external_lock.json"
AMEND_PATH = REPORTS / "final_external_protocol_amendment_A1.md"

EXPECTED_LOCK_HASH = "67B21BF5D360CA70B76033A3BD324624940169B3A0692E4E31A69639EA279255".lower()
EXPECTED_AMEND_HASH = "DC4DC1DC513DBBF7A01242A54CACE8F49FB6AC78718F82C2FCB63BAEB5C1C363".lower()

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(1<<20)
            if not b:
                break
            h.update(b)
    return h.hexdigest().lower()

def run_checks():
    if not LOCK_PATH.exists():
        sys.exit(f"Lock file missing: {LOCK_PATH}")
    with LOCK_PATH.open() as f:
        lock = json.load(f)
    
    if lock["status"] != "DEVELOPMENT_FROZEN_EXTERNAL_AUTHORIZED":
        sys.exit(f"Status is not DEVELOPMENT_FROZEN_EXTERNAL_AUTHORIZED: {lock['status']}")
    
    if not lock.get("external_expression_loading_authorized"):
        sys.exit("external_expression_loading_authorized is not true.")
    if not lock.get("external_performance_evaluation_authorized"):
        sys.exit("external_performance_evaluation_authorized is not true.")
        
    actual_lock_hash = sha256(LOCK_PATH)
    if actual_lock_hash != EXPECTED_LOCK_HASH:
        sys.exit(f"Lock hash mismatch.\nExpected: {EXPECTED_LOCK_HASH}\nGot:      {actual_lock_hash}")
        
    actual_amend_hash = sha256(AMEND_PATH)
    if actual_amend_hash != EXPECTED_AMEND_HASH:
        sys.exit(f"Amendment hash mismatch.\nExpected: {EXPECTED_AMEND_HASH}\nGot:      {actual_amend_hash}")
        
    for k, v in lock["artifact_hashes"].items():
        p = ROOT / v["path"]
        if not p.exists():
            sys.exit(f"Artifact missing: {p}")
        h = sha256(p)
        if h != v["sha256"]:
            sys.exit(f"Hash mismatch for {k}: expected {v['sha256']}, got {h}")
            
    for fm in lock["frozen_models"]:
        method = fm["method"]
        p = ROOT / "results" / "final_development" / f"{method.lower()}_final_model.json"
        with p.open() as f:
            j = json.load(f)
        for field in ["lambda_portability", "m", "C"]:
            if field in fm and fm[field] is not None:
                if not np.isclose(float(fm[field]), float(j[field]), atol=1e-5):
                    sys.exit(f"Hyperparameter mismatch for {method}: {field} {fm[field]} vs {j[field]}")

    return lock

def load_model(method):
    method_lower = method.lower()
    df = pd.read_csv(ROOT / "results" / "final_development" / f"{method_lower}_final_model.csv.gz")
    with open(ROOT / "results" / "final_development" / f"{method_lower}_final_model.json") as f:
        meta = json.load(f)
    return df, meta

def prepare_and_check_structural_availability(lock, models, external_accessions, manifest):
    with open(ROOT / "results" / "qc" / "gene_intersection.txt") as f:
        common_genes = [line.strip() for line in f if line.strip()]
        
    analysis_inputs = []
    
    for acc in external_accessions:
        expr_path = ROOT / "data" / "processed" / "rank_normal_scores" / f"{acc}.pkl.gz"
        expr = pd.read_pickle(expr_path).T
        
        analyses = []
        if acc == "GSE44076":
            primary_df = manifest[(manifest["accession"] == acc) & manifest["primary_include"]]
            stress_df = manifest[(manifest["accession"] == acc) & (manifest["phenotype"].isin(["tumor", "normal_healthy"]))]
            analyses.append(("primary", primary_df))
            analyses.append(("stress_test", stress_df))
        else:
            analyses.append(("primary", manifest[(manifest["accession"] == acc) & manifest["primary_include"]]))
            
        for analysis_type, df in analyses:
            if len(df) == 0:
                continue
                
            df = df[df["phenotype"].isin(["tumor", "normal_adjacent", "normal_healthy", "normal_other"])].copy()
            y_true = (df["phenotype"] == "tumor").astype(int).values
            sample_ids = df["sample_id"].values
            participant_ids = df["participant_id"].values
            
            # Fill missing participant IDs with sample IDs
            participant_ids = np.where(pd.isna(participant_ids) | (participant_ids == ""), sample_ids, participant_ids)
            
            if len(sample_ids) != len(set(sample_ids)):
                sys.exit(f"Duplicate sample IDs in {acc} {analysis_type}")
                
            if pd.isna(participant_ids).any() or (participant_ids == "").any():
                sys.exit(f"Missing participant IDs in {acc} {analysis_type} after fallback")
                
            if len(df["phenotype"].dropna()) != len(df):
                sys.exit(f"Missing phenotypes in {acc} {analysis_type}")
                
            if analysis_type == "primary" and len(np.unique(y_true)) < 2:
                sys.exit(f"Primary cohort {acc} lacks both classes")
                
            missing = [s for s in sample_ids if s not in expr.index]
            if missing:
                if acc == "GSE44076" and analysis_type == "stress_test":
                    sys.path.insert(0, str(ROOT / "scripts"))
                    from common import (
                        read_series_expression,
                        read_platform_annotation,
                        collapse_probes_to_genes,
                        rank_normal_scores,
                        RAW_GEO,
                        RAW_PLATFORM,
                    )
                    
                    reg = pd.read_csv(CONFIG / "datasets.csv")
                    platform = reg[reg["accession"] == acc]["platform"].values[0]
                    
                    raw_matrix_path = RAW_GEO / f"{acc}_series_matrix.txt.gz"
                    if not raw_matrix_path.exists():
                        sys.exit(f"Recovery raw matrix {raw_matrix_path} missing for stress test.")
                        
                    mat = read_series_expression(raw_matrix_path, use_samples=missing)
                    annot = read_platform_annotation(RAW_PLATFORM / f"{platform}.annot.gz")
                    gene, _ = collapse_probes_to_genes(mat, annot)
                    gene = gene.reindex(common_genes)
                    ranked = rank_normal_scores(gene).T
                    
                    expr = pd.concat([expr, ranked])
                else:
                    sys.exit(f"Missing samples {missing} in {acc} {analysis_type} and no recovery route exists.")
            
            expr_sub = expr.loc[sample_ids]
            for method, (m_df, m_meta) in models.items():
                missing_genes = [g for g in m_df['gene'] if g not in expr_sub.columns]
                if missing_genes:
                    sys.exit(f"Missing required model genes in {acc} {analysis_type}: {missing_genes}")
                    
            analysis_inputs.append({
                "acc": acc,
                "analysis_type": analysis_type,
                "df": df,
                "y_true": y_true,
                "participant_ids": participant_ids,
                "sample_ids": sample_ids,
                "expr_sub": expr_sub
            })
            
    return analysis_inputs

def score_samples(expr_df, model_df, model_meta):
    X = expr_df[model_df['gene']].values
    means = model_df['training_mean'].values
    sds = model_df['training_sd'].values
    coefs = model_df['coefficient_standardized'].values
    intercept = model_meta['intercept_standardized']
    
    Z = (X - means) / sds
    raw_scores = Z @ coefs + intercept
    
    cal_intercept = model_meta['calibration_intercept']
    cal_slope = model_meta['calibration_slope']
    calibrated_scores = expit(raw_scores * cal_slope + cal_intercept)
    
    return calibrated_scores

def eval_metrics(y_true, y_pred, threshold=0.5):
    if len(np.unique(y_true)) < 2:
        return {
            'AUROC': np.nan, 'AUPRC': np.nan, 'balanced_Brier': np.nan,
            'calibration_intercept': np.nan, 'calibration_slope': np.nan,
            'sensitivity_at_0.5': np.nan, 'specificity_at_0.5': np.nan
        }
        
    auc = roc_auc_score(y_true, y_pred)
    auprc = average_precision_score(y_true, y_pred)
    brier = brier_score_loss(y_true, y_pred)
    
    eps = 1e-15
    y_pred_clipped = np.clip(y_pred, eps, 1-eps)
    log_odds = logit(y_pred_clipped)
    
    if np.var(log_odds) > 1e-8:
        try:
            lr = LogisticRegression(penalty=None, solver='lbfgs')
            lr.fit(log_odds.reshape(-1, 1), y_true)
            cal_slope_eval = lr.coef_[0][0]
            cal_intercept_eval = lr.intercept_[0]
        except Exception:
            cal_slope_eval = np.nan
            cal_intercept_eval = np.nan
    else:
        cal_slope_eval = np.nan
        cal_intercept_eval = np.nan
        
    y_pred_bin = (y_pred >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_bin, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    
    return {
        'AUROC': auc,
        'AUPRC': auprc,
        'balanced_Brier': brier,
        'calibration_intercept': cal_intercept_eval,
        'calibration_slope': cal_slope_eval,
        'sensitivity_at_0.5': sens,
        'specificity_at_0.5': spec
    }

def get_bootstrap_metrics(y_true, y_pred, participant_ids, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    unique_participants = np.unique(participant_ids)
    
    aucs = []
    for _ in range(n_boot):
        boot_participants = rng.choice(unique_participants, size=len(unique_participants), replace=True)
        boot_idx = []
        for p in boot_participants:
            boot_idx.extend(np.where(participant_ids == p)[0])
        
        y_true_b = y_true[boot_idx]
        y_pred_b = y_pred[boot_idx]
        
        if len(np.unique(y_true_b)) == 2:
            aucs.append(roc_auc_score(y_true_b, y_pred_b))
        else:
            aucs.append(np.nan)
            
    aucs = np.array(aucs)
    aucs = aucs[~np.isnan(aucs)]
    if len(aucs) == 0:
        return np.nan, np.nan, np.nan
        
    return np.mean(aucs), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--recovery-preflight", action="store_true")
    args = parser.parse_args()

    lock = run_checks()
    
    if OUT.exists() and any(OUT.iterdir()) and not (args.preflight or args.recovery_preflight):
        sys.exit(f"Output directory {OUT} is not empty. Failsafe activated to prevent overwrite.")
    
    if args.preflight:
        print("FINAL EXTERNAL VALIDATION PREFLIGHT PASSED")
        sys.exit(0)

    manifest = pd.read_csv(MANIFEST_PATH, sep='\t')
    external_accessions = lock["external_accessions"]
    
    models = {}
    for method in ["PAFS", "M_ONLY", "COHORT_ADJ", "ELASTICNET"]:
        models[method] = load_model(method)
        
    # Structural check phase: process all data dependencies before any metrics are generated.
    analysis_inputs = prepare_and_check_structural_availability(lock, models, external_accessions, manifest)
    
    if args.recovery_preflight:
        print("EXTERNAL VALIDATION RECOVERY PREFLIGHT PASSED")
        sys.exit(0)
        
    OUT.mkdir(parents=True, exist_ok=True)

    predictions_list = []
    primary_metrics_list = []
    secondary_metrics_list = []
    bootstrap_list = []
    
    sample_counts = {}
    class_counts = {}
    participant_counts = {}
    
    BOOT_SEED = 20260924
    BOOT_REPS = 2000
    
    for inputs in analysis_inputs:
        acc = inputs["acc"]
        analysis_type = inputs["analysis_type"]
        df = inputs["df"]
        y_true = inputs["y_true"]
        participant_ids = inputs["participant_ids"]
        sample_ids = inputs["sample_ids"]
        expr_sub = inputs["expr_sub"]
        
        preds = {}
        for method, (m_df, m_meta) in models.items():
            preds[method] = score_samples(expr_sub, m_df, m_meta)
            
        if not np.allclose(preds["PAFS"], preds["M_ONLY"], atol=1e-7):
            sys.exit(f"Identity check failed for {acc} {analysis_type}: PAFS and M_ONLY predictions differ.")
            
        for method, (m_df, m_meta) in models.items():
            m_res = eval_metrics(y_true, preds[method])
            row = {
                "cohort": acc,
                "analysis": analysis_type,
                "method": method,
                **m_res
            }
            if method == "PAFS":
                primary_metrics_list.append(row)
            else:
                secondary_metrics_list.append(row)
            
            if method == "PAFS":
                b_mean, b_low, b_high = get_bootstrap_metrics(y_true, preds[method], participant_ids, n_boot=BOOT_REPS, seed=BOOT_SEED)
                bootstrap_list.append({
                    "cohort": acc,
                    "analysis": analysis_type,
                    "AUROC": m_res["AUROC"],
                    "boot_mean": b_mean,
                    "boot_2.5": b_low,
                    "boot_97.5": b_high
                })
                
        for i, sid in enumerate(sample_ids):
            predictions_list.append({
                "cohort": acc,
                "analysis": analysis_type,
                "sample_id": sid,
                "true_class": y_true[i],
                "PAFS": preds["PAFS"][i],
                "M_ONLY": preds["M_ONLY"][i],
                "COHORT_ADJ": preds["COHORT_ADJ"][i],
                "ELASTICNET": preds["ELASTICNET"][i]
            })
            
        sample_counts[f"{acc}_{analysis_type}"] = len(df)
        class_counts[f"{acc}_{analysis_type}"] = {"tumor": int(np.sum(y_true)), "normal": int(np.sum(1-y_true))}
        participant_counts[f"{acc}_{analysis_type}"] = len(np.unique(participant_ids))
            
    primary_df = pd.DataFrame(primary_metrics_list)
    primary_cohorts = primary_df[primary_df["analysis"] == "primary"]
    macro_auroc = primary_cohorts["AUROC"].mean()
    
    pd.DataFrame(primary_metrics_list).to_csv(OUT / "external_primary_metrics.csv", index=False)
    pd.DataFrame(secondary_metrics_list).to_csv(OUT / "external_secondary_metrics.csv", index=False)
    pd.DataFrame(predictions_list).to_csv(OUT / "external_predictions.csv.gz", index=False, compression="gzip")
    pd.DataFrame(bootstrap_list).to_csv(OUT / "external_bootstrap_primary.csv", index=False)
    
    id_check = []
    pred_df = pd.DataFrame(predictions_list)
    for (acc, analysis), grp in pred_df.groupby(["cohort", "analysis"]):
        diff = np.abs(grp["PAFS"] - grp["M_ONLY"])
        id_check.append({
            "cohort": acc,
            "analysis": analysis,
            "max_abs_diff": float(diff.max()),
            "mean_abs_diff": float(diff.mean())
        })
    pd.DataFrame(id_check).to_csv(OUT / "external_identity_check.csv", index=False)
    
    with open(OUT / "external_execution_manifest.json", "w") as f:
        json.dump({
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python_version": sys.version,
            "script_sha256": sha256(Path(__file__)),
            "lock_sha256": EXPECTED_LOCK_HASH,
            "amendment_sha256": EXPECTED_AMEND_HASH,
            "artifact_hashes": lock["artifact_hashes"],
            "external_cohorts": external_accessions,
            "sample_counts": sample_counts,
            "class_counts": class_counts,
            "participant_counts": participant_counts,
            "bootstrap_seed": BOOT_SEED,
            "bootstrap_replicates": BOOT_REPS,
            "primary_estimand": lock["primary_external_estimand"],
            "secondary_endpoints": lock["secondary_external_endpoints"],
            "confirmation": "No fitting/tuning was performed on external data. Frozen models used as-is.",
            "macro_mean_primary_AUROC": macro_auroc
        }, f, indent=4)
        
    with open(REPORTS / "final_external_validation.md", "w") as f:
        f.write("# Final External Validation\n\n")
        f.write(f"Macro-mean AUROC (Primary analyses, PAFS): {macro_auroc:.4f}\n")
        
    print("External validation complete.")

if __name__ == "__main__":
    main()
