"""
Export per-gene SHAP values for ALL LOEUF protein-coding genes in the DB.

Expanded version of export_gene_shap.py:
  - Loads LOEUF / pLoF / expression features from humanproof.db for all ~17.7K genes
  - Runs inference with the existing trained XGBoost models (no retraining)
  - Exports SHAP for every DB gene whose ENSG appears in the Open Targets parquet
  - Saves to data/safety_model_output/gene_shap_all.json

Usage:
    ~/miniconda3/bin/python data/export_shap_all_genes.py
"""

import json
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import shap
import xgboost as xgb

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "data"))

import re

from safety_model import (  # noqa: E402
    load_opentargets,
    log_transform_pvals,
    get_feature_cols,
    OUT_DIR,
    OUT_DIR_DRUGGED,
    DB_PATH,
)
from export_gene_shap import (  # noqa: E402
    READABLE_FEATURES,
    get_feature_group,
    _sigmoid,
)

# Note: OUT_DIR is still imported for the output path (gene_shap_all.json)


# ─── Expanded feature loaders (read all genes from DB) ────────────────────────

def load_expanded_loeuf() -> pd.DataFrame:
    """LOEUF scores for all DB genes → keyed by targetId (Ensembl)."""
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql(
        "SELECT gene_symbol, ensembl_id AS targetId, loeuf_score "
        "FROM gene_dosage_sensitivity WHERE ensembl_id != ''",
        conn,
    )
    conn.close()
    print(f"  Expanded LOEUF: {len(df):,} genes")
    return df


def load_expanded_plof() -> pd.DataFrame:
    """Gene-level pLoF aggregates for all DB genes → keyed by gene_symbol.

    Uses phenotype_category (lowercase, e.g. 'cardiovascular') to match the
    feature names the models were trained with.
    """
    conn = sqlite3.connect(DB_PATH)

    # Overall gene-level stats
    base = pd.read_sql("""
        SELECT
            gene_symbol,
            MIN(p_value)                                      AS plof_min_pval,
            SUM(CASE WHEN p_value < 0.00000005 THEN 1 ELSE 0 END) AS plof_n_sig,
            MAX(ABS(beta))                                    AS plof_max_abs_beta,
            COUNT(*)                                          AS plof_n_phenotypes
        FROM plof_associations
        GROUP BY gene_symbol
    """, conn)

    # Per-organ stats — use phenotype_category (lowercase) to match model features
    organ = pd.read_sql("""
        SELECT
            gene_symbol,
            phenotype_category,
            MIN(p_value)     AS min_p,
            MAX(ABS(beta))   AS max_beta
        FROM plof_associations
        GROUP BY gene_symbol, phenotype_category
    """, conn)
    conn.close()

    # Pivot organs → wide: ("min_p","cardiovascular") → "plof_cardiovascular_min_p"
    # All phenotype categories from the DB are included — no hardcoded filter.
    organ_wide = organ.pivot_table(
        index="gene_symbol",
        columns="phenotype_category",
        values=["min_p", "max_beta"],
    )
    organ_wide.columns = [
        f"plof_{col[1]}_{col[0]}" for col in organ_wide.columns
    ]
    organ_wide = organ_wide.reset_index()

    result = base.merge(organ_wide, on="gene_symbol", how="left")
    print(f"  Expanded pLoF: {len(result):,} genes")
    return result


def _ct_col(ct: str) -> str:
    """Sanitize a cell-type label into a valid column name prefix."""
    return "expr_ct_" + re.sub(r"[^a-z0-9]+", "_", ct.lower()).strip("_")


def load_expanded_expression() -> pd.DataFrame:
    """Cell-type and organ-level expression features for all DB genes → by gene_symbol.

    Uses ALL cell types present in the DB — no hardcoded filter.
    Column naming: expr_ct_<sanitized_cell_type> (mean) and ..._pct (% expressing).
    """
    conn = sqlite3.connect(DB_PATH)
    expr = pd.read_sql(
        "SELECT gene_symbol, cell_type, organ, mean_expression, pct_expressed "
        "FROM expression_summary",
        conn,
    )
    conn.close()

    expr["col"] = expr["cell_type"].map(_ct_col)

    # Mean expression pivot (one column per cell type)
    mean_pivot = (
        expr.pivot_table(
            index="gene_symbol", columns="col",
            values="mean_expression", aggfunc="mean",
        )
        .reset_index()
    )
    mean_pivot.columns.name = None

    # Percent expressing pivot
    expr["col_pct"] = expr["col"] + "_pct"
    pct_pivot = (
        expr.pivot_table(
            index="gene_symbol", columns="col_pct",
            values="pct_expressed", aggfunc="mean",
        )
        .reset_index()
    )
    pct_pivot.columns.name = None

    # Organ-level max expression (unchanged)
    organ_max = (
        expr
        .groupby(["gene_symbol", "organ"])["mean_expression"]
        .max()
        .unstack("organ")
        .add_prefix("expr_organ_")
        .reset_index()
    )

    result = mean_pivot.merge(pct_pivot, on="gene_symbol", how="outer")
    result = result.merge(organ_max, on="gene_symbol", how="left")
    n_expr_cols = len(result.columns) - 1
    print(f"  Expanded expression: {len(result):,} genes, {n_expr_cols} expression features")
    return result


def build_expanded_feature_matrix() -> pd.DataFrame:
    """Feature matrix for ALL DB genes (OT as master, expanded LOEUF/pLoF/expr)."""
    print("Loading Open Targets features (master, ~78K genes)...")
    ot = load_opentargets()

    print("Loading expanded LOEUF features...")
    loeuf = load_expanded_loeuf()

    print("Loading expanded pLoF features...")
    plof = load_expanded_plof()

    print("Loading expanded expression features...")
    expr = load_expanded_expression()

    # 1. Merge LOEUF by targetId (ENSG) — gives us gene_symbol for DB genes
    df = ot.merge(
        loeuf[["targetId", "gene_symbol", "loeuf_score"]],
        on="targetId",
        how="left",
    )

    # 2. Merge pLoF by gene_symbol (from LOEUF merge)
    df = df.merge(plof, on="gene_symbol", how="left", suffixes=("", "_plof"))

    # 3. Merge expression by gene_symbol
    df = df.merge(expr, on="gene_symbol", how="left", suffixes=("", "_expr"))

    print(f"\nExpanded feature matrix: {df.shape[0]:,} rows × {df.shape[1]} cols")
    n_with_sym = df["gene_symbol"].notna().sum()
    print(f"  Rows with DB gene_symbol: {n_with_sym:,}")
    return df


# ─── Main export ───────────────────────────────────────────────────────────────

def export_shap_all_genes() -> None:
    print("=" * 60)
    print("Export SHAP — ALL DB protein-coding genes")
    print("=" * 60)

    # ── Build feature matrix ──────────────────────────────────────
    df = build_expanded_feature_matrix()

    pval_cols    = [c for c in df.columns if "min_p" in c]
    df           = log_transform_pvals(df, pval_cols)
    feature_cols = get_feature_cols(df)
    print(f"\nFeature columns ({len(feature_cols)}): {feature_cols[:5]} ...")

    # ── Load ENSG → gene_symbol from DB (canonical mapping) ──────
    conn = sqlite3.connect(DB_PATH)
    db_genes = pd.read_sql(
        "SELECT gene_symbol, ensembl_id FROM gene_dosage_sensitivity WHERE ensembl_id != ''",
        conn,
    )
    conn.close()
    ensg_to_symbol: dict[str, str] = dict(zip(db_genes["ensembl_id"], db_genes["gene_symbol"]))
    print(f"\nDB gene map: {len(ensg_to_symbol):,} ENSG → symbol pairs")

    # ── Load Model B (drugged-only; applied to ALL genes) ─────────
    print("\nLoading XGBoost Model B...")
    model_b = xgb.XGBClassifier()
    model_b.load_model(OUT_DIR_DRUGGED / "model.json")

    # ── SHAP explainer ────────────────────────────────────────────
    print("Creating SHAP explainer...")
    explainer_b = shap.TreeExplainer(model_b)
    # expected_value may be a scalar or array (per-class) depending on SHAP version
    _ev = explainer_b.expected_value
    _ev_arr = np.array(_ev).flatten()
    base_b = float(_ev_arr[1] if len(_ev_arr) > 1 else _ev_arr[0])

    # ── Load OOF scores for drugged genes ─────────────────────────
    preds_b = pd.read_csv(OUT_DIR_DRUGGED / "predictions_oof.csv")
    score_map_b = dict(zip(preds_b["targetId"], preds_b["safety_score_oof"]))
    label_map_b = dict(zip(preds_b["targetId"], preds_b["safety_label"]))

    # ── Filter to DB genes ────────────────────────────────────────
    db_subset = df[df["targetId"].isin(ensg_to_symbol)].copy()
    print(f"\nTotal DB genes for SHAP: {len(db_subset):,}")

    # ── Batch-compute SHAP for all DB genes using Model B ─────────
    X_all = db_subset[feature_cols].values.astype(float)
    print(f"Computing SHAP for {len(db_subset):,} DB genes (Model B)...")
    _sv_raw = explainer_b.shap_values(X_all)
    # shap_values() may return a list [neg_class, pos_class] or a single 2D array
    sv_all = _sv_raw[1] if isinstance(_sv_raw, list) else _sv_raw

    result: dict = {}

    for idx, (_, row) in enumerate(db_subset.iterrows()):
        ensg     = row["targetId"]
        gene_sym = ensg_to_symbol[ensg]
        sv       = sv_all[idx]

        # Clinical metadata
        is_drugged     = bool(pd.notna(row["maxClinicalTrialPhase"]))
        clinical_phase = float(row["maxClinicalTrialPhase"]) if is_drugged else None
        hs_val         = row.get("hasSafetyEvent", float("nan"))
        if pd.isna(hs_val):
            has_safety_event = None
        elif float(hs_val) == -1.0:
            has_safety_event = True
        else:
            has_safety_event = False

        # Score: OOF for drugged (honest), predict_proba for undrugged (extrapolated)
        if is_drugged and ensg in score_map_b:
            raw_score = score_map_b[ensg]
            label     = label_map_b.get(ensg, "drugged_safe")
        else:
            raw_score = _sigmoid(base_b + float(sv.sum()))
            label     = "unlabeled"

        features = []
        for feat, sv_val in zip(feature_cols, sv):
            raw_val = row[feat]
            features.append({
                "name":          feat,
                "label":         READABLE_FEATURES.get(feat, feat),
                "group":         get_feature_group(feat),
                "shap_value":    round(float(sv_val), 6),
                "feature_value": None if pd.isna(raw_val) else round(float(raw_val), 6),
            })
        features.sort(key=lambda x: abs(x["shap_value"]), reverse=True)

        result[gene_sym] = {
            "ensembl_id":       ensg,
            "safety_score":     round(float(raw_score), 4),
            "base_value":       round(base_b, 6),
            "model":            "B",
            "is_drugged":       is_drugged,
            "clinical_phase":   clinical_phase,
            "has_safety_event": has_safety_event,
            "safety_label":     label,
            "features":         features,
        }

    n_drugged   = sum(1 for v in result.values() if v["is_drugged"])
    n_undrugged = sum(1 for v in result.values() if not v["is_drugged"])
    n_safety    = sum(1 for v in result.values() if v["has_safety_event"] is True)
    print(f"\n  Drugged genes (OOF score)      : {n_drugged:,}")
    print(f"  Undrugged genes (extrapolated) : {n_undrugged:,}")
    print(f"  Genes with safety event        : {n_safety:,}")
    print(f"  Total exported                 : {len(result):,}")

    # ── Save ──────────────────────────────────────────────────────
    out_path = OUT_DIR / "gene_shap_all.json"
    with open(out_path, "w") as f:
        json.dump(result, f, separators=(",", ":"))
    print(f"\nSaved → {out_path}  ({out_path.stat().st_size / 1024 / 1024:.1f} MB)")

    # Verify SHROOM3
    if "SHROOM3" in result:
        s = result["SHROOM3"]
        print(f"\nSHROOM3: score={s['safety_score']}, model={s['model']}, "
              f"ENSG={s['ensembl_id']}")


if __name__ == "__main__":
    export_shap_all_genes()
