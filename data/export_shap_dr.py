"""
Export per-gene SHAP values for ALL DB genes using the DR+PU final model.

Uses the deployed m_DR(X) regressor from safety_model_dr.py — a single model
that produces valid, selection-bias-corrected safety scores for all 17,745
protein-coding genes, whether drugged or not.

Outputs:
  data/safety_model_output/dr/gene_shap_dr.json   (replaces gene_shap_all.json)

Each gene entry:
  safety_score      : DR safety score ∈ [0, 1] (higher = more concern)
  base_value        : SHAP base value (mean prediction in training set)
  model             : "DR"
  is_drugged        : bool
  clinical_phase    : float | null
  has_safety_event  : bool | null
  safety_label      : "positive" | "drugged_safe" | "unlabeled"
  features          : list of {name, label, group, shap_value, feature_value}
                      sorted by |shap_value| descending

Usage:
    ~/miniconda3/bin/python data/export_shap_dr.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import shap
import xgboost as xgb

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "data"))

from safety_model import (  # noqa: E402
    get_feature_cols,
    log_transform_pvals,
    OUT_DIR,
)
from export_shap_all_genes import build_expanded_feature_matrix  # noqa: E402
from export_gene_shap import get_feature_group, get_readable_label  # noqa: E402

OUT_DIR_DR  = OUT_DIR / "dr"
MODEL_PATH  = OUT_DIR_DR / "model_final.json"
PREDS_PATH  = OUT_DIR_DR / "predictions.csv"
OUT_JSON    = OUT_DIR_DR / "gene_shap_dr.json"


def export_shap_dr() -> None:
    print("=" * 60)
    print("Export SHAP — DR+PU model — ALL DB genes")
    print("=" * 60)

    # ── Feature matrix ────────────────────────────────────────────
    print("\nBuilding expanded feature matrix...")
    df_raw    = build_expanded_feature_matrix()
    pval_cols = [c for c in df_raw.columns if "min_p" in c]
    df_raw    = log_transform_pvals(df_raw, pval_cols)

    # Restrict to DB genes
    df = df_raw[df_raw["gene_symbol"].notna()].copy().reset_index(drop=True)
    print(f"DB genes: {len(df):,}")

    feature_cols = get_feature_cols(df)
    print(f"Feature columns: {len(feature_cols)}")

    # ── Load DR model + pre-computed scores ───────────────────────
    print("\nLoading DR final model...")
    model = xgb.XGBRegressor()
    model.load_model(MODEL_PATH)

    print("Loading DR predictions...")
    preds    = pd.read_csv(PREDS_PATH)
    score_map = dict(zip(preds["targetId"], preds["safety_score_dr"]))
    label_map = dict(zip(preds["targetId"], preds["safety_label"]))

    # ── SHAP explainer ────────────────────────────────────────────
    print("Creating SHAP TreeExplainer...")
    explainer  = shap.TreeExplainer(model)
    base_value = float(np.array(explainer.expected_value).flatten()[0])
    print(f"  Base value (mean DR score): {base_value:.4f}")

    # ── Batch SHAP for all DB genes ───────────────────────────────
    X_all = df[feature_cols].values.astype(float)
    print(f"Computing SHAP for {len(df):,} genes...")
    sv_raw = explainer.shap_values(X_all)
    # XGBRegressor returns a 2D array directly
    sv_all = sv_raw if not isinstance(sv_raw, list) else sv_raw[0]

    # ── Build output dict ─────────────────────────────────────────
    result: dict = {}
    ensg_to_sym = dict(zip(df["targetId"], df["gene_symbol"]))

    for idx, (_, row) in enumerate(df.iterrows()):
        ensg     = row["targetId"]
        gene_sym = ensg_to_sym[ensg]
        sv       = sv_all[idx]

        # Clinical metadata
        is_drugged       = bool(pd.notna(row["maxClinicalTrialPhase"]))
        clinical_phase   = float(row["maxClinicalTrialPhase"]) if is_drugged else None
        hs_val           = row.get("hasSafetyEvent", float("nan"))
        if pd.isna(hs_val):
            has_safety_event = None
        elif float(hs_val) == -1.0:
            has_safety_event = True
        else:
            has_safety_event = False

        # DR score (always from pre-computed predictions.csv)
        safety_score = float(score_map.get(ensg, float(np.clip(base_value + sv.sum(), 0, 1))))
        safety_label = label_map.get(ensg, "unlabeled")

        features = []
        for feat, sv_val in zip(feature_cols, sv):
            raw_val = row[feat]
            features.append({
                "name":          feat,
                "label":         get_readable_label(feat),
                "group":         get_feature_group(feat),
                "shap_value":    round(float(sv_val), 6),
                "feature_value": None if pd.isna(raw_val) else round(float(raw_val), 6),
            })
        features.sort(key=lambda x: abs(x["shap_value"]), reverse=True)

        result[gene_sym] = {
            "ensembl_id":       ensg,
            "safety_score":     round(safety_score, 4),
            "base_value":       round(base_value, 6),
            "model":            "DR",
            "is_drugged":       is_drugged,
            "clinical_phase":   clinical_phase,
            "has_safety_event": has_safety_event,
            "safety_label":     safety_label,
            "features":         features,
        }

    n_drugged   = sum(1 for v in result.values() if v["is_drugged"])
    n_undrugged = sum(1 for v in result.values() if not v["is_drugged"])
    n_safety    = sum(1 for v in result.values() if v["has_safety_event"] is True)
    print(f"\n  Drugged genes   : {n_drugged:,}")
    print(f"  Undrugged genes : {n_undrugged:,}")
    print(f"  Safety events   : {n_safety:,}")
    print(f"  Total exported  : {len(result):,}")

    # ── Save ──────────────────────────────────────────────────────
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, separators=(",", ":"))
    size_mb = OUT_JSON.stat().st_size / 1024 / 1024
    print(f"\nSaved → {OUT_JSON}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    export_shap_dr()
