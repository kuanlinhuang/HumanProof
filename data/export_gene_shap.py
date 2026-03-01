"""
Export per-gene SHAP values for the 108 HumanProof target genes.

Produces: data/safety_model_output/gene_shap_108.json

Each gene entry contains:
  safety_score  : HumanProof safety score (0–1; higher = more safety concern)
  base_value    : SHAP base (expected model output in log-odds space)
  model         : "B" (drugged-only OOF, recommended) or "A" (all-labeled)
  features      : list of {name, label, group, shap_value, feature_value}
                  sorted by |shap_value| descending

Usage:
    ~/miniconda3/bin/python data/export_gene_shap.py
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
    build_feature_matrix,
    get_feature_cols,
    log_transform_pvals,
    GENES,
    OUT_DIR,
    OUT_DIR_DRUGGED,
)

# ─── Human-readable feature labels ────────────────────────────────────────────

READABLE_FEATURES: dict[str, str] = {
    # Open Targets
    "isInMembrane":                        "Membrane protein",
    "isSecreted":                          "Secreted protein",
    "hasPocket":                           "Binding pocket predicted",
    "hasLigand":                           "Known ligand",
    "hasSmallMoleculeBinder":              "Small molecule binder",
    "geneticConstraint":                   "Genetic constraint (OT)",
    "paralogMaxIdentityPercentage":        "Paralog identity %",
    "mouseOrthologMaxIdentityPercentage":  "Mouse ortholog identity %",
    "isCancerDriverGene":                  "Cancer driver gene",
    "hasTEP":                              "Target Enabling Package",
    "mouseKOScore":                        "Mouse KO phenotype score",
    "hasHighQualityChemicalProbes":        "High-quality chemical probes",
    "tissueSpecificity":                   "Tissue specificity",
    "tissueDistribution":                  "Tissue distribution",
    # LOEUF
    "loeuf_score":                         "LOEUF constraint score",
    # pLoF global
    "plof_min_pval":                       "pLoF: min −log₁₀(p)",
    "plof_n_sig":                          "pLoF: genome-wide sig. hits",
    "plof_max_abs_beta":                   "pLoF: max |β| (all pheno)",
    "plof_n_phenotypes":                   "pLoF: N phenotypes",
    # pLoF organs
    "plof_cardiovascular_min_p":           "pLoF: Cardiovascular −log₁₀p",
    "plof_cardiovascular_max_beta":        "pLoF: Cardiovascular |β|",
    "plof_hepatic_min_p":                  "pLoF: Hepatic −log₁₀p",
    "plof_hepatic_max_beta":               "pLoF: Hepatic |β|",
    "plof_neurological_min_p":             "pLoF: Neurological −log₁₀p",
    "plof_neurological_max_beta":          "pLoF: Neurological |β|",
    "plof_renal_min_p":                    "pLoF: Renal −log₁₀p",
    "plof_renal_max_beta":                 "pLoF: Renal |β|",
    "plof_respiratory_min_p":              "pLoF: Respiratory −log₁₀p",
    "plof_respiratory_max_beta":           "pLoF: Respiratory |β|",
    "plof_hematologic_min_p":              "pLoF: Hematologic −log₁₀p",
    "plof_hematologic_max_beta":           "pLoF: Hematologic |β|",
    "plof_musculoskeletal_min_p":          "pLoF: Musculoskeletal −log₁₀p",
    "plof_musculoskeletal_max_beta":       "pLoF: Musculoskeletal |β|",
    # Organ-level max expression
    "expr_organ_Brain":                    "Brain expression (max)",
    "expr_organ_Heart":                    "Heart expression (max)",
    "expr_organ_Liver":                    "Liver expression (max)",
    "expr_organ_Kidney":                   "Kidney expression (max)",
    "expr_organ_Lung":                     "Lung expression (max)",
    "expr_organ_Muscle":                   "Muscle expression (max)",
    "expr_organ_Bone marrow":              "Bone marrow expression (max)",
    "expr_organ_Spleen":                   "Spleen expression (max)",
    "expr_organ_Pancreas":                 "Pancreas expression (max)",
    "expr_organ_Skin":                     "Skin expression (max)",
    "expr_organ_Intestine":                "Intestine expression (max)",
    "expr_organ_Adrenal":                  "Adrenal expression (max)",
    "expr_organ_Breast":                   "Breast expression (max)",
    "expr_organ_Lymph node":               "Lymph node expression (max)",
}


def get_feature_group(feat: str) -> str:
    if feat.startswith("expr_"):
        return "expression"
    if feat.startswith("plof_") or feat == "loeuf_score":
        return "genetics"
    return "ot"


def get_readable_label(feat: str) -> str:
    """Return a human-readable label for any feature name.

    For static features, looks up READABLE_FEATURES.
    For dynamic cell-type features (expr_ct_*) and pLoF organ features
    (plof_*_min_p / plof_*_max_beta), generates a label from the column name.
    """
    if feat in READABLE_FEATURES:
        return READABLE_FEATURES[feat]
    if feat.startswith("expr_ct_"):
        base = feat[len("expr_ct_"):]
        if base.endswith("_pct"):
            ct = base[:-4].replace("_", " ").title()
            return f"{ct} (% expressing)"
        ct = base.replace("_", " ").title()
        return f"{ct} (expression)"
    if feat.startswith("plof_") and feat.endswith("_min_p"):
        organ = feat[len("plof_"):-len("_min_p")].replace("_", " ").title()
        return f"pLoF: {organ} \u2212log\u2081\u2080p"
    if feat.startswith("plof_") and feat.endswith("_max_beta"):
        organ = feat[len("plof_"):-len("_max_beta")].replace("_", " ").title()
        return f"pLoF: {organ} |\u03b2|"
    return feat


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + np.exp(-x))


def export_gene_shap() -> None:
    print("=" * 60)
    print("Export per-gene SHAP values — 108 HumanProof target genes")
    print("=" * 60)

    # ── Load full feature matrix ──────────────────────────────────
    print("\nBuilding feature matrix (this may take ~30 s)...")
    df = build_feature_matrix()

    pval_cols    = [c for c in df.columns if "min_p" in c]
    df           = log_transform_pvals(df, pval_cols)
    feature_cols = get_feature_cols(df)
    print(f"  Feature columns: {len(feature_cols)}")

    # ── Partition into drugged / undrugged ────────────────────────
    drugged   = df[df["maxClinicalTrialPhase"].notna()].reset_index(drop=True)
    undrugged = df[df["maxClinicalTrialPhase"].isna()].reset_index(drop=True)

    # ── Load models ───────────────────────────────────────────────
    print("\nLoading XGBoost models...")
    model_b = xgb.XGBClassifier()
    model_b.load_model(OUT_DIR_DRUGGED / "model.json")

    model_a = xgb.XGBClassifier()
    model_a.load_model(OUT_DIR / "model.json")

    # ── SHAP explainers ───────────────────────────────────────────
    print("Creating SHAP explainers...")
    explainer_b = shap.TreeExplainer(model_b)
    explainer_a = shap.TreeExplainer(model_a)

    # ── Load OOF predictions for safety scores ────────────────────
    preds_b = pd.read_csv(OUT_DIR_DRUGGED / "predictions_oof.csv")
    preds_a = pd.read_csv(OUT_DIR / "predictions.csv")

    score_map_b = dict(zip(preds_b["targetId"], preds_b["safety_score_oof"]))
    score_map_a = dict(zip(preds_a["targetId"], preds_a["safety_score"]))
    label_map_b = dict(zip(preds_b["targetId"], preds_b["safety_label"]))
    label_map_a = dict(zip(preds_a["targetId"], preds_a["safety_label"]))

    # ── Compute SHAP for each of our 108 genes ────────────────────
    result: dict = {}
    n_model_b = 0
    n_model_a = 0

    for gene_symbol, ensg in GENES.items():
        mask_b = drugged["targetId"] == ensg
        mask_a = undrugged["targetId"] == ensg

        if mask_b.any():
            # Preferred: Model B (drugged-only, honest evaluation)
            row      = drugged[mask_b].iloc[0]
            X        = row[feature_cols].values.reshape(1, -1).astype(float)
            sv       = explainer_b.shap_values(X)[0]
            base_val = float(explainer_b.expected_value)
            score    = float(score_map_b.get(ensg, _sigmoid(base_val + sv.sum())))
            label    = label_map_b.get(ensg, "drugged_safe")
            model    = "B"
            n_model_b += 1

        elif mask_a.any():
            # Fallback: Model A (undrugged genes)
            row      = undrugged[mask_a].iloc[0]
            X        = row[feature_cols].values.reshape(1, -1).astype(float)
            sv       = explainer_a.shap_values(X)[0]
            base_val = float(explainer_a.expected_value)
            score    = float(score_map_a.get(ensg, _sigmoid(base_val + sv.sum())))
            label    = label_map_a.get(ensg, "unlabeled")
            model    = "A"
            n_model_a += 1

        else:
            print(f"  WARNING: {gene_symbol} ({ensg}) not found in feature matrix, skipping")
            continue

        # Build per-feature list, sorted by |shap_value|
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

        result[gene_symbol] = {
            "ensembl_id":   ensg,
            "safety_score": round(score, 4),
            "base_value":   round(base_val, 6),
            "model":        model,
            "safety_label": label,
            "features":     features,
        }

    print(f"\n  Genes using Model B (drugged-only): {n_model_b}")
    print(f"  Genes using Model A (extrapolated) : {n_model_a}")
    print(f"  Total exported                     : {len(result)}")

    out_path = OUT_DIR / "gene_shap_108.json"
    with open(out_path, "w") as f:
        json.dump(result, f, separators=(",", ":"))

    print(f"\nSaved → {out_path}  ({out_path.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    export_gene_shap()
