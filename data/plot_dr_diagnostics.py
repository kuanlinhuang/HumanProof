"""
Diagnostic plots for the DR+PU safety prediction model.

Generates:
  data/safety_model_output/dr/figures/
    01_propensity_overlap.png     — P(drugged|X) distribution by group
    02_score_distribution.png     — DR score by safety label
    03_pseudo_outcome.png         — Ỹ distribution for drugged vs undrugged
    04_feature_importance.png     — Top 20 features by mean |SHAP| (from DR model)

Usage:
    ~/miniconda3/bin/python data/plot_dr_diagnostics.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import shap
import xgboost as xgb

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "data"))

from safety_model import get_feature_cols, log_transform_pvals, OUT_DIR
from export_shap_all_genes import build_expanded_feature_matrix

OUT_DIR_DR  = OUT_DIR / "dr"
FIG_DIR     = OUT_DIR_DR / "figures"
FIG_DIR.mkdir(exist_ok=True)

COLORS = {
    "positive":     "#e53e3e",   # red
    "drugged_safe": "#805ad5",   # purple
    "unlabeled":    "#718096",   # gray
    "drugged":      "#6b46c1",
    "undrugged":    "#4a5568",
}


def load_data():
    print("Loading feature matrix and DR predictions...")
    df_raw    = build_expanded_feature_matrix()
    pval_cols = [c for c in df_raw.columns if "min_p" in c]
    df_raw    = log_transform_pvals(df_raw, pval_cols)
    df        = df_raw[df_raw["gene_symbol"].notna()].copy().reset_index(drop=True)

    preds     = pd.read_csv(OUT_DIR_DR / "predictions.csv")
    # Rename to avoid conflict with OT's own safety_label column
    preds_sub = preds[["targetId", "safety_score_dr", "safety_label", "pi_hat", "m_hat", "Y_tilde"]].copy()
    preds_sub = preds_sub.rename(columns={"safety_label": "dr_label"})
    df        = df.merge(preds_sub, on="targetId", how="left")
    feature_cols = get_feature_cols(df)
    return df, feature_cols


def plot_propensity_overlap(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))

    d_pi = df.loc[df["is_drugged"] == 1, "pi_hat"].dropna()
    u_pi = df.loc[df["is_drugged"] == 0, "pi_hat"].dropna()

    bins = np.linspace(0, 1, 41)
    ax.hist(u_pi, bins=bins, alpha=0.7, color=COLORS["undrugged"],
            label=f"Undrugged (n={len(u_pi):,})", density=True)
    ax.hist(d_pi, bins=bins, alpha=0.7, color=COLORS["drugged"],
            label=f"Drugged (n={len(d_pi):,})", density=True)

    ax.set_xlabel("Propensity Score  P(S=1 | X)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("Propensity Overlap: Drugged vs Undrugged Genes", fontsize=13, fontweight="bold")
    ax.legend(framealpha=0.9)

    median_d = d_pi.median()
    median_u = u_pi.median()
    ax.axvline(median_d, color=COLORS["drugged"],   linestyle="--", lw=1.5,
               label=f"Drugged median: {median_d:.3f}")
    ax.axvline(median_u, color=COLORS["undrugged"], linestyle="--", lw=1.5,
               label=f"Undrugged median: {median_u:.3f}")
    ax.legend(framealpha=0.9, fontsize=9)

    fig.tight_layout()
    out = FIG_DIR / "01_propensity_overlap.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out}")


def plot_score_distribution(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5))

    bins = np.linspace(0, 1, 41)
    label_order = ["positive", "drugged_safe", "unlabeled"]
    label_names = {
        "positive":     "Safety event reported",
        "drugged_safe": "Drugged, no event",
        "unlabeled":    "Undrugged (novel)",
    }
    for lbl in label_order:
        sub = df.loc[df["dr_label"] == lbl, "safety_score_dr"].dropna()
        ax.hist(sub, bins=bins, alpha=0.75, color=COLORS[lbl],
                label=f"{label_names[lbl]} (n={len(sub):,})")

    ax.set_xlabel("DR Safety Score  m_DR(X)", fontsize=12)
    ax.set_ylabel("Gene Count", fontsize=12)
    ax.set_title("DR+PU Safety Score Distribution — All 17,745 DB Genes", fontsize=13, fontweight="bold")
    ax.legend(framealpha=0.9)
    ax.set_yscale("log")

    fig.tight_layout()
    out = FIG_DIR / "02_score_distribution.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out}")


def plot_pseudo_outcomes(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: Ỹ distribution
    ax = axes[0]
    d_yt = df.loc[df["is_drugged"] == 1, "Y_tilde"].dropna()
    u_yt = df.loc[df["is_drugged"] == 0, "Y_tilde"].dropna()

    bins_d = np.linspace(-1, 2, 61)
    ax.hist(u_yt.clip(-1, 2), bins=bins_d, alpha=0.7, color=COLORS["undrugged"],
            label=f"Undrugged Ỹ = m̂  (n={len(u_yt):,})")
    ax.hist(d_yt.clip(-1, 2), bins=bins_d, alpha=0.7, color=COLORS["drugged"],
            label=f"Drugged AIPW Ỹ  (n={len(d_yt):,})")
    ax.axvline(0, color="black", lw=0.8, linestyle=":")
    ax.axvline(1, color="black", lw=0.8, linestyle=":")
    ax.set_xlabel("Pseudo-outcome  Ỹ", fontsize=11)
    ax.set_ylabel("Gene Count (log)", fontsize=11)
    ax.set_title("AIPW Pseudo-outcomes", fontsize=12, fontweight="bold")
    ax.set_yscale("log")
    ax.legend(fontsize=9)

    # Right: DR score vs naive m̂ scatter (drugged only)
    ax2 = axes[1]
    d = df[df["is_drugged"] == 1].dropna(subset=["safety_score_dr", "m_hat"])
    colors_scatter = [
        COLORS["positive"] if lbl == "positive" else COLORS["drugged_safe"]
        for lbl in d["dr_label"]
    ]
    ax2.scatter(d["m_hat"], d["safety_score_dr"], c=colors_scatter,
                alpha=0.5, s=12, linewidths=0)
    ax2.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="y = x")
    ax2.set_xlabel("Naive m̂ (outcome model)", fontsize=11)
    ax2.set_ylabel("m_DR(X) (final DR score)", fontsize=11)
    ax2.set_title("DR Score vs Naive m̂ — Drugged Genes", fontsize=12, fontweight="bold")

    pos_patch  = mpatches.Patch(color=COLORS["positive"],     label="Safety event")
    safe_patch = mpatches.Patch(color=COLORS["drugged_safe"], label="No event")
    ax2.legend(handles=[pos_patch, safe_patch], fontsize=9)

    fig.tight_layout()
    out = FIG_DIR / "03_pseudo_outcome.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out}")


def plot_feature_importance(df: pd.DataFrame, feature_cols: list[str]) -> None:
    print("  Computing SHAP for feature importance plot...")
    model = xgb.XGBRegressor()
    model.load_model(OUT_DIR_DR / "model_final.json")

    X = df[feature_cols].values.astype(float)
    explainer   = shap.TreeExplainer(model)
    sv          = explainer.shap_values(X)
    mean_abs    = np.abs(sv).mean(axis=0)

    importance = pd.DataFrame({"feature": feature_cols, "mean_abs_shap": mean_abs})
    importance = importance.sort_values("mean_abs_shap", ascending=False).head(20)

    from export_gene_shap import get_feature_group, get_readable_label
    importance["label"] = importance["feature"].map(get_readable_label)
    importance["group"] = importance["feature"].map(get_feature_group)

    group_colors = {"ot": "#3182ce", "genetics": "#38a169", "expression": "#d69e2e"}
    colors = importance["group"].map(group_colors).tolist()

    fig, ax = plt.subplots(figsize=(9, 7))
    bars = ax.barh(range(len(importance)), importance["mean_abs_shap"].values,
                   color=colors, alpha=0.85)
    ax.set_yticks(range(len(importance)))
    ax.set_yticklabels(importance["label"].values, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Mean |SHAP value|  (DR model)", fontsize=11)
    ax.set_title("Top 20 Features — DR+PU Model", fontsize=13, fontweight="bold")

    legend_patches = [
        mpatches.Patch(color="#3182ce", label="Open Targets"),
        mpatches.Patch(color="#38a169", label="Genetics (LOEUF/pLoF)"),
        mpatches.Patch(color="#d69e2e", label="Expression"),
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc="lower right")

    fig.tight_layout()
    out = FIG_DIR / "04_feature_importance.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved → {out}")


def main() -> None:
    print("=" * 60)
    print("DR+PU Model — Diagnostic Plots")
    print("=" * 60)

    df, feature_cols = load_data()

    print("\nPlot 1: Propensity overlap")
    plot_propensity_overlap(df)

    print("Plot 2: Score distribution")
    plot_score_distribution(df)

    print("Plot 3: Pseudo-outcomes")
    plot_pseudo_outcomes(df)

    print("Plot 4: Feature importance")
    plot_feature_importance(df, feature_cols)

    print(f"\nAll figures → {FIG_DIR}")


if __name__ == "__main__":
    main()
