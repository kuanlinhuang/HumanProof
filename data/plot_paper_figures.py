"""
Publication figures for HumanProof DR+PU safety model.

Outputs saved to figures/ (project root):
  fig1_roc_prc.png            — AUROC and AUPRC curves
  fig2_score_distribution.png — DR score distribution by safety label
  fig3_feature_importance.png — Top 20 features by mean |SHAP| (sampled)
  fig4_calibration.png        — Calibration curve on drugged genes (Step 3 OOF)

Usage:
    ~/miniconda3/bin/python data/plot_paper_figures.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_curve,
    precision_recall_curve,
    average_precision_score,
    roc_auc_score,
)
from sklearn.calibration import calibration_curve
import xgboost as xgb
import shap

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "data"))

from safety_model import get_feature_cols, log_transform_pvals, OUT_DIR   # noqa: E402
from export_shap_all_genes import build_expanded_feature_matrix            # noqa: E402
from export_gene_shap import get_feature_group, get_readable_label         # noqa: E402

OUT_DIR_DR = OUT_DIR / "dr"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ── Colour palette ─────────────────────────────────────────────────────────────
C = {
    "positive":     "#e53e3e",   # red
    "drugged_safe": "#805ad5",   # purple
    "unlabeled":    "#718096",   # gray
    "oof":          "#3182ce",   # blue  — Step 3 OOF
    "dr":           "#dd6b20",   # orange — Step 5 DR (in-sample)
    "random":       "#a0aec0",   # light gray
}
GROUP_COLORS = {"ot": "#3182ce", "genetics": "#38a169", "expression": "#d69e2e"}

plt.rcParams.update({
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "font.size":         11,
    "axes.titlesize":    13,
    "axes.labelsize":    12,
})


# ── Data loaders ───────────────────────────────────────────────────────────────

def load_predictions() -> pd.DataFrame:
    df = pd.read_csv(OUT_DIR_DR / "predictions.csv")
    df["true_label"] = (df["hasSafetyEvent"] == -1.0).astype(int)
    return df


def load_feature_matrix():
    print("  Building expanded feature matrix (may take ~1 min)…")
    df_raw = build_expanded_feature_matrix()
    pval_cols = [c for c in df_raw.columns if "min_p" in c]
    df_raw = log_transform_pvals(df_raw, pval_cols)
    df = df_raw[df_raw["gene_symbol"].notna()].copy().reset_index(drop=True)
    feature_cols = get_feature_cols(df)
    return df, feature_cols


# ── Figure 1: ROC + PRC ────────────────────────────────────────────────────────

def plot_roc_prc(preds: pd.DataFrame) -> None:
    drugged = preds[preds["is_drugged"] == 1].copy()
    Y = drugged["true_label"].values
    pos_rate = Y.mean()

    # Step 3: m_hat is the OOF outcome-model prediction for drugged genes
    # Step 5: safety_score_dr is the final DR regressor (in-sample on drugged)
    curves = [
        (drugged["m_hat"].values,          "Outcome model — OOF (Step 3)", C["oof"], "-"),
        (drugged["safety_score_dr"].values, "DR final model (Step 5)",      C["dr"],  "--"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Left: ROC ──
    ax = axes[0]
    for scores, label, color, ls in curves:
        fpr, tpr, _ = roc_curve(Y, scores)
        auroc = roc_auc_score(Y, scores)
        ax.plot(fpr, tpr, color=color, lw=2.2, ls=ls,
                label=f"{label}\nAUROC = {auroc:.3f}")
    ax.plot([0, 1], [0, 1], color=C["random"], lw=1, ls=":", label="Random classifier")
    ax.fill_between([0, 1], [0, 1], alpha=0.05, color="gray")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve", fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9, loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.text(0.98, 0.08, f"n drugged = {len(drugged):,}\nn positive = {int(Y.sum())}",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5, color="gray")

    # ── Right: PRC ──
    ax = axes[1]
    for scores, label, color, ls in curves:
        prec, rec, _ = precision_recall_curve(Y, scores)
        ap = average_precision_score(Y, scores)
        ax.plot(rec, prec, color=color, lw=2.2, ls=ls,
                label=f"{label}\nAUPRC = {ap:.3f}")
    ax.axhline(pos_rate, color=C["random"], lw=1, ls=":",
               label=f"Random (AUPRC ≈ {pos_rate:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curve", fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9, loc="upper right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.text(0.02, 0.08,
            "Step 3 OOF = honest held-out evaluation\nStep 5 = in-sample diagnostic",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8, color="gray", style="italic")

    fig.suptitle("HumanProof DR+PU Model — Discriminative Performance\n"
                 "(evaluated on drugged genes with known safety outcomes)",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    out = FIG_DIR / "fig1_roc_prc.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


# ── Figure 2: Score distribution ───────────────────────────────────────────────

def plot_score_distribution(preds: pd.DataFrame) -> None:
    label_order = ["positive", "drugged_safe", "unlabeled"]
    label_names = {
        "positive":     "Safety event\n(n={n})",
        "drugged_safe": "Drugged, no event\n(n={n})",
        "unlabeled":    "Undrugged / novel\n(n={n})",
    }

    data_by_label = {}
    for lbl in label_order:
        sub = preds.loc[preds["safety_label"] == lbl, "safety_score_dr"].dropna().values
        data_by_label[lbl] = sub

    fig, ax = plt.subplots(figsize=(9, 5.5))

    positions = [1, 2, 3]
    vp = ax.violinplot(
        [data_by_label[l] for l in label_order],
        positions=positions,
        widths=0.65,
        showmedians=True,
        showextrema=False,
    )
    for body, lbl in zip(vp["bodies"], label_order):
        body.set_facecolor(C[lbl])
        body.set_alpha(0.75)
    vp["cmedians"].set_color("black")
    vp["cmedians"].set_linewidth(1.5)

    # Overlay jittered strip for small groups
    for pos, lbl in zip(positions, label_order):
        vals = data_by_label[lbl]
        if len(vals) <= 600:
            jitter = np.random.default_rng(0).uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(pos + jitter, vals, s=6, color=C[lbl], alpha=0.55, linewidths=0)

    # Threshold lines
    ax.axhline(0.70, color="#c53030", lw=1.2, ls="--", alpha=0.7, label="High threshold (0.70)")
    ax.axhline(0.50, color="#c05621", lw=1.2, ls=":",  alpha=0.7, label="Moderate threshold (0.50)")

    tick_labels = [
        label_names[l].format(n=f"{len(data_by_label[l]):,}")
        for l in label_order
    ]
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_ylabel("DR Safety Score  m_DR(X)")
    ax.set_title("HumanProof Safety Score Distribution by Label",
                 fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.set_ylim(-0.05, 1.05)

    # Median annotation
    for pos, lbl in zip(positions, label_order):
        med = np.median(data_by_label[lbl])
        ax.text(pos, med + 0.03, f"med={med:.2f}", ha="center", va="bottom",
                fontsize=8, color="black", fontweight="bold")

    fig.tight_layout()
    out = FIG_DIR / "fig2_score_distribution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


# ── Figure 3: Feature importance (SHAP) ───────────────────────────────────────

def plot_feature_importance(df: pd.DataFrame, feature_cols: list[str],
                            n_sample: int = 3000) -> None:
    print(f"  Computing SHAP on {min(n_sample, len(df)):,} genes…")
    model = xgb.XGBRegressor()
    model.load_model(OUT_DIR_DR / "model_final.json")

    rng = np.random.default_rng(42)
    idx = rng.choice(len(df), size=min(n_sample, len(df)), replace=False)
    X_sample = df.iloc[idx][feature_cols].values.astype(float)

    explainer = shap.TreeExplainer(model)
    sv = explainer.shap_values(X_sample)
    mean_abs = np.abs(sv).mean(axis=0)

    importance = pd.DataFrame({"feature": feature_cols, "mean_abs_shap": mean_abs})
    importance = importance.sort_values("mean_abs_shap", ascending=False).head(20)
    importance["label"] = importance["feature"].map(get_readable_label)
    importance["group"] = importance["feature"].map(get_feature_group)
    colors = importance["group"].map(GROUP_COLORS).fillna("#a0aec0").tolist()

    fig, ax = plt.subplots(figsize=(9, 7.5))
    ax.barh(range(len(importance)), importance["mean_abs_shap"].values,
            color=colors, alpha=0.85, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(importance)))
    ax.set_yticklabels(importance["label"].values, fontsize=9.5)
    ax.invert_yaxis()
    ax.set_xlabel("Mean |SHAP value|  (DR model)")
    ax.set_title(f"Top 20 Features — DR+PU Model\n"
                 f"(SHAP on {min(n_sample, len(df)):,}-gene sample)",
                 fontweight="bold")

    legend_patches = [
        mpatches.Patch(color=GROUP_COLORS["ot"],         label="Open Targets"),
        mpatches.Patch(color=GROUP_COLORS["genetics"],   label="Genetics (LOEUF / pLoF)"),
        mpatches.Patch(color=GROUP_COLORS["expression"], label="Expression (CellxGene)"),
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc="lower right")

    fig.tight_layout()
    out = FIG_DIR / "fig3_feature_importance.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


# ── Figure 4: Calibration ──────────────────────────────────────────────────────

def plot_calibration(preds: pd.DataFrame) -> None:
    """Calibration of Step 3 OOF outcome model on drugged genes."""
    drugged = preds[preds["is_drugged"] == 1].dropna(subset=["m_hat"])
    Y = drugged["true_label"].values
    scores = drugged["m_hat"].values

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # ── Left: Calibration curve ──
    ax = axes[0]
    for n_bins, ls, alpha in [(10, "-", 0.9), (20, "--", 0.6)]:
        try:
            frac_pos, mean_pred = calibration_curve(Y, scores, n_bins=n_bins, strategy="quantile")
            ax.plot(mean_pred, frac_pos, color=C["oof"], lw=2, ls=ls, alpha=alpha,
                    label=f"{n_bins} bins (quantile)")
        except Exception:
            pass

    ax.plot([0, 1], [0, 1], color=C["random"], lw=1, ls=":", label="Perfect calibration")
    ax.set_xlabel("Mean Predicted Probability")
    ax.set_ylabel("Fraction of Positives (Observed)")
    ax.set_title("Calibration — Outcome Model OOF\n(drugged genes, Step 3)",
                 fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")

    # ── Right: Score histogram for drugged genes ──
    ax2 = axes[1]
    bins = np.linspace(0, 1, 31)
    pos_scores = scores[Y == 1]
    neg_scores = scores[Y == 0]
    ax2.hist(neg_scores, bins=bins, alpha=0.7, color=C["drugged_safe"],
             label=f"No safety event (n={len(neg_scores):,})", density=True)
    ax2.hist(pos_scores, bins=bins, alpha=0.8, color=C["positive"],
             label=f"Safety event (n={len(pos_scores):,})", density=True)
    ax2.set_xlabel("Outcome Model Score  m̂(X)  [OOF]")
    ax2.set_ylabel("Density")
    ax2.set_title("Score Separation — Drugged Genes Only\n(Step 3 OOF)",
                  fontweight="bold")
    ax2.legend(fontsize=9, framealpha=0.9)

    fig.suptitle("Model Calibration and Score Separation (Step 3 OOF Outcome Model)",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    out = FIG_DIR / "fig4_calibration.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("HumanProof — Generating Publication Figures")
    print(f"Output directory: {FIG_DIR}")
    print("=" * 60)

    print("\n[1/4] ROC + PRC curves")
    preds = load_predictions()
    plot_roc_prc(preds)

    print("\n[2/4] Score distribution")
    plot_score_distribution(preds)

    print("\n[3/4] Feature importance (SHAP)")
    df, feature_cols = load_feature_matrix()
    plot_feature_importance(df, feature_cols)

    print("\n[4/4] Calibration")
    plot_calibration(preds)

    print(f"\nDone. All figures → {FIG_DIR}/")


if __name__ == "__main__":
    main()
