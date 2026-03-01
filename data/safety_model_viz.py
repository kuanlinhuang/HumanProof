"""
HumanProof safety model visualizations.

Produces four publication-quality figures saved to data/safety_model_output/figures/:

  fig1_performance.png      ROC + PR curves (Model B OOF)
  fig2_shap_importance.png  SHAP feature importance bar chart (Model B)
  fig3_humanproof_score.png HumanProof safety score for all 108 target genes
  fig4_score_distribution.png Score density: confirmed safety event vs no event

Usage:
    python data/safety_model_viz.py
"""

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve, precision_recall_curve

# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).parent.parent
OUT_BASE  = ROOT / "data" / "safety_model_output"
OUT_B     = OUT_BASE / "drugged_only"
FIG_DIR   = OUT_BASE / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ─── Palette ──────────────────────────────────────────────────────────────────
C_POS     = "#C0392B"   # confirmed safety event
C_SAFE    = "#27AE60"   # drugged, no safety event
C_UNK     = "#95A5A6"   # undrugged / unlabeled
C_ACCENT  = "#2980B9"   # blue accent for curves

READABLE_FEATURES = {
    "mouseKOScore":                       "Mouse KO phenotype score",
    "geneticConstraint":                  "Genetic constraint (OT)",
    "mouseOrthologMaxIdentityPercentage": "Mouse ortholog identity %",
    "hasSmallMoleculeBinder":             "Small molecule binder",
    "hasHighQualityChemicalProbes":       "High-quality chemical probes",
    "paralogMaxIdentityPercentage":       "Paralog identity %",
    "hasLigand":                          "Known ligand",
    "tissueSpecificity":                  "Tissue specificity",
    "hasPocket":                          "Binding pocket predicted",
    "tissueDistribution":                 "Tissue distribution",
    "isInMembrane":                       "Membrane protein",
    "isSecreted":                         "Secreted protein",
    "plof_respiratory_min_p":             "pLoF: Respiratory −log₁₀p",
    "expr_organ_Heart":                   "Heart expression (max)",
    "isCancerDriverGene":                 "Cancer driver gene",
    "plof_musculoskeletal_max_beta":      "pLoF: Musculoskeletal |β|",
    "expr_lymph_tcell":                   "T cell expression",
    "plof_n_phenotypes":                  "pLoF: N phenotypes",
    "loeuf_score":                        "LOEUF score",
    "expr_bm_erythroblasts":              "Erythroblast expression",
}

plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
})


# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

def load_data():
    preds_b  = pd.read_csv(OUT_B / "predictions_oof.csv")
    preds_a  = pd.read_csv(OUT_BASE / "predictions.csv")
    shap_b   = pd.read_csv(OUT_B / "shap_summary.csv")
    return preds_b, preds_a, shap_b


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: ROC + PR curves
# ─────────────────────────────────────────────────────────────────────────────

def fig1_performance(preds_b: pd.DataFrame) -> None:
    y    = preds_b["true_label"].values
    prob = preds_b["safety_score_oof"].values

    fpr, tpr, _  = roc_curve(y, prob)
    prec, rec, _ = precision_recall_curve(y, prob)
    auc  = roc_auc_score(y, prob)
    ap   = average_precision_score(y, prob)
    prev = y.mean()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle("Model B — Drugged-only  |  Out-of-fold performance (n = 1,564)",
                 fontsize=13, fontweight="bold", y=1.01)

    # ROC
    ax = axes[0]
    ax.plot(fpr, tpr, color=C_ACCENT, lw=2.2, label=f"OOF AUROC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="#BDC3C7", lw=1.2)
    ax.set_xlabel("False Positive Rate", fontsize=11)
    ax.set_ylabel("True Positive Rate", fontsize=11)
    ax.set_title("ROC Curve", fontweight="bold")
    ax.legend(loc="lower right", frameon=False, fontsize=11)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)

    # Precision-Recall
    ax = axes[1]
    ax.plot(rec, prec, color=C_POS, lw=2.2, label=f"OOF AUPRC = {ap:.3f}")
    ax.axhline(prev, ls="--", color="#BDC3C7", lw=1.2,
               label=f"Random baseline = {prev:.2f}")
    ax.set_xlabel("Recall", fontsize=11)
    ax.set_ylabel("Precision", fontsize=11)
    ax.set_title("Precision–Recall Curve", fontweight="bold")
    ax.legend(loc="upper right", frameon=False, fontsize=11)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.08)

    fig.tight_layout()
    path = FIG_DIR / "fig1_performance.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: SHAP feature importance
# ─────────────────────────────────────────────────────────────────────────────

def fig2_shap(shap_b: pd.DataFrame, n_top: int = 15) -> None:
    top = shap_b.head(n_top).copy()
    top["label"] = top["feature"].map(READABLE_FEATURES).fillna(top["feature"])

    # Colour bars by feature group
    def bar_color(feat):
        if feat.startswith("expr_"):
            return "#8E44AD"   # purple: expression
        if feat.startswith("plof_"):
            return "#E67E22"   # orange: pLoF
        if feat == "loeuf_score":
            return "#E67E22"
        return C_ACCENT        # blue: OT features

    colors = [bar_color(f) for f in top["feature"]]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    y_pos = np.arange(len(top))
    bars = ax.barh(y_pos, top["mean_abs_shap"], color=colors, height=0.65, edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top["label"], fontsize=10.5)
    ax.invert_yaxis()
    ax.set_xlabel("Mean |SHAP value| (impact on model output)", fontsize=11)
    ax.set_title("Feature Importance — Model B (Drugged-only)\nTop 15 features by SHAP",
                 fontweight="bold", fontsize=12)

    legend_patches = [
        mpatches.Patch(color=C_ACCENT,   label="Open Targets features"),
        mpatches.Patch(color="#E67E22",   label="pLoF / LOEUF"),
        mpatches.Patch(color="#8E44AD",   label="Cell-type expression"),
    ]
    ax.legend(handles=legend_patches, frameon=False, fontsize=10, loc="lower right")
    ax.set_xlim(0, top["mean_abs_shap"].max() * 1.18)

    # Value labels on bars
    for bar, val in zip(bars, top["mean_abs_shap"]):
        ax.text(val + 0.005, bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}", va="center", fontsize=9, color="#555")

    fig.tight_layout()
    path = FIG_DIR / "fig2_shap_importance.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: HumanProof safety score — all 108 target genes
# ─────────────────────────────────────────────────────────────────────────────

def fig3_humanproof_score(preds_b: pd.DataFrame, preds_a: pd.DataFrame) -> None:
    """Lollipop plot showing HumanProof safety score for all 108 genes.

    Drugged genes (73) use Model B OOF scores (reliable).
    Undrugged genes (35) use Model A scores (extrapolated, shown distinctly).
    """
    # ── Model B genes ────────────────────────────────────────────
    b108 = preds_b[preds_b["gene_symbol"].notna()].copy()
    b108 = b108.rename(columns={"safety_score_oof": "score"})
    b108["model"]  = "B (drugged-only)"
    b108["source"] = "Model B — OOF score"

    # ── Model A undrugged genes ───────────────────────────────────
    a_unk = preds_a[
        preds_a["gene_symbol"].notna() & (preds_a["safety_label"] == "unlabeled")
    ].copy()
    a_unk = a_unk.rename(columns={"safety_score": "score"})
    a_unk["model"]        = "A (extrapolated)"
    a_unk["safety_label"] = "unlabeled"
    a_unk["source"]       = "Model A — extrapolated"

    # Remaining 108 positives not in drugged (e.g. PTEN, VHL, BRCA1 with hasSafetyEvent but no trial)
    a_pos = preds_a[
        preds_a["gene_symbol"].notna() & (preds_a["safety_label"] == "positive") &
        ~preds_a["gene_symbol"].isin(b108["gene_symbol"])
    ].copy()
    a_pos = a_pos.rename(columns={"safety_score": "score"})
    a_pos["model"]  = "A (extrapolated)"
    a_pos["source"] = "Model A — extrapolated"

    # Combine and sort
    all_genes = pd.concat([b108, a_unk, a_pos], ignore_index=True)
    all_genes = all_genes.sort_values("score", ascending=True).reset_index(drop=True)

    # ── Plot ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 13))

    color_map  = {"positive": C_POS, "drugged_safe": C_SAFE, "unlabeled": C_UNK}
    marker_map = {"B (drugged-only)": "o", "A (extrapolated)": "D"}
    alpha_map  = {"B (drugged-only)": 1.0, "A (extrapolated)": 0.55}

    for i, row in all_genes.iterrows():
        col   = color_map.get(row["safety_label"], C_UNK)
        mkr   = marker_map[row["model"]]
        alpha = alpha_map[row["model"]]
        ms    = 8 if row["model"] == "B (drugged-only)" else 7

        ax.hlines(i, 0, row["score"], color="#E8E8E8", lw=1.0, zorder=1)
        ax.plot(row["score"], i, marker=mkr, color=col, ms=ms,
                alpha=alpha, zorder=3, markeredgewidth=0.4,
                markeredgecolor="white" if row["model"] == "B (drugged-only)" else col)

    ax.set_yticks(range(len(all_genes)))
    ax.set_yticklabels(all_genes["gene_symbol"], fontsize=8.5)
    ax.set_xlabel("HumanProof Safety Score", fontsize=12, labelpad=8)
    ax.set_xlim(-0.03, 1.0)
    ax.set_title(
        "HumanProof Safety Score — 108 Target Genes\n"
        "Model B (drugged) · Model A extrapolated (undrugged)",
        fontweight="bold", fontsize=12,
    )

    # Threshold line at 0.5
    ax.axvline(0.5, ls="--", color="#BDC3C7", lw=1.0, alpha=0.8, zorder=0)
    ax.text(0.505, len(all_genes) - 1, "0.5", fontsize=8.5, color="#999", va="top")

    # Legend
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_POS,
               markersize=9, label="Confirmed safety event"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_SAFE,
               markersize=9, label="Drugged — no safety event"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_UNK,
               markersize=9, label="Not yet drugged (unlabeled)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#555",
               markersize=9, label="● Model B OOF score (reliable)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#555",
               markersize=8, alpha=0.55, label="◆ Model A score (extrapolated)"),
    ]
    ax.legend(handles=legend_elements, frameon=True, framealpha=0.95,
              fontsize=9.5, loc="lower right", edgecolor="#DDD")

    fig.tight_layout()
    path = FIG_DIR / "fig3_humanproof_score.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Score distribution — positive vs drugged_safe
# ─────────────────────────────────────────────────────────────────────────────

def fig4_score_distribution(preds_b: pd.DataFrame) -> None:
    pos  = preds_b.loc[preds_b["safety_label"] == "positive",  "safety_score_oof"].values
    safe = preds_b.loc[preds_b["safety_label"] == "drugged_safe", "safety_score_oof"].values

    bins = np.linspace(0, 1, 26)

    fig, ax = plt.subplots(figsize=(8, 4.2))

    ax.hist(safe, bins=bins, color=C_SAFE, alpha=0.65, label=f"No safety event (n={len(safe):,})",
            edgecolor="white", linewidth=0.5)
    ax.hist(pos,  bins=bins, color=C_POS,  alpha=0.75, label=f"Confirmed safety event (n={len(pos):,})",
            edgecolor="white", linewidth=0.5)

    # Median lines
    ax.axvline(np.median(pos),  color=C_POS,  ls="--", lw=1.6, alpha=0.9,
               label=f"Median positive = {np.median(pos):.2f}")
    ax.axvline(np.median(safe), color=C_SAFE, ls="--", lw=1.6, alpha=0.9,
               label=f"Median safe = {np.median(safe):.2f}")

    ax.set_xlabel("HumanProof Safety Score (OOF)", fontsize=12)
    ax.set_ylabel("Number of targets", fontsize=12)
    ax.set_title(
        "Score Distribution — Drugged Targets (Model B)\n"
        "Confirmed safety event vs no recorded event",
        fontweight="bold", fontsize=12,
    )
    ax.legend(frameon=False, fontsize=10)
    ax.set_xlim(0, 1)

    fig.tight_layout()
    path = FIG_DIR / "fig4_score_distribution.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("HumanProof Safety Model — Visualizations")
    print("=" * 60)

    preds_b, preds_a, shap_b = load_data()

    print("\nFig 1: ROC + PR curves...")
    fig1_performance(preds_b)

    print("Fig 2: SHAP feature importance...")
    fig2_shap(shap_b)

    print("Fig 3: HumanProof score — 108 genes...")
    fig3_humanproof_score(preds_b, preds_a)

    print("Fig 4: Score distribution...")
    fig4_score_distribution(preds_b)

    print(f"\nAll figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()
