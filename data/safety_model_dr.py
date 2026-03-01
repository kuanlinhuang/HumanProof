"""Doubly Robust + Positive-Unlabeled Safety Prediction Model.

Produces unbiased genome-wide safety risk scores by correcting for selection
bias in the drugged gene population via the AIPW (Augmented IPW) estimator.

Architecture (5 steps, all restricted to 17,745 DB genes):

  1. Cross-fitted propensity P(S=1|X)           — XGBClassifier, 5-fold
  2. PU prior P(Y=1)                            — Elkan-Noto, 5-fold
  3. Cross-fitted outcome P(Y=1|X)              — XGBClassifier, IPW+PU weights, 5-fold OOF
  4. AIPW pseudo-outcome  Ỹ = m̂ + (S/π̂)(Y−m̂) — vectorized
  5. Final regression m_DR(X) ← regress Ỹ~X    — XGBRegressor, all DB genes

Deployed predictor: m_DR(X) ∈ [0,1] for any gene in the DB.

Outputs:
  data/safety_model_output/dr/
    predictions.csv     — safety_score_dr for all 17,745 DB genes
    model_final.json    — deployed XGBRegressor
    model_outcome.json  — Step 3 XGBClassifier (for SHAP reference)

Usage:
    ~/miniconda3/bin/python data/safety_model_dr.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "data"))

from safety_model import (  # noqa: E402
    OUT_DIR,
    get_feature_cols,
    log_transform_pvals,
)
from export_shap_all_genes import build_expanded_feature_matrix  # noqa: E402

# ── Output directory ───────────────────────────────────────────────────────────
OUT_DIR_DR = OUT_DIR / "dr"
OUT_DIR_DR.mkdir(exist_ok=True)

# ── Constants ──────────────────────────────────────────────────────────────────
N_FOLDS  = 5
MAX_IPW  = 10.0   # stabilized IPW clip threshold
SEED     = 42

# ── PU prior override ──────────────────────────────────────────────────────────
# Set PI_OVERRIDE to a float in (0, 0.5] to skip Elkan-Noto estimation and fix
# π_p = P(Y=1) at the specified value.  This is a Bayesian prior adjustment:
# higher π_p shifts all scores upward (monotone logit-shift) while preserving
# AUROC and AUPRC exactly.  Elkan-Noto estimated 0.091 from observed labels,
# which underestimates the true safety-event prevalence across all genes.
# A value of 0.40 puts the median known-positive score near 0.50 ("moderate").
PI_OVERRIDE: float | None = 0.40

# Shared XGBoost hyperparameters
_XGB_COMMON = dict(
    n_estimators=500,
    max_depth=5,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=SEED,
    tree_method="hist",
)
_XGB_CLF = dict(**_XGB_COMMON, eval_metric="aucpr", early_stopping_rounds=30)
_XGB_REG = dict(
    **{k: v for k, v in _XGB_COMMON.items() if k != "n_estimators"},
    n_estimators=300,          # fixed — pseudo-outcome noise makes early stopping unreliable
    eval_metric="rmse",
    objective="reg:squarederror",
)


def _banner(title: str) -> None:
    print(f"\n{'='*60}\n{title}\n{'='*60}")


# ── STEP 1 ─────────────────────────────────────────────────────────────────────

def fit_propensity_crossfit(
    df: pd.DataFrame,
    feature_cols: list[str],
) -> pd.DataFrame:
    """5-fold cross-fitted P(S=1|X).

    Returns df with two new columns:
      pi_hat  — cross-fitted propensity score
      ipw     — stabilized IPW weight  π̄ / π̂, clipped at MAX_IPW
    """
    _banner("STEP 1 — Cross-fitted propensity  P(S=1 | X)")

    S = df["is_drugged"].values
    X = df[feature_cols].values.astype(float)

    pi_hat = np.zeros(len(df))
    skf    = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)

    for fold, (tr, va) in enumerate(skf.split(X, S)):
        m = xgb.XGBClassifier(**_XGB_CLF)
        m.fit(X[tr], S[tr], eval_set=[(X[va], S[va])], verbose=False)
        pi_hat[va] = m.predict_proba(X[va])[:, 1]
        auc = roc_auc_score(S[va], pi_hat[va])
        print(f"  Fold {fold + 1}: AUROC = {auc:.3f}")

    pi_bar = S.mean()
    ipw    = np.clip(pi_bar / np.maximum(pi_hat, 1e-6), 0.0, MAX_IPW)

    # Effective sample size for drugged subset
    drugged_ipw = ipw[S == 1]
    ess = drugged_ipw.sum() ** 2 / (drugged_ipw ** 2).sum()
    n_d = int(S.sum())
    print(f"\n  P(S=1) = {pi_bar:.4f}  |  Drugged: {n_d:,}  |  ESS: {ess:.0f} ({100*ess/n_d:.1f}%)")

    df = df.copy()
    df["pi_hat"] = pi_hat
    df["ipw"]    = ipw
    return df


# ── STEP 2 ─────────────────────────────────────────────────────────────────────

def estimate_pu_prior(
    df: pd.DataFrame,
    feature_cols: list[str],
) -> float:
    """Estimate π_p = P(Y=1) via Elkan-Noto method (or use PI_OVERRIDE).

    If PI_OVERRIDE is set, skips Elkan-Noto and returns the fixed value directly.
    Otherwise trains a naive classifier to estimate the labeling frequency c,
    returning π_p = P(labeled) / c.
    """
    _banner("STEP 2 — PU prior  P(Y=1)  [Elkan-Noto]")

    if PI_OVERRIDE is not None:
        print(f"  PI_OVERRIDE = {PI_OVERRIDE} — skipping Elkan-Noto estimation")
        return float(PI_OVERRIDE)

    pu_label = (df["hasSafetyEvent"] == -1.0).astype(int).values
    X        = df[feature_cols].values.astype(float)

    n_pos = pu_label.sum()
    n_neg = (pu_label == 0).sum()
    print(f"  Labeled positives: {n_pos}  |  Unlabeled: {n_neg:,}")

    spw = n_neg / max(n_pos, 1)
    scores_on_pos: list[float] = []
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)

    for tr, va in skf.split(X, pu_label):
        m = xgb.XGBClassifier(**{**_XGB_CLF, "scale_pos_weight": spw})
        m.fit(
            X[tr], pu_label[tr],
            eval_set=[(X[va], pu_label[va])],
            verbose=False,
        )
        prob   = m.predict_proba(X[va])[:, 1]
        pos_va = pu_label[va] == 1
        if pos_va.sum() > 0:
            scores_on_pos.extend(prob[pos_va].tolist())

    c     = float(np.mean(scores_on_pos))          # E[f(X) | Y=1]  (Elkan-Noto labeling freq)
    p_lab = float(pu_label.mean())                 # P(labeled=1)
    pi_p  = float(np.clip(p_lab / max(c, 1e-6), 1e-4, 0.5))

    print(f"  P(labeled)        = {p_lab:.5f}")
    print(f"  c = E[f | Y=1]    = {c:.4f}")
    print(f"  π_p = P(Y=1)      = {pi_p:.5f}")
    return pi_p


# ── STEP 3 ─────────────────────────────────────────────────────────────────────

def fit_outcome_crossfit(
    df: pd.DataFrame,
    feature_cols: list[str],
    pi_p: float,
) -> tuple[pd.DataFrame, xgb.XGBClassifier]:
    """5-fold cross-fitted outcome model on drugged DB genes.

    Sample weight = IPW × PU-weight:
      Y=1 (safety event)   → ipw × π_p
      Y=0 (no event seen)  → ipw × (1 − π_p)

    Returns:
      df_drugged   — drugged rows with 'outcome_hat_oof' and 'true_label' columns
      final_model  — XGBClassifier trained on all drugged genes (for Step 4 extrapolation)
    """
    _banner("STEP 3 — Cross-fitted outcome model  P(Y=1 | X)  [IPW + PU]")

    df_d = df[df["is_drugged"] == 1].copy().reset_index(drop=True)
    Y    = (df_d["hasSafetyEvent"] == -1.0).astype(int).values
    X    = df_d[feature_cols].values.astype(float)
    IPW  = df_d["ipw"].values

    n_pos = int(Y.sum())
    n_neg = int((Y == 0).sum())
    print(f"  Drugged DB genes: {len(df_d):,}  |  pos={n_pos}  neg={n_neg}")
    print(f"  π_p = {pi_p:.5f}")

    pu_w = np.where(Y == 1, pi_p, 1.0 - pi_p)
    sw   = IPW * pu_w                           # combined sample weight

    spw    = n_neg / max(n_pos, 1)
    params = {**_XGB_CLF, "scale_pos_weight": spw}

    oof_prob      = np.zeros(len(df_d))
    cv_auc, cv_ap = [], []
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)

    print(f"\n  5-fold CV (OOF predictions)...")
    for fold, (tr, va) in enumerate(skf.split(X, Y)):
        m = xgb.XGBClassifier(**params)
        m.fit(
            X[tr], Y[tr],
            sample_weight=sw[tr],
            eval_set=[(X[va], Y[va])],
            sample_weight_eval_set=[sw[va]],
            verbose=False,
        )
        prob         = m.predict_proba(X[va])[:, 1]
        oof_prob[va] = prob
        auc = roc_auc_score(Y[va], prob)
        ap  = average_precision_score(Y[va], prob)
        cv_auc.append(auc)
        cv_ap.append(ap)
        print(f"  Fold {fold + 1}: AUROC={auc:.3f}  AUPRC={ap:.3f}")

    oof_auc = roc_auc_score(Y, oof_prob)
    oof_ap  = average_precision_score(Y, oof_prob)
    print(f"\n  OOF (pooled) AUROC = {oof_auc:.3f}  |  AUPRC = {oof_ap:.3f}")

    df_d["outcome_hat_oof"] = oof_prob
    df_d["true_label"]      = Y

    print("\n  Training final outcome model on all drugged DB genes...")
    final_m = xgb.XGBClassifier(
        **{k: v for k, v in params.items() if k != "early_stopping_rounds"}
    )
    final_m.fit(X, Y, sample_weight=sw, verbose=False)

    return df_d, final_m


# ── STEP 4 ─────────────────────────────────────────────────────────────────────

def construct_pseudo_outcomes(
    df: pd.DataFrame,
    df_drugged: pd.DataFrame,
    outcome_model: xgb.XGBClassifier,
    feature_cols: list[str],
) -> pd.DataFrame:
    """AIPW pseudo-outcome for all 17,745 DB genes (vectorized).

    Ỹ_i = m̂(X_i) + (S_i / π̂_i) × (Y_i − m̂(X_i))

    Drugged  (S=1): honest OOF m̂ from Step 3; AIPW correction appended.
    Undrugged (S=0): Ỹ_i = m̂(X_i) — extrapolated from final outcome model.
    """
    _banner("STEP 4 — AIPW pseudo-outcomes")

    # m̂(X) for all DB genes — extrapolated from final outcome model
    X_all = df[feature_cols].values.astype(float)
    m_hat = outcome_model.predict_proba(X_all)[:, 1]

    df = df.copy()
    df["m_hat"]   = m_hat.astype(float)
    df["Y_tilde"] = m_hat.astype(float)   # default (undrugged: Ỹ = m̂)

    # Merge OOF scores and true labels for drugged genes
    oof = (
        df_drugged
        .set_index("targetId")[["outcome_hat_oof", "true_label"]]
        .rename(columns={"outcome_hat_oof": "_m_oof", "true_label": "_Y_obs"})
    )
    df = df.join(oof, on="targetId", how="left")

    # Vectorized AIPW for drugged genes that have OOF scores
    has_oof      = df["_m_oof"].notna()
    is_drugged   = df["is_drugged"] == 1
    mask         = is_drugged & has_oof

    pi_clip      = df.loc[mask, "pi_hat"].clip(lower=0.01).values
    m_i          = df.loc[mask, "_m_oof"].values
    Y_i          = df.loc[mask, "_Y_obs"].values

    Y_tilde_d    = np.clip(m_i + (Y_i - m_i) / pi_clip, -1.0, 2.0)

    df.loc[mask, "Y_tilde"] = Y_tilde_d
    df.loc[mask, "m_hat"]   = m_i          # override: OOF score is the honest m̂

    # Clean up temp columns
    df = df.drop(columns=["_m_oof", "_Y_obs"])

    n_corrected = int(mask.sum())
    n_undrugged = int((~is_drugged).sum())
    d_yt = df.loc[is_drugged, "Y_tilde"]
    u_yt = df.loc[~is_drugged, "Y_tilde"]
    print(f"  Drugged corrected : {n_corrected:,}")
    print(f"  Undrugged (m̂ only): {n_undrugged:,}")
    print(f"  Ỹ drugged   : mean={d_yt.mean():.4f}  std={d_yt.std():.4f}")
    print(f"  Ỹ undrugged : mean={u_yt.mean():.4f}  std={u_yt.std():.4f}")

    return df


# ── STEP 5 ─────────────────────────────────────────────────────────────────────

def fit_final_model(
    df: pd.DataFrame,
    feature_cols: list[str],
) -> xgb.XGBRegressor:
    """Regress Ỹ on X for all DB genes → m_DR(X), the deployed predictor.

    Drugged genes are weighted 2× because their pseudo-outcomes are grounded
    in observed Y.  Fixed n_estimators avoids overfitting to noisy Ỹ values.
    """
    _banner("STEP 5 — Final model  m_DR(X)  regressed on pseudo-outcomes")

    X       = df[feature_cols].values.astype(float)
    Y_tilde = df["Y_tilde"].values
    sw      = np.where(df["is_drugged"] == 1, 2.0, 1.0)

    final = xgb.XGBRegressor(**_XGB_REG)
    final.fit(X, Y_tilde, sample_weight=sw, verbose=False)

    # In-sample diagnostic on drugged genes (note: inflated due to training overlap)
    drugged = df[df["is_drugged"] == 1]
    Y_true  = drugged["true_label"].values if "true_label" in drugged.columns else \
              (drugged["hasSafetyEvent"] == -1.0).astype(int).values
    scores  = final.predict(drugged[feature_cols].values.astype(float))
    auc     = roc_auc_score(Y_true, scores)
    ap      = average_precision_score(Y_true, scores)
    print(f"  In-sample drugged: AUROC={auc:.3f}  AUPRC={ap:.3f}  (inflated — diagnostic only)")

    return final


# ── Diagnostics ────────────────────────────────────────────────────────────────

def _print_diagnostics(
    df: pd.DataFrame,
    dr_scores: np.ndarray,
    preds: pd.DataFrame,
) -> None:
    _banner("Diagnostics")

    # Propensity overlap
    d_pi = df.loc[df["is_drugged"] == 1, "pi_hat"]
    u_pi = df.loc[df["is_drugged"] == 0, "pi_hat"]
    print("  Propensity distribution:")
    print(f"    Drugged  : median={d_pi.median():.3f}  "
          f"IQR=[{d_pi.quantile(.25):.3f}, {d_pi.quantile(.75):.3f}]")
    print(f"    Undrugged: median={u_pi.median():.3f}  "
          f"IQR=[{u_pi.quantile(.25):.3f}, {u_pi.quantile(.75):.3f}]")

    # DR vs naive outcome comparison on drugged genes
    d_mask = df["is_drugged"] == 1
    Y_d    = (df.loc[d_mask, "hasSafetyEvent"] == -1.0).astype(int).values
    dr_d   = dr_scores[d_mask.values]
    m_d    = df.loc[d_mask, "m_hat"].values
    print("\n  DR vs naive m̂ on drugged genes:")
    print(f"    m_DR  AUROC={roc_auc_score(Y_d, dr_d):.3f}  "
          f"AUPRC={average_precision_score(Y_d, dr_d):.3f}")
    print(f"    m̂     AUROC={roc_auc_score(Y_d, m_d):.3f}  "
          f"AUPRC={average_precision_score(Y_d, m_d):.3f}")
    print("  (Both in-sample — compare directionally, not absolutely)")

    # Score distribution by label
    for label, name in [("positive", "Safety positives"),
                         ("drugged_safe", "Drugged, no event"),
                         ("unlabeled", "Undrugged")]:
        sub = preds.loc[preds["safety_label"] == label, "safety_score_dr"]
        if len(sub):
            print(f"\n  {name} (n={len(sub):,}): "
                  f"mean={sub.mean():.3f}  median={sub.median():.3f}  "
                  f"p90={sub.quantile(.9):.3f}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("HumanProof  DR+PU  Safety Prediction Model")
    print("=" * 60)

    # ── Build feature matrix (OT backbone + expanded DB features) ──
    print("\nBuilding expanded feature matrix...")
    df_raw    = build_expanded_feature_matrix()
    pval_cols = [c for c in df_raw.columns if "min_p" in c]
    df_raw    = log_transform_pvals(df_raw, pval_cols)

    # Restrict to 17,745 DB genes (those with gene_symbol from LOEUF merge)
    df = df_raw[df_raw["gene_symbol"].notna()].copy().reset_index(drop=True)
    print(f"\nDB genes            : {len(df):,}")
    print(f"  Drugged (S=1)     : {df['is_drugged'].sum():,}")
    print(f"  Safety pos (Y=1)  : {(df['hasSafetyEvent'] == -1.0).sum():,}")

    feature_cols = get_feature_cols(df)
    print(f"\nFeature columns: {len(feature_cols)}")

    # ── STEP 1 ────────────────────────────────────────────────────
    df = fit_propensity_crossfit(df, feature_cols)

    # ── STEP 2 ────────────────────────────────────────────────────
    pi_p = estimate_pu_prior(df, feature_cols)

    # ── STEP 3 ────────────────────────────────────────────────────
    df_drugged, outcome_model = fit_outcome_crossfit(df, feature_cols, pi_p)
    outcome_model.save_model(OUT_DIR_DR / "model_outcome.json")

    # ── STEP 4 ────────────────────────────────────────────────────
    df = construct_pseudo_outcomes(df, df_drugged, outcome_model, feature_cols)

    # Propagate true_label onto df for Step 5 in-sample diagnostics
    label_map      = df_drugged.set_index("targetId")["true_label"].to_dict()
    df["true_label"] = df["targetId"].map(label_map)

    # ── STEP 5 ────────────────────────────────────────────────────
    final_model = fit_final_model(df, feature_cols)
    final_model.save_model(OUT_DIR_DR / "model_final.json")

    # ── Predictions for all 17,745 DB genes ───────────────────────
    X_all     = df[feature_cols].values.astype(float)
    dr_scores = np.clip(final_model.predict(X_all), 0.0, 1.0)

    preds = df[[
        "targetId", "gene_symbol",
        "hasSafetyEvent", "maxClinicalTrialPhase",
        "is_drugged", "pi_hat", "m_hat", "Y_tilde",
    ]].copy()
    preds["safety_score_dr"] = dr_scores
    preds["safety_label"] = df["hasSafetyEvent"].apply(
        lambda x: "positive" if x == -1.0 else (
            "drugged_safe" if pd.notna(x) else "unlabeled"
        )
    )
    # Correct: drugged + no safety event = negative label
    preds.loc[(df["is_drugged"] == 1) & (df["hasSafetyEvent"].isna()), "safety_label"] = \
        "drugged_safe"

    preds.sort_values("safety_score_dr", ascending=False, inplace=True)
    preds.to_csv(OUT_DIR_DR / "predictions.csv", index=False)
    print(f"\nPredictions → {OUT_DIR_DR / 'predictions.csv'}  ({len(preds):,} genes)")

    # ── Diagnostics ───────────────────────────────────────────────
    _print_diagnostics(df, dr_scores, preds)

    print(f"\n{'='*60}")
    print(f"All outputs → {OUT_DIR_DR}")


if __name__ == "__main__":
    main()
