"""
Safety prediction model for biological targets.

Two complementary models are trained and evaluated:

MODEL A — All-labeled (semi-supervised / PU-like)
  Positive  : hasSafetyEvent == -1  (945 targets with confirmed safety events)
  Negative  : drugged target (maxClinicalTrialPhase notna) + hasSafetyEvent NaN (1,149)
  Unlabeled : not drugged, no safety record — predicted but NOT trained on
  Limitation: predictions for undrugged genes are extrapolations far outside the
              training distribution and should be interpreted with caution.

MODEL B — Drugged-only (focused, honest evaluation)
  Restricted to genes that have reached clinical trials (1,564 targets).
  Positive  : hasSafetyEvent == -1 within drugged subset (415)
  Negative  : hasSafetyEvent NaN within drugged subset (1,149)
  Predictions use out-of-fold (OOF) probabilities from 5-fold CV so that every
  drugged gene gets a score estimated from a model that never trained on it.
  This is the recommended model for clinical-stage target assessment.

Features (shared):
  A) Open Targets platform features (target_prioritisation)
  B) LOEUF genetic constraint
  C) pLoF burden: gene-level aggregates from Genebass
  D) Cell-type expression: per-organ max/mean from CellxGene

Usage:
    python data/safety_model.py

Outputs:
  data/safety_model_output/            Model A artefacts
    predictions.csv                    safety score for all OT genes (78K)
    shap_summary.csv                   mean |SHAP| per feature
    shap_values.npy                    full SHAP matrix for labeled set
    model.json                         XGBoost model

  data/safety_model_output/drugged_only/   Model B artefacts
    predictions_oof.csv                OOF safety scores for all 1,564 drugged genes
    shap_summary.csv                   mean |SHAP| per feature
    shap_values.npy                    full SHAP matrix
    model.json                         final XGBoost model (trained on all 1,564)
"""

import pickle
import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from sklearn.metrics import (
    average_precision_score,
    classification_report,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold

# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT             = Path(__file__).parent.parent
DATA_DIR         = ROOT / "data"
BACKEND          = ROOT / "backend"
OT_FILE          = DATA_DIR / "opentargets_target_prioritisation.parquet"
LOEUF_FILE       = DATA_DIR / "LOEUF_scores.csv.gz"
PLOF_PKL         = DATA_DIR / "genebass_pLoF_filtered.pkl"
DB_PATH          = BACKEND / "humanproof.db"
OUT_DIR          = DATA_DIR / "safety_model_output"
OUT_DIR_DRUGGED  = OUT_DIR / "drugged_only"
OUT_DIR.mkdir(exist_ok=True)
OUT_DIR_DRUGGED.mkdir(exist_ok=True)

# ─── Gene map (symbol → ensembl) ──────────────────────────────────────────────
GENES = {
    "TP53": "ENSG00000141510", "BRCA1": "ENSG00000012048", "BRCA2": "ENSG00000139618",
    "EGFR": "ENSG00000146648", "ERBB2": "ENSG00000141736", "KRAS": "ENSG00000133703",
    "BRAF": "ENSG00000157764", "PIK3CA": "ENSG00000121879", "PTEN": "ENSG00000171862",
    "MYC": "ENSG00000136997", "RB1": "ENSG00000139687", "APC": "ENSG00000134982",
    "VHL": "ENSG00000134086", "WT1": "ENSG00000184937", "NF1": "ENSG00000196712",
    "NF2": "ENSG00000186575", "RET": "ENSG00000165731", "MET": "ENSG00000105976",
    "ALK": "ENSG00000171094", "FGFR1": "ENSG00000077782", "FGFR2": "ENSG00000066468",
    "FGFR3": "ENSG00000068078", "PDGFRA": "ENSG00000134853", "KIT": "ENSG00000157404",
    "ABL1": "ENSG00000097007", "JAK2": "ENSG00000096968", "CDK4": "ENSG00000135446",
    "CDK6": "ENSG00000105810", "CCND1": "ENSG00000110092", "MDM2": "ENSG00000135679",
    "BCL2": "ENSG00000171791", "BAX": "ENSG00000087088", "CASP3": "ENSG00000164305",
    "TNF": "ENSG00000232810", "IL6": "ENSG00000136244", "IL1B": "ENSG00000125538",
    "IL2": "ENSG00000109471", "IL10": "ENSG00000136634", "IFNG": "ENSG00000111537",
    "TGFB1": "ENSG00000105329", "VEGFA": "ENSG00000112715", "FLT1": "ENSG00000102755",
    "KDR": "ENSG00000128052", "PDCD1": "ENSG00000188389", "CD274": "ENSG00000120217",
    "CTLA4": "ENSG00000163599", "CD28": "ENSG00000178562", "FOXP3": "ENSG00000049768",
    "CD4": "ENSG00000010610", "CD8A": "ENSG00000153563", "HLA-A": "ENSG00000206503",
    "B2M": "ENSG00000166710", "PCNA": "ENSG00000132646", "TOP2A": "ENSG00000131747",
    "ESR1": "ENSG00000091831", "AR": "ENSG00000169083", "ERBB3": "ENSG00000065361",
    "ERBB4": "ENSG00000178568", "IGF1R": "ENSG00000140443", "INSR": "ENSG00000171105",
    "SRC": "ENSG00000197122", "RAF1": "ENSG00000132155", "MAP2K1": "ENSG00000169032",
    "MAPK1": "ENSG00000100030", "MAPK3": "ENSG00000102882", "AKT1": "ENSG00000142208",
    "MTOR": "ENSG00000198793", "NOTCH1": "ENSG00000148400", "WNT1": "ENSG00000125084",
    "SHH": "ENSG00000164690", "PTCH1": "ENSG00000185920", "SMO": "ENSG00000128602",
    "STAT3": "ENSG00000168610", "NFKB1": "ENSG00000109320", "RELA": "ENSG00000173039",
    "CREBBP": "ENSG00000005339", "EP300": "ENSG00000100393", "HDAC1": "ENSG00000116478",
    "DNMT1": "ENSG00000130816", "TET2": "ENSG00000168769", "IDH1": "ENSG00000138413",
    "IDH2": "ENSG00000182054", "ARID1A": "ENSG00000117713", "SMAD4": "ENSG00000141646",
    "CDKN2A": "ENSG00000147889", "CDKN1A": "ENSG00000124762", "GAPDH": "ENSG00000111640",
    "ACTB": "ENSG00000075624", "ALB": "ENSG00000163631", "INS": "ENSG00000254647",
    "GCG": "ENSG00000115263", "EPO": "ENSG00000130427", "THPO": "ENSG00000090534",
    "CSF3": "ENSG00000108342", "APOE": "ENSG00000130203", "APP": "ENSG00000142192",
    "MAPT": "ENSG00000186868", "SNCA": "ENSG00000145335", "HTT": "ENSG00000197386",
    "SOD1": "ENSG00000142168", "CFTR": "ENSG00000001626", "DMD": "ENSG00000198947",
    "F8": "ENSG00000185010", "HBB": "ENSG00000244734", "HBA1": "ENSG00000206172",
    "FGA": "ENSG00000171560", "CYP3A4": "ENSG00000160868", "CYP2D6": "ENSG00000100197",
}
ENSG_TO_SYMBOL = {v: k for k, v in GENES.items()}

PLOF_ORGANS = [
    "cardiovascular", "hepatic", "neurological", "renal",
    "respiratory", "hematologic", "musculoskeletal",
]  # legacy constant — only used by the 108-gene Model A/B loader below


# ─────────────────────────────────────────────────────────────────────────────
# Feature loaders
# ─────────────────────────────────────────────────────────────────────────────

def load_opentargets() -> pd.DataFrame:
    """Load Open Targets target_prioritisation features."""
    print("Loading Open Targets data...")
    df = pd.read_parquet(OT_FILE)

    # Encode hasSafetyEvent: -1 → 1 (positive), NaN → NaN (unknown)
    df["safety_label"] = df["hasSafetyEvent"].apply(
        lambda x: 1 if x == -1.0 else (0 if pd.notna(x) else np.nan)
    )
    # Trusted negative: has been drugged but no safety event
    df["is_drugged"]   = df["maxClinicalTrialPhase"].notna().astype(int)
    df["trusted_neg"]  = (df["is_drugged"] == 1) & (df["hasSafetyEvent"].isna())

    # Recode isCancerDriverGene: -1 → 1
    df["isCancerDriverGene"] = df["isCancerDriverGene"].apply(
        lambda x: 1.0 if x == -1.0 else (0.0 if pd.notna(x) else np.nan)
    )

    print(f"  {len(df):,} total targets")
    print(f"  Safety positives: {(df['hasSafetyEvent'] == -1.0).sum()}")
    print(f"  Trusted negatives (drugged, no safety): {df['trusted_neg'].sum()}")
    return df


def load_loeuf_features() -> pd.DataFrame:
    """Load LOEUF per gene → map to Ensembl ID."""
    print("Loading LOEUF features...")
    loeuf_df = pd.read_csv(LOEUF_FILE)
    loeuf_map = dict(zip(loeuf_df["gene"], loeuf_df["LOEUF"]))

    rows = []
    for symbol, ensg in GENES.items():
        rows.append({
            "targetId":    ensg,
            "gene_symbol": symbol,
            "loeuf_score": loeuf_map.get(symbol, np.nan),
        })
    return pd.DataFrame(rows)


def load_plof_features() -> pd.DataFrame:
    """Aggregate Genebass pLoF data to gene-level features.

    Features per gene:
      plof_min_pval        : minimum p-value across all phenotypes
      plof_n_sig           : number of genome-wide significant (p<5e-8) associations
      plof_max_abs_beta    : maximum |beta| across phenotypes
      plof_{organ}_min_p   : minimum p-value for each major organ system
      plof_{organ}_max_beta: strongest beta for each organ system
    """
    print("Loading pLoF features...")
    df = pd.read_pickle(PLOF_PKL)
    df = df[df["gene"].isin(set(GENES.keys()))].copy()
    df = df.dropna(subset=["Pvalue"])
    df["BETA_Burden"] = df["BETA_Burden"].fillna(0)

    # Classify into organ systems (simplified inline)
    organ_keywords = {
        "cardiovascular": ["heart","cardiac","coronary","hypertension","blood pressure",
                           "stroke","vascular","arrhythmia","cardiomyopathy","aortic"],
        "hepatic":        ["liver","hepat","bilirubin","alt","ast","ggt","cirrhosis"],
        "neurological":   ["brain","cerebral","alzheimer","parkinson","dementia",
                           "epilepsy","schizophrenia","depression","anxiety","cognitive"],
        "renal":          ["kidney","renal","creatinine","glomerular","nephr","urea"],
        "respiratory":    ["lung","pulmonary","asthma","copd","fev1","fvc","respiratory"],
        "hematologic":    ["haemoglobin","hemoglobin","platelet","erythrocyte",
                           "red blood cell","white blood cell","anaemia","anemia"],
        "musculoskeletal":["bone","fracture","osteopor","arthritis","joint","muscle",
                           "grip strength","height","muscular"],
    }

    def get_organ(desc):
        d = str(desc).lower()
        for organ, kws in organ_keywords.items():
            if any(kw in d for kw in kws):
                return organ
        return "other"

    df["organ"] = df["pheno_description"].apply(get_organ)

    records = []
    for gene, grp in df.groupby("gene"):
        ensg = GENES[gene]
        row = {"targetId": ensg, "gene_symbol": gene}
        row["plof_min_pval"]     = grp["Pvalue"].min()
        row["plof_n_sig"]        = (grp["Pvalue"] < 5e-8).sum()
        row["plof_max_abs_beta"] = grp["BETA_Burden"].abs().max()
        row["plof_n_phenotypes"] = len(grp)

        for organ in PLOF_ORGANS:
            sub = grp[grp["organ"] == organ]
            row[f"plof_{organ}_min_p"]    = sub["Pvalue"].min() if len(sub) else np.nan
            row[f"plof_{organ}_max_beta"] = sub["BETA_Burden"].abs().max() if len(sub) else np.nan

        records.append(row)

    return pd.DataFrame(records)


def _ct_col(ct: str) -> str:
    """Sanitize a cell-type label into a valid column name."""
    return "expr_ct_" + re.sub(r"[^a-z0-9]+", "_", ct.lower()).strip("_")


def load_expression_features() -> pd.DataFrame:
    """Load per-gene expression for ALL cell types in the DB.

    Uses dynamic pivot — no hardcoded cell-type filter.
    Column naming: expr_ct_<sanitized_cell_type> (mean) and ..._pct (% expressing).
    """
    print("Loading expression features...")

    if not DB_PATH.exists():
        raise FileNotFoundError(f"DB not found: {DB_PATH}. Run backend/load_real_data.py first.")

    conn = sqlite3.connect(DB_PATH)
    expr = pd.read_sql(
        "SELECT gene_symbol, cell_type, organ, mean_expression, pct_expressed "
        "FROM expression_summary",
        conn,
    )
    conn.close()

    # Restrict to the 108 target genes
    expr = expr[expr["gene_symbol"].isin(set(GENES.keys()))].copy()

    expr["col"] = expr["cell_type"].map(_ct_col)

    mean_pivot = (
        expr.pivot_table(
            index="gene_symbol", columns="col",
            values="mean_expression", aggfunc="mean",
        )
        .reset_index()
    )
    mean_pivot.columns.name = None

    expr["col_pct"] = expr["col"] + "_pct"
    pct_pivot = (
        expr.pivot_table(
            index="gene_symbol", columns="col_pct",
            values="pct_expressed", aggfunc="mean",
        )
        .reset_index()
    )
    pct_pivot.columns.name = None

    organ_max = (
        expr
        .groupby(["gene_symbol", "organ"])["mean_expression"]
        .max()
        .unstack("organ")
        .add_prefix("expr_organ_")
        .reset_index()
    )

    base_df = mean_pivot.merge(pct_pivot, on="gene_symbol", how="outer")
    base_df = base_df.merge(organ_max, on="gene_symbol", how="left")

    # Map gene_symbol → targetId (Ensembl)
    base_df["targetId"] = base_df["gene_symbol"].map(GENES)
    return base_df[base_df["targetId"].notna()].copy()


# ─────────────────────────────────────────────────────────────────────────────
# Build feature matrix
# ─────────────────────────────────────────────────────────────────────────────

def build_feature_matrix() -> pd.DataFrame:
    ot       = load_opentargets()
    loeuf    = load_loeuf_features()
    plof     = load_plof_features()
    expr     = load_expression_features()

    # Merge: OT is the master table (all 78K genes)
    df = ot.copy()
    df = df.merge(loeuf.drop(columns=["gene_symbol"]), on="targetId", how="left")
    df = df.merge(plof.drop(columns=["gene_symbol"]),  on="targetId", how="left")
    df = df.merge(expr.drop(columns=["gene_symbol"]),  on="targetId", how="left")

    # Add gene symbol where known
    df["gene_symbol"] = df["targetId"].map(ENSG_TO_SYMBOL)

    print(f"\nFull feature matrix: {df.shape}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Model training
# ─────────────────────────────────────────────────────────────────────────────

OT_FEATURES = [
    "isInMembrane", "isSecreted", "hasPocket", "hasLigand", "hasSmallMoleculeBinder",
    "geneticConstraint", "paralogMaxIdentityPercentage", "mouseOrthologMaxIdentityPercentage",
    "isCancerDriverGene", "hasTEP", "mouseKOScore", "hasHighQualityChemicalProbes",
    "tissueSpecificity", "tissueDistribution",
]

LOEUF_FEATURES = ["loeuf_score"]

PLOF_BASE_FEATURES = [
    "plof_min_pval", "plof_n_sig", "plof_max_abs_beta", "plof_n_phenotypes",
]
def get_feature_cols(df: pd.DataFrame) -> list[str]:
    """Return feature columns that exist in df.

    OT and LOEUF features are checked by name.
    pLoF base features are checked by name; pLoF organ features (plof_*_min_p /
    plof_*_max_beta) are discovered dynamically from df.columns — all phenotype
    categories in the DB are included, no hardcoded organ filter.
    Expression features (expr_ct_* and expr_organ_*) are also discovered dynamically.
    """
    static = OT_FEATURES + LOEUF_FEATURES + PLOF_BASE_FEATURES
    plof_organ_cols = sorted(
        c for c in df.columns
        if c.startswith("plof_") and c not in set(PLOF_BASE_FEATURES)
    )
    expr_cols = sorted(
        c for c in df.columns
        if c.startswith("expr_ct_") or c.startswith("expr_organ_")
    )
    return [c for c in static if c in df.columns] + plof_organ_cols + expr_cols


def log_transform_pvals(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Replace p-value columns with -log10(p), clipping at 50."""
    df = df.copy()
    for col in cols:
        if col in df.columns:
            df[col] = -np.log10(df[col].clip(lower=1e-300)).clip(upper=50)
    return df


def train_and_evaluate(df: pd.DataFrame) -> tuple[xgb.XGBClassifier, list[str], pd.DataFrame]:
    """Train XGBoost safety classifier."""

    # ── Label construction ───────────────────────────────────────
    # Positives: hasSafetyEvent == -1
    # Trusted negatives: drugged (any clinical phase) + no safety event
    # Unlabeled: rest (excluded from training)
    pos_mask = df["hasSafetyEvent"] == -1.0
    neg_mask = df["trusted_neg"]
    labeled  = pos_mask | neg_mask

    df_labeled = df[labeled].copy()
    df_labeled["y"] = pos_mask[labeled].astype(int)

    print(f"\nLabeled training set: {len(df_labeled):,}")
    print(f"  Positives: {df_labeled['y'].sum()}")
    print(f"  Negatives: {(df_labeled['y']==0).sum()}")

    # ── Features ─────────────────────────────────────────────────
    pval_cols = [c for c in df_labeled.columns if "min_p" in c]
    df_labeled = log_transform_pvals(df_labeled, pval_cols)

    feature_cols = get_feature_cols(df_labeled)
    print(f"\nFeatures: {len(feature_cols)}")

    X = df_labeled[feature_cols].values
    y = df_labeled["y"].values

    # Scale pos weight for imbalance
    n_neg = (y == 0).sum()
    n_pos = (y == 1).sum()
    spw   = n_neg / n_pos
    print(f"  scale_pos_weight: {spw:.2f}")

    # ── Cross-validation ─────────────────────────────────────────
    cv_auc, cv_ap = [], []
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    model_params = dict(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        scale_pos_weight=spw,
        random_state=42,
        eval_metric="aucpr",
        early_stopping_rounds=30,
        tree_method="hist",
    )

    print("\n5-fold cross-validation...")
    for fold, (tr, va) in enumerate(skf.split(X, y)):
        m = xgb.XGBClassifier(**model_params)
        m.fit(
            X[tr], y[tr],
            eval_set=[(X[va], y[va])],
            verbose=False,
        )
        prob = m.predict_proba(X[va])[:, 1]
        auc  = roc_auc_score(y[va], prob)
        ap   = average_precision_score(y[va], prob)
        cv_auc.append(auc)
        cv_ap.append(ap)
        print(f"  Fold {fold+1}: AUROC={auc:.3f}  AUPRC={ap:.3f}")

    print(f"\n  Mean AUROC : {np.mean(cv_auc):.3f} ± {np.std(cv_auc):.3f}")
    print(f"  Mean AUPRC : {np.mean(cv_ap):.3f} ± {np.std(cv_ap):.3f}")

    # ── Final model on all labeled data ──────────────────────────
    print("\nTraining final model on all labeled data...")
    final_model = xgb.XGBClassifier(**{k: v for k, v in model_params.items()
                                        if k not in ("early_stopping_rounds",)})
    final_model.fit(X, y, verbose=False)

    return final_model, feature_cols, df_labeled


def predict_all(
    model: xgb.XGBClassifier,
    df: pd.DataFrame,
    feature_cols: list[str],
) -> pd.DataFrame:
    """Generate safety score predictions for all genes."""
    pval_cols = [c for c in feature_cols if "min_p" in c]
    df_pred   = log_transform_pvals(df, pval_cols)

    X_all = df_pred[feature_cols].values
    proba = model.predict_proba(X_all)[:, 1]

    out = df[["targetId", "gene_symbol", "hasSafetyEvent", "maxClinicalTrialPhase",
              "is_drugged"]].copy()
    out["safety_score"]  = proba
    out["safety_label"]  = df["hasSafetyEvent"].apply(
        lambda x: "positive" if x == -1.0 else ("drugged_safe" if pd.notna(x) else "unlabeled")
    )
    # Correct: drugged+no-safety = negative label
    out.loc[df["trusted_neg"], "safety_label"] = "drugged_safe"
    out.sort_values("safety_score", ascending=False, inplace=True)
    return out


def shap_analysis(
    model: xgb.XGBClassifier,
    df_labeled: pd.DataFrame,
    feature_cols: list[str],
    out_dir: Path,
) -> pd.DataFrame:
    """Compute SHAP values and return mean |SHAP| summary."""
    print("\nComputing SHAP values...")
    pval_cols = [c for c in feature_cols if "min_p" in c]
    df_s = log_transform_pvals(df_labeled, pval_cols)
    X    = df_s[feature_cols].values

    explainer   = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)

    np.save(out_dir / "shap_values.npy", shap_values)

    summary = pd.DataFrame({
        "feature":       feature_cols,
        "mean_abs_shap": np.abs(shap_values).mean(axis=0),
    }).sort_values("mean_abs_shap", ascending=False)

    print("\nTop 20 features by mean |SHAP|:")
    print(summary.head(20).to_string(index=False))
    return summary


# ─────────────────────────────────────────────────────────────────────────────
# Model B: drugged-only
# ─────────────────────────────────────────────────────────────────────────────

def train_drugged_only(df: pd.DataFrame) -> tuple[xgb.XGBClassifier, list[str], pd.DataFrame]:
    """Train XGBoost restricted to clinically-tested targets only.

    Population: genes with maxClinicalTrialPhase notna (1,564 targets).
    Label     : hasSafetyEvent == -1 → 1 (positive);  NaN → 0 (negative).
    Predictions are out-of-fold (OOF) probabilities from 5-fold CV so every
    gene receives a score from a model that never trained on it.
    """
    drugged = df[df["maxClinicalTrialPhase"].notna()].copy()
    drugged["y"] = (drugged["hasSafetyEvent"] == -1.0).astype(int)

    n_pos = drugged["y"].sum()
    n_neg = (drugged["y"] == 0).sum()
    print(f"\n{'='*60}")
    print("MODEL B — Drugged-only")
    print(f"{'='*60}")
    print(f"  Drugged targets : {len(drugged):,}")
    print(f"  Positives (safety event): {n_pos}")
    print(f"  Negatives (no event)    : {n_neg}")

    pval_cols   = [c for c in drugged.columns if "min_p" in c]
    drugged     = log_transform_pvals(drugged, pval_cols)
    feature_cols= get_feature_cols(drugged)
    print(f"  Features: {len(feature_cols)}")

    X = drugged[feature_cols].values
    y = drugged["y"].values

    spw = n_neg / n_pos
    print(f"  scale_pos_weight: {spw:.2f}")

    model_params = dict(
        n_estimators=500,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        scale_pos_weight=spw,
        random_state=42,
        eval_metric="aucpr",
        early_stopping_rounds=30,
        tree_method="hist",
    )

    # ── OOF predictions (5-fold) ─────────────────────────────────
    oof_proba = np.zeros(len(drugged))
    cv_auc, cv_ap = [], []
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    print("\n5-fold CV (OOF predictions)...")
    for fold, (tr, va) in enumerate(skf.split(X, y)):
        m = xgb.XGBClassifier(**model_params)
        m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])], verbose=False)
        prob = m.predict_proba(X[va])[:, 1]
        oof_proba[va] = prob
        auc = roc_auc_score(y[va], prob)
        ap  = average_precision_score(y[va], prob)
        cv_auc.append(auc)
        cv_ap.append(ap)
        print(f"  Fold {fold+1}: AUROC={auc:.3f}  AUPRC={ap:.3f}")

    print(f"\n  Mean AUROC : {np.mean(cv_auc):.3f} ± {np.std(cv_auc):.3f}")
    print(f"  Mean AUPRC : {np.mean(cv_ap):.3f} ± {np.std(cv_ap):.3f}")

    # OOF AUROC / AUPRC on full drugged set
    oof_auc = roc_auc_score(y, oof_proba)
    oof_ap  = average_precision_score(y, oof_proba)
    print(f"\n  OOF (pooled) AUROC: {oof_auc:.3f}")
    print(f"  OOF (pooled) AUPRC: {oof_ap:.3f}")

    # ── Final model on all drugged data (for SHAP) ───────────────
    print("\nTraining final model on all drugged genes...")
    final_model = xgb.XGBClassifier(**{k: v for k, v in model_params.items()
                                        if k not in ("early_stopping_rounds",)})
    final_model.fit(X, y, verbose=False)

    # ── OOF predictions output ────────────────────────────────────
    out = drugged[["targetId", "gene_symbol", "hasSafetyEvent",
                   "maxClinicalTrialPhase", "y"]].copy()
    out = out.rename(columns={"y": "true_label"})
    out["safety_score_oof"] = oof_proba
    out["safety_label"] = out["true_label"].map(
        {1: "positive", 0: "drugged_safe"}
    )
    out = out.sort_values("safety_score_oof", ascending=False)

    return final_model, feature_cols, drugged, out


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("HumanProof Safety Prediction Model")
    print("=" * 60)

    df = build_feature_matrix()

    # ── MODEL A: all-labeled (semi-supervised) ────────────────────
    print(f"\n{'='*60}")
    print("MODEL A — All-labeled (semi-supervised)")
    print(f"{'='*60}")
    model_a, feat_a, df_labeled_a = train_and_evaluate(df)

    model_a.save_model(OUT_DIR / "model.json")
    print(f"\nModel A saved → {OUT_DIR / 'model.json'}")

    preds_a = predict_all(model_a, df, feat_a)
    preds_a.to_csv(OUT_DIR / "predictions.csv", index=False)
    print(f"Predictions saved → {OUT_DIR / 'predictions.csv'}  ({len(preds_a):,} genes)")

    shap_a = shap_analysis(model_a, df_labeled_a, feat_a, OUT_DIR)
    shap_a.to_csv(OUT_DIR / "shap_summary.csv", index=False)
    print(f"SHAP summary saved → {OUT_DIR / 'shap_summary.csv'}")

    # ── MODEL B: drugged-only (honest evaluation) ─────────────────
    model_b, feat_b, df_drugged, preds_b = train_drugged_only(df)

    model_b.save_model(OUT_DIR_DRUGGED / "model.json")
    print(f"\nModel B saved → {OUT_DIR_DRUGGED / 'model.json'}")

    preds_b.to_csv(OUT_DIR_DRUGGED / "predictions_oof.csv", index=False)
    print(f"OOF predictions saved → {OUT_DIR_DRUGGED / 'predictions_oof.csv'}  ({len(preds_b):,} drugged genes)")

    shap_b = shap_analysis(model_b, df_drugged, feat_b, OUT_DIR_DRUGGED)
    shap_b.to_csv(OUT_DIR_DRUGGED / "shap_summary.csv", index=False)
    print(f"SHAP summary saved → {OUT_DIR_DRUGGED / 'shap_summary.csv'}")

    # ── Summary for our 108 genes ─────────────────────────────────
    print(f"\n{'='*60}")
    print("Our 108 genes — Model B (drugged-only OOF scores)")
    print(f"{'='*60}")
    our_b = preds_b[preds_b["gene_symbol"].notna()].copy()
    print(our_b[["gene_symbol","safety_score_oof","safety_label",
                 "maxClinicalTrialPhase"]].to_string(index=False))

    undrugged_108 = [g for g in GENES if g not in our_b["gene_symbol"].values]
    if undrugged_108:
        print(f"\n  Not in drugged set ({len(undrugged_108)} genes — Model A scores only):")
        our_a = preds_a[preds_a["gene_symbol"].isin(undrugged_108)].copy()
        print(our_a[["gene_symbol","safety_score","safety_label"]].to_string(index=False))

    print(f"\n{'='*60}")
    print(f"All outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
