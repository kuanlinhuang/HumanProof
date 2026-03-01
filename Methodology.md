# HumanProof — Scientific Methodology

This document describes the data acquisition, processing, feature engineering, and statistical modelling used by HumanProof to produce genome-wide drug target safety risk scores.

---

## Table of Contents

1. [Gene Universe](#1-gene-universe)
2. [Data Sources and Processing](#2-data-sources-and-processing)
   - 2.1 [CZ CellxGene Census — Expression](#21-cz-cellxgene-census--expression)
   - 2.2 [Genebass UK Biobank — Protein-Truncating Variant Burden Tests](#22-genebass-uk-biobank--protein-truncating-variant-burden-tests)
   - 2.3 [LOEUF Constraint Scores (gnomAD)](#23-loeuf-constraint-scores-gnomaD)
   - 2.4 [Open Targets Platform 25.12](#24-open-targets-platform-2512)
3. [Feature Engineering](#3-feature-engineering)
   - 3.1 [Expression Features](#31-expression-features)
   - 3.2 [pLoF Burden Features](#32-plof-burden-features)
   - 3.3 [LOEUF Feature](#33-loeuf-feature)
   - 3.4 [Open Targets Features](#34-open-targets-features)
   - 3.5 [Feature Matrix Assembly](#35-feature-matrix-assembly)
   - 3.6 [p-value Transformation](#36-p-value-transformation)
4. [DR+PU Safety Prediction Model](#4-drpu-safety-prediction-model)
   - 4.1 [Problem Formulation](#41-problem-formulation)
   - 4.2 [Step 1 — Cross-fitted Propensity Model P(S=1|X)](#42-step-1--cross-fitted-propensity-model-ps1x)
   - 4.3 [Step 2 — PU Prior Estimation via Elkan–Noto](#43-step-2--pu-prior-estimation-via-elkannoto)
   - 4.4 [Step 3 — Cross-fitted Outcome Model P(Y=1|X)](#44-step-3--cross-fitted-outcome-model-py1x)
   - 4.5 [Step 4 — AIPW Pseudo-outcome Construction](#45-step-4--aipw-pseudo-outcome-construction)
   - 4.6 [Step 5 — Final DR Regression m_DR(X)](#46-step-5--final-dr-regression-m_drx)
   - 4.7 [SHAP Interpretability](#47-shap-interpretability)
   - 4.8 [Statistical Guarantees](#48-statistical-guarantees)
5. [Label Definitions](#5-label-definitions)
6. [Model Hyperparameters](#6-model-hyperparameters)
7. [Implementation Details](#7-implementation-details)
8. [References](#8-references)

---

## 1. Gene Universe

The analysis is anchored to the set of **protein-coding genes** for which constraint information is available from gnomAD, as distributed via the LOEUF score file (`data/LOEUF_scores.csv.gz`). After removing genes with missing LOEUF and collapsing duplicate symbols to the most-constrained entry, this yields **~19,155 unique protein-coding genes** as the master universe.

A working subset of **17,745 genes** is used for modelling: the intersection of the LOEUF universe with genes present in the CZ CellxGene expression matrices, resolved via the CellxGene `gene_metadata.csv` file (which maps matrix column names to Ensembl IDs). For genes represented by multiple Ensembl IDs in the matrix (suffix `_ENSGXXXXXXXX`), the first encountered mapping is used.

The model training and scoring population is further intersected with the Open Targets Platform `target_prioritisation` parquet (78,725 targets), which provides the outcome variable (`hasSafetyEvent`) and clinical development status (`maxClinicalTrialPhase`). All 17,745 DB genes are present in the OT parquet.

---

## 2. Data Sources and Processing

### 2.1 CZ CellxGene Census — Expression

**Source:** CZ CellxGene single-cell census, pre-aggregated per cell type across all normal human tissue donors. Files: `data/cellxgene/`.

**Pseudobulk aggregation (pre-computed):**
The census provides two aggregation matrices over all donors of the same tissue and cell type:

| File | Metric | Description |
|------|--------|-------------|
| `celltype_log1p_mean_expression.csv.gz` | `mean_expression` | log1p of mean raw UMI count per (cell type, gene) across all cells of that type |
| `celltype_fraction_expressing.csv.gz` | `pct_expressed` | Fraction of cells with raw count > 0 |
| `celltype_metadata.csv` | `n_cells` | Total number of cells of that type across all sampled tissues |

The matrices are wide-format (rows = cell types, columns = gene symbols). These are **not** computed by HumanProof — they are downloaded as pre-computed aggregations from the CellxGene census.

`mean_expression = log1p(mean_count)` where `mean_count` is the arithmetic mean raw UMI count over all cells of a given type. log1p compression reduces dynamic range and is standard in single-cell analysis.

**Cell type selection:** 60 biologically diverse cell types across 14 tissue compartments, curated in `CELLXGENE_TISSUE_MAP` in `backend/load_real_data.py`. Cell types were selected for biological relevance to drug safety (covering major toxicity-relevant organs: heart, liver, kidney, brain, lung) and data quality (total cells ≥ ~5,000 in the census). All 60 are included as features with no downstream filtering — see §3.1.

| Tissue | Cell types included |
|--------|---------------------|
| Brain | Neurons, glutamatergic neurons, GABAergic neurons, interneurons, oligodendrocytes, astrocytes, microglia |
| Heart | Cardiomyocytes, cardiac endothelial, cardiac fibroblasts |
| Liver | Hepatocytes, Kupffer cells, hepatic stellate cells, cholangiocytes |
| Lung | Type I/II pneumocytes, alveolar macrophages, ciliated epithelium, club cells, lung macrophages |
| Kidney | Proximal tubule, tubular cells, podocytes, collecting duct, distal tubule |
| Intestine | Enterocytes, goblet cells, enteroendocrine cells, intestinal stem cells |
| Bone marrow | HSCs, erythroblasts, erythrocytes, megakaryocytes, neutrophils, plasma cells |
| Lymph node | T cells, CD4+ T cells, CD8+ T cells, Tregs, B cells, NK cells, dendritic cells |
| Spleen | Macrophages, monocytes |
| Skin | Keratinocytes, melanocytes, skin fibroblasts, dermal fibroblasts |
| Pancreas | Beta cells, alpha cells, acinar cells, ductal cells |
| Muscle | Skeletal muscle fibers, skeletal muscle cells, satellite cells |
| Adrenal | Adrenal cortical cells, chromaffin cells |
| Breast | Mammary epithelial, luminal epithelial, mammary fibroblasts |

**Data quality handling:** Non-finite values (NaN, Inf) in `mean_expression` are set to 0.0. `pct_expressed` is clipped to [0, 1].

**Database storage:** Records are stored long-format in the `expression_summary` SQLite table: one row per (gene × cell type), with `gene_symbol`, `ensembl_id`, `cell_type`, `tissue`, `organ`, `mean_expression`, `pct_expressed`, `n_cells`.

---

### 2.2 Genebass UK Biobank — Protein-Truncating Variant Burden Tests

**Source:** Genebass (genebass.org), exome-wide burden test results from the UK Biobank (~400,000 whole-exome sequenced participants).

**Test methodology:** Each association is a gene-based rare-variant burden test for protein-truncating variants (pLoF: stop-gained, frameshift, splice-site-disrupting). Results include:
- `Pvalue`: combined burden/SKAT omnibus p-value
- `Pvalue_Burden`: burden-only p-value
- `Pvalue_SKAT`: SKAT p-value
- `BETA_Burden`: burden test effect size (β), expressed in units of phenotype SD per minor allele
- `SE_Burden`: standard error of β

**Filtering and truncation:** The full Genebass pLoF file (`genebass_pLoF_filtered.pkl`) contains results filtered to a subset of genes. For DB construction, the top 200 associations per gene (by ascending `Pvalue`) are retained to keep the database size manageable (~18,092 genes × 200 = ~3.6M rows). Associations with missing `Pvalue` are dropped; Inf/−Inf values in effect sizes are set to NaN.

**Phenotype classification:** Each phenotype description (free-text from UK Biobank) is classified into one of 19 categories using a rule-based keyword matcher (`classify_phenotype()` in `backend/load_real_data.py`):

| Category | Example keywords |
|----------|-----------------|
| cardiovascular | heart, cardiac, hypertension, stroke, arrhythmia |
| metabolic | diabetes, cholesterol, BMI, obesity, thyroid |
| hematologic | haemoglobin, platelet, red blood cell, anaemia |
| neurological | brain, Alzheimer, Parkinson, depression, cognitive |
| hepatic | liver, ALT, AST, bilirubin, cirrhosis |
| renal | kidney, creatinine, glomerular, urea |
| respiratory | lung, asthma, COPD, FEV1, FVC |
| musculoskeletal | bone, fracture, arthritis, grip strength, height |
| cancer | cancer, carcinoma, lymphoma, leukaemia |
| immunologic | autoimmune, CRP, immunoglobulin, allergy |
| ophthalmologic | eye, retinal, glaucoma, macular |
| dermatologic | skin, acne, eczema, hair, pigment |
| reproductive | prostate, ovarian, fertility, pregnancy |
| gastrointestinal | gastric, bowel, colon, reflux, IBS |
| audiologic | hearing, tinnitus, cochlear |
| anthropometric | arm span, leg, comparative body size |
| medication | medication, treatment, prescribed, tablet |
| procedural | operation, surgery, procedure |
| other | unclassified |

**Database storage:** Records stored in `plof_associations` table: `gene_symbol`, `phenotype`, `phenotype_category` (lowercase, e.g. `cardiovascular`), `organ_system` (title-case), `p_value`, `p_value_burden`, `p_value_skat`, `beta`, `se`, `direction` (`loss` if β < 0, else `gain`).

---

### 2.3 LOEUF Constraint Scores (gnomAD)

**Source:** Loss-of-Function Observed/Expected Upper-bound Fraction, distributed via [grr.iossifovlab.com](https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/) (`data/LOEUF_scores.csv.gz`). Derived from the gnomAD v4 dataset.

**Definition:** LOEUF is the 90th-percentile upper bound of the ratio of observed to expected pLoF variants in gnomAD population sequencing data. A lower LOEUF indicates stronger purifying selection against heterozygous loss-of-function variants, reflecting gene essentiality:
- **LOEUF < 0.35** → `critical` (highly intolerant; equivalent to gnomAD pLI ≳ 0.9)
- **0.35 ≤ LOEUF < 0.70** → `high`
- **0.70 ≤ LOEUF < 1.00** → `moderate`
- **LOEUF ≥ 1.00** → `low` (tolerant to LoF)

**Pre-processing:** Genes with missing LOEUF are dropped. For the ~42 duplicated gene symbols in the file, the minimum (most constrained) LOEUF is retained per symbol. This yields ~19,155 unique entries.

**Ensembl ID mapping:** LOEUF scores are keyed by gene symbol. Ensembl IDs are filled from the CellxGene `gene_universe` mapping where available; genes absent from CellxGene receive an empty Ensembl ID but are still included in the LOEUF-only DB table.

**Note:** pLI and missense z-scores are not available from this source; these fields are stored as `0.0` placeholders in the DB.

---

### 2.4 Open Targets Platform 25.12

**Source:** Open Targets Platform release 25.12, `target_prioritisation` dataset. File: `data/opentargets_target_prioritisation.parquet` (Snappy-compressed Parquet, ~908 KB).

**Coverage:** 78,725 human gene targets × 17 curated features.

**Encoding conventions:**
- Binary druggability flags (`isInMembrane`, `isSecreted`, `hasPocket`, `hasLigand`, `hasSmallMoleculeBinder`, `hasTEP`, `hasHighQualityChemicalProbes`): `1.0` = true, `0.0` = false, `NaN` = not assessed.
- `hasSafetyEvent`: Open Targets convention uses `−1.0` to indicate a confirmed safety event. `NaN` means no safety event recorded (not confirmed safe). This asymmetric encoding reflects the positive-unlabeled nature of the label — absence of evidence is not evidence of absence.
- `isCancerDriverGene`: also uses `−1.0` = true; recoded to `1.0` during feature loading.
- `maxClinicalTrialPhase`: `0.25/0.50/0.75/1.0` encoding for Phases I–IV respectively; `NaN` = gene never entered clinical trials.
- Continuous features (`geneticConstraint`, `mouseKOScore`, `paralogMaxIdentityPercentage`, `mouseOrthologMaxIdentityPercentage`, `tissueSpecificity`, `tissueDistribution`): real-valued, possibly `NaN`.

**Role in the model:** The OT parquet serves as the backbone of the feature matrix (all 78,725 targets). The `hasSafetyEvent` and `maxClinicalTrialPhase` columns define the outcome label and training population respectively. All 14 OT features listed in §3.4 are used as model inputs.

---

## 3. Feature Engineering

### 3.1 Expression Features

From the `expression_summary` DB table, two feature types are derived per cell type:

**Per-cell-type mean expression** (`expr_ct_{cell_type}`):
Mean log1p expression for that gene in that cell type, directly from `mean_expression`. Column names are sanitized: `"expr_ct_" + re.sub(r"[^a-z0-9]+", "_", cell_type.lower()).strip("_")`. For example, `"cardiac muscle cell"` → `expr_ct_cardiac_muscle_cell`, `"CD4-positive, alpha-beta T cell"` → `expr_ct_cd4_positive_alpha_beta_t_cell`.

**Per-cell-type % expressing** (`expr_ct_{cell_type}_pct`):
Fraction of cells expressing the gene (count > 0) in that cell type. Column name = mean column + `_pct`.

**Organ-level max expression** (`expr_organ_{organ}`):
Maximum `mean_expression` across all cell types within the same organ (14 organs). Captures peak expression within each organ regardless of which cell type is most relevant.

This yields **120 cell-type features** (60 mean + 60 pct) + **14 organ features** = **134 expression features** total.

All expression feature columns are **discovered dynamically** from the DataFrame columns using `get_feature_cols()` — no hardcoded cell-type list exists downstream of the DB.

---

### 3.2 pLoF Burden Features

From the `plof_associations` DB table, gene-level aggregates are computed with SQL:

**Global aggregates (4 features):**

| Feature | Computation |
|---------|-------------|
| `plof_min_pval` | `MIN(p_value)` across all phenotypes |
| `plof_n_sig` | Count of associations with `p_value < 5×10⁻⁸` |
| `plof_max_abs_beta` | `MAX(ABS(beta))` across all phenotypes |
| `plof_n_phenotypes` | Total number of phenotype associations |

**Per-category aggregates (19 categories × 2 = 38 features):**

For each `phenotype_category` in the DB:

| Feature | Computation |
|---------|-------------|
| `plof_{category}_min_p` | `MIN(p_value)` within category |
| `plof_{category}_max_beta` | `MAX(ABS(beta))` within category |

Genes with no associations in a given category receive `NaN` (handled by XGBoost as missing). These columns are also discovered dynamically from DataFrame columns via `get_feature_cols()`.

Total pLoF features: **4 global + 38 per-category = 42**.

---

### 3.3 LOEUF Feature

A single feature `loeuf_score` (continuous, 0.03–2.0; lower = more constrained). Missing values (genes absent from gnomAD) are left as `NaN`.

---

### 3.4 Open Targets Features

14 features from the OT `target_prioritisation` parquet:

| Feature | Type | Description |
|---------|------|-------------|
| `isInMembrane` | binary | Protein localised to cell membrane |
| `isSecreted` | binary | Secreted protein |
| `hasPocket` | binary | Druggable pocket predicted |
| `hasLigand` | binary | Known endogenous or exogenous ligand |
| `hasSmallMoleculeBinder` | binary | Small molecule binder reported |
| `geneticConstraint` | continuous | OT-normalised genetic constraint score |
| `paralogMaxIdentityPercentage` | continuous | Maximum sequence identity to any paralogue |
| `mouseOrthologMaxIdentityPercentage` | continuous | Maximum mouse orthologue identity |
| `isCancerDriverGene` | binary | Annotated cancer driver gene (recoded −1→1) |
| `hasTEP` | binary | Target Enabling Package exists |
| `mouseKOScore` | continuous | Aggregate phenotypic severity from mouse knockout |
| `hasHighQualityChemicalProbes` | binary | Validated chemical probes available |
| `tissueSpecificity` | continuous | Degree of tissue-specific expression (OT score) |
| `tissueDistribution` | continuous | Breadth of tissue expression (OT score) |

---

### 3.5 Feature Matrix Assembly

The feature matrix is assembled by joining all four sources, using the Open Targets parquet as the master (backbone):

```
1. Load OT parquet → 78,725 rows
2. LEFT JOIN LOEUF by targetId (Ensembl)  → adds gene_symbol for ~17,745 DB genes
3. LEFT JOIN pLoF by gene_symbol           → adds 42 pLoF columns
4. LEFT JOIN expression by gene_symbol     → adds 134 expression columns
5. Restrict to rows with gene_symbol ≠ NaN → 17,745 DB genes
```

Total: **191 features** per gene. XGBoost handles missing values natively (learned split direction for NaN).

---

### 3.6 p-value Transformation

All p-value features (`plof_min_pval`, `plof_{category}_min_p`) are transformed to:

```
feature ← −log₁₀(p)  clipped to [0, 50]
```

with a floor of `1×10⁻³⁰⁰` before log to avoid underflow. This puts more weight on genome-wide significant signals (−log₁₀(5×10⁻⁸) ≈ 7.3) and compresses the tail. Transformation is applied in-place before any model fitting.

---

## 4. DR+PU Safety Prediction Model

### 4.1 Problem Formulation

**Notation:**
- X ∈ ℝ¹⁹¹ — gene features
- S ∈ {0, 1} — selection indicator: S=1 if gene has entered clinical trials (`maxClinicalTrialPhase` not null), S=0 otherwise
- Y ∈ {0, 1} — true safety liability (Y=1 if the gene has a confirmed drug safety event)
- Observed label: Y is observable only when S=1; for S=0 genes Y is latent

Two structural biases must be corrected simultaneously:

**1. Selection bias:** Genes that reach clinical trials (S=1) are a non-random, druggability-enriched subset. The propensity model AUROC ~0.91 confirms that features strongly predict selection. A naive model trained only on S=1 genes would learn to predict druggability, not safety liability.

**2. Positive-unlabeled (PU) structure:** Among S=1 genes, absence of a recorded safety event (`hasSafetyEvent` = NaN) does not mean the gene is safe — it may simply not have been observed or reported. Treating unlabeled genes as true negatives introduces systematic false-negative bias (underestimates risk for novel targets). The Elkan–Noto method estimates the fraction of all genes with true safety liability.

The DR+PU model addresses both simultaneously through a 5-step pipeline.

**Training population:** 17,745 DB genes; 1,506 drugged (S=1); 409 confirmed safety positives (Y=1 within S=1).

---

### 4.2 Step 1 — Cross-fitted Propensity Model P(S=1|X)

**Goal:** Estimate the probability that a gene enters clinical trials given its features.

**Method:** 5-fold stratified cross-fitting. In each fold, an XGBClassifier is trained on the 4/5 training partition of **all 17,745 DB genes** and predicts on the 1/5 held-out partition. The held-out predictions are assembled into a complete cross-fitted propensity vector π̂(X_i) such that every gene's score was estimated by a model that never saw it.

**Stabilized IPW weight:**

```
w_i = π̄ / max(π̂(X_i), ε),    ε = 1×10⁻⁶
```

where π̄ = mean(S) ≈ 0.085 (overall drugging rate among DB genes). Weights are clipped at **10×** (MAX_IPW = 10) to prevent extreme weights from dominating.

- Weights **up-weight** surprising drug targets (low π̂ but S=1) — genes that entered trials despite appearing undruggable.
- Weights **down-weight** predictable targets (high π̂ and S=1) — genes that entered trials for obvious druggability reasons.

**Effective sample size (ESS)** for the drugged subset:

```
ESS = (Σ w_i)² / Σ w_i²   (for S=1 genes only)
```

ESS ≈ 22% of the 1,506 drugged genes, confirming strong selection bias that would corrupt an unweighted model.

**Performance:** Propensity AUROC ~0.91 (mean across 5 folds).

---

### 4.3 Step 2 — PU Prior Estimation via Elkan–Noto

**Goal:** Estimate π_p = P(Y=1), the population prevalence of true safety liability across all protein-coding genes.

**Method (Elkan & Noto, 2008):** Train a naive classifier f that predicts the "labeled positive" indicator:

```
pu_label_i = 1  if  hasSafetyEvent_i == −1.0
             0  otherwise  (including unlabeled genes)
```

using all 17,745 DB genes, with 5-fold stratified CV and class-balance correction (`scale_pos_weight` = n_unlabeled / n_labeled). The mean classifier score on held-out labeled positives estimates the **labeling frequency**:

```
c = E[f(X) | Y=1] ≈ mean f(X_i)  over held-out labeled positives
```

The prior is then:

```
π_p = P(labeled=1) / c
```

where `P(labeled=1)` = proportion of all DB genes with a confirmed safety event.

**Observed values:** P(labeled=1) ≈ 0.052, c ≈ 0.576, **π_p ≈ 0.090**.

Interpretation: approximately 9% of all protein-coding genes are estimated to have true safety liability if pharmacologically inhibited — compared to the ~5% that have an observed recorded event among genes ever tested clinically.

π_p is clipped to [10⁻⁴, 0.50] to prevent degenerate estimates.

---

### 4.4 Step 3 — Cross-fitted Outcome Model P(Y=1|X)

**Goal:** Estimate the conditional probability of a safety event given features, corrected for both selection bias (via IPW) and label noise (via PU reweighting).

**Training population:** 1,506 drugged DB genes (S=1) only.

**Combined sample weight:**

```
w_i = IPW_i × PU_weight_i

where:
  PU_weight_i = π_p          if Y_i = 1  (confirmed safety event)
  PU_weight_i = (1 − π_p)    if Y_i = 0  (no safety event recorded — treated as noisy negative)
```

The PU weight down-weights "negative" drugged genes by (1−π_p) ≈ 0.91, reflecting that ~9% of them likely have true safety liability that was not observed. Positive labels are up-weighted by π_p.

**Additional class balance:** `scale_pos_weight` = n_neg / n_pos ≈ 2.68 (409 positives, 1,097 negatives within drugged set).

**Cross-fitting:** 5-fold stratified CV produces OOF predictions m̂_oof(X_i) for each drugged gene — scored by a model that never trained on it. This eliminates self-prediction bias and ensures the OOF scores are honest estimates of held-out performance.

A final outcome model is then trained on **all 1,506 drugged genes** (using the same weights, without early stopping) to enable extrapolation to undrugged genes in Step 4.

**Performance (OOF, pooled):** AUROC = 0.691, AUPRC = 0.426.

---

### 4.5 Step 4 — AIPW Pseudo-outcome Construction

**Goal:** Construct a doubly-robust corrected target Ỹ for every gene that absorbs both the selection bias and the OOF correction.

**Formula (AIPW estimator):**

```
Ỹ_i  =  m̂(X_i)  +  (S_i / π̂_i) × (Y_i − m̂(X_i))
```

**For drugged genes (S_i = 1):**

```
Ỹ_i  =  m̂_oof(X_i)  +  (Y_i − m̂_oof(X_i)) / max(π̂_i, 0.01)
```

- m̂(X_i) = OOF score from Step 3 (honest, not self-predicted)
- π̂_i clipped at lower bound 0.01 to prevent extreme correction
- Y_i ∈ {0, 1} binary safety label
- The IPW correction removes residual confounding: genes with unexpectedly high or low Y relative to their predicted m̂ are corrected upward or downward, weighted by how surprising their selection was

**For undrugged genes (S_i = 0):**

```
Ỹ_i  =  m̂(X_i)
```

- No IPW correction (S=0 so the second term vanishes)
- m̂(X_i) = prediction from the final outcome model extrapolated to undrugged genes

**Clipping:** Ỹ is clipped to [−1, 2] to limit the influence of extreme corrections from near-zero propensities.

**Empirical distribution:**
- Drugged genes Ỹ: mean ≈ 0.17, std ≈ 1.00
- Undrugged genes Ỹ: mean ≈ 0.02, std ≈ 0.06

The high variance in drugged-gene pseudo-outcomes reflects the genuine uncertainty in safety liability corrected for selection.

**Double robustness property:** The AIPW estimator is consistent if **either** the propensity model (Step 1) **or** the outcome model (Step 3) is correctly specified — a much weaker assumption than requiring both to be correct.

---

### 4.6 Step 5 — Final DR Regression m_DR(X)

**Goal:** Train a single, smooth predictor on all 17,745 genes that maps features to the AIPW-corrected pseudo-outcomes. This is the deployed model that scores every gene in the UI.

**Model:** XGBRegressor with squared-error objective, trained on all 17,745 DB genes with:
- **Target:** Ỹ from Step 4
- **Sample weights:** 2× for drugged genes (S=1), 1× for undrugged genes (S=0)
  - Drugged genes carry real label information Y through the AIPW correction; their pseudo-outcomes are anchored to observed events
  - Undrugged genes contribute only extrapolated m̂ — useful as "soft structure" but less informative than real labels
- **Fixed number of trees:** 300 (early stopping is not used for the final regressor, as optimising on noisy pseudo-outcomes is unreliable as an early-stopping criterion)

**Output transformation:** Predictions are clipped to [0, 1] to form the reported `safety_score_dr` for each gene.

**Properties of m_DR(X):**
- Produces a valid, de-biased risk score for any gene — drugged or not — with no model-switching logic
- SHAP values have a consistent, genome-wide interpretation across all 17,745 genes
- The regressor smooths over the noisy individual pseudo-outcomes, providing a more stable estimate than using Ỹ directly

---

### 4.7 SHAP Interpretability

SHAP (SHapley Additive exPlanations) values are computed for m_DR(X) using `shap.TreeExplainer`, which provides exact Shapley values for tree-based models in O(TLD) time (T = trees, L = leaves, D = depth).

For each gene, the SHAP decomposition:

```
m_DR(X_i) ≈ E[m_DR(X)] + Σ_j φ_j(X_i)
```

where φ_j(X_i) is the SHAP value for feature j: the average marginal contribution of feature j to the model output for gene i across all possible feature coalitions.

**Batch computation:** SHAP values are computed for all 17,745 DB genes simultaneously using vectorised TreeExplainer. Results are stored in `gene_shap_dr.json` (~455 MB), with per-gene records sorted by |φ_j| descending.

**Feature groups for display:**
- `ot` — Open Targets platform features (14)
- `genetics` — LOEUF constraint + all pLoF features (43 total)
- `expression` — all `expr_ct_*` and `expr_organ_*` features (134 total)

---

### 4.8 Statistical Guarantees

**Double robustness (Robins et al., 1994; Scharfstein et al., 1999):** The AIPW estimator satisfies

```
E[Ỹ − Y_true] → 0
```

if either the propensity model π̂(X) or the outcome model m̂(X) is consistent — not necessarily both. This makes the estimator robust to moderate misspecification in one of the two models.

**Cross-fitting (Chernozhukov et al., 2018):** Both the propensity (Step 1) and outcome (Step 3) models use 5-fold cross-fitting. Without cross-fitting, the in-sample residual (Y − m̂) could correlate with propensity estimation noise, introducing second-order bias that would break the double-robustness guarantee. Cross-fitting ensures that the residual used in the AIPW correction is evaluated on held-out data, eliminating this correlation.

**PU consistency (Elkan & Noto, 2008):** The Elkan–Noto prior correction ensures the outcome model's loss function is calibrated to the true positive prevalence π_p rather than the observed labeling rate. Without this correction, treating all unlabeled negatives as true negatives would systematically underestimate risk for undrugged genes.

**IPW stability:** Stabilized weights (π̄ / π̂ rather than 1/π̂) and hard clipping at 10× prevent any single gene from dominating the weighted training distribution.

---

## 5. Label Definitions

| Label | Condition | Count (DB genes) |
|-------|-----------|-----------------|
| `positive` | `hasSafetyEvent == −1.0` (OT convention) | ~409 (within drugged DB genes) |
| `drugged_safe` | `maxClinicalTrialPhase` not null AND `hasSafetyEvent` is null | ~1,097 |
| `unlabeled` | `maxClinicalTrialPhase` is null | ~16,239 |

The `drugged_safe` label is an **assumed negative** — the gene has been tested clinically and no safety event has been recorded. It is not a confirmed-safe label; the PU weighting in Step 3 acknowledges that ~9% of these genes may still have true safety liability.

---

## 6. Model Hyperparameters

All XGBoost models share a common base configuration (`_XGB_COMMON`):

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `n_estimators` | 500 (classifiers) / 300 (final regressor) | 500 with early stopping for classifiers; fixed 300 for regressor |
| `max_depth` | 5 | Limits overfitting; sufficient for tabular biological data |
| `learning_rate` | 0.05 | Conservative shrinkage for 500-tree ensembles |
| `subsample` | 0.8 | Row subsampling per tree; reduces variance |
| `colsample_bytree` | 0.8 | Column subsampling per tree; reduces variance and speeds training |
| `tree_method` | `hist` | Histogram-based algorithm; faster for large feature sets |
| `random_state` | 42 | Reproducibility |
| `eval_metric` | `aucpr` (classifiers) / `rmse` (regressor) | AUPRC is the appropriate metric for imbalanced binary classification |
| `early_stopping_rounds` | 30 (classifiers only) | Patience of 30 rounds on AUPRC |

**Class imbalance correction (classifiers only):** `scale_pos_weight = n_neg / n_pos` is set for all XGBClassifier models. For the propensity model, the full DB gene set is balanced by the 85:15 undrugged:drugged ratio. For the outcome and PU prior models, the drugged-gene positive:negative ratio (~1:2.7) is used.

**Cross-validation:** `StratifiedKFold(n_splits=5, shuffle=True, random_state=42)` for all models.

---

## 7. Implementation Details

| Component | Version |
|-----------|---------|
| Python | 3.13 (Miniconda) |
| XGBoost | 3.2.0 |
| SHAP | 0.47 |
| scikit-learn | 1.8.0 |
| pandas | 2.x |
| pyarrow | 23.x |
| SQLite | bundled with Python |

**Scripts:**

| Script | Purpose |
|--------|---------|
| `backend/load_real_data.py` | Builds `humanproof.db` from CellxGene, Genebass, LOEUF |
| `data/safety_model.py` | Feature loaders, `get_feature_cols()`, helper utilities |
| `data/export_shap_all_genes.py` | Builds expanded feature matrix for all 17,745 DB genes |
| `data/safety_model_dr.py` | Executes 5-step DR+PU pipeline; outputs `predictions.csv`, `model_final.json` |
| `data/export_shap_dr.py` | Batch SHAP computation for all 17,745 genes; outputs `gene_shap_dr.json` |
| `data/plot_dr_diagnostics.py` | 4 diagnostic figures (propensity, scores, pseudo-outcomes, importance) |

**Feature column discovery:** All feature columns are discovered dynamically at runtime from `df.columns` via `get_feature_cols(df)`. There are no hardcoded lists for cell types or pLoF categories downstream of the DB — adding a new cell type to `CELLXGENE_TISSUE_MAP` or a new phenotype category automatically propagates through all downstream feature engineering and modelling steps.

---

## 8. References

1. **Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J.** (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68.

2. **Elkan, C., & Noto, K.** (2008). Learning classifiers from only positive and unlabeled data. *Proceedings of the 14th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, 213–220.

3. **Kiryo, R., Niu, G., du Plessis, M. C., & Sugiyama, M.** (2017). Positive-unlabeled learning with non-negative risk estimator. *Advances in Neural Information Processing Systems*, 30.

4. **Robins, J. M., Rotnitzky, A., & Zhao, L. P.** (1994). Estimation of regression coefficients when some regressors are not always observed. *Journal of the American Statistical Association*, 89(427), 846–866.

5. **CZ CellxGene Consortium** (2023). CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *bioRxiv*.

6. **Genebass:** Karczewski, K. J., et al. (2022). Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. *Cell Genomics*, 2(9), 100168.

7. **gnomAD LOEUF:** Karczewski, K. J., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*, 581(7809), 434–443.

8. **Open Targets Platform:** Ochoa, D., et al. (2023). The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. *Nucleic Acids Research*, 51(D1), D1353–D1359.

9. **XGBoost:** Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. *Proceedings of the 22nd ACM SIGKDD*, 785–794.

10. **SHAP:** Lundberg, S. M., & Lee, S.-I. (2017). A unified approach to interpreting model predictions. *Advances in Neural Information Processing Systems*, 30.
