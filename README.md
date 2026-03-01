# HumanProof

**In Silico Drug Target Safety Platform**

> Given a biologic sequence, which human proteins does it engage, in which cells do those engagements matter, and what does human genetics already tell us happens when those genes are disrupted?

HumanProof reframes safety prediction from "simulate an animal" to "interrogate human biology directly" -- a fundamentally stronger epistemic position for regulatory purposes.

---

## Overview

HumanProof is a full-stack web application that predicts the safety profile of biologic therapeutics (antibodies, nanobodies, peptides) by integrating:

- **Binding prediction** against the human proteome
- **Cell-type expression mapping** across 14 tissues and 60 cell types (real CZ CellxGene census data)
- **Loss-of-function genetics** (pLOF phenotype associations from Genebass UK Biobank burden tests)
- **Gene dosage sensitivity** (LOEUF constraint scores from gnomAD via grr.iossifovlab.com)
- **Target prioritisation features** (Open Targets Platform 25.12: druggability, mouse KO phenotypes, genetic constraint, tissue specificity)
- **ML safety prediction** (XGBoost model trained on observed drug safety events as weak labels, with SHAP interpretability)

The platform uses **real biological data** for all **~19,000 protein-coding genes** across four evidence layers. The binding predictor is currently a mock implementation; the architecture supports a seamless swap to GPU-accelerated structure prediction (Boltz-2).

---

## Quick Start

### Prerequisites

- Python >= 3.11
- Node.js >= 18
- npm >= 9

### 1. Clone and set up the backend

```bash
cd backend
python -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e .
pip install greenlet          # required for async SQLAlchemy
```

### 2. Download required data files

The database loader requires three data files (already present in `data/`):

| File | Description | Source |
|------|-------------|--------|
| `data/genebass_pLoF_filtered.pkl` | Genebass UK Biobank pLoF burden-test results | Provided |
| `data/cellxgene/` | CZ CellxGene census expression aggregations | Provided |
| `data/LOEUF_scores.csv.gz` | LOEUF constraint scores | [grr.iossifovlab.com](https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/) |
| `data/opentargets_target_prioritisation.parquet` | Open Targets Platform 25.12 target prioritisation | [ftp.ebi.ac.uk](http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/target_prioritisation/) |

To re-download LOEUF scores:
```bash
curl -L "https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/LOEUF_scores.csv.gz" \
     -o data/LOEUF_scores.csv.gz
```

To re-download Open Targets target prioritisation:
```bash
curl -L "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/target_prioritisation/part-00000-cad91f39-c3ab-4d9a-9a62-407309b45590-c000.snappy.parquet" \
     -o data/opentargets_target_prioritisation.parquet
```

### 3. Build the database from real data

```bash
cd backend
python load_real_data.py
```

This creates `humanproof.db` (SQLite) with **real data** for all ~19,000 protein-coding genes:
- **~1.06M expression records** -- log1p mean expression + fraction expressing from CZ CellxGene census across 14 tissues and 60 cell types (17,745 genes)
- **~3.6M pLoF associations** -- Genebass UK Biobank burden-test results (18,092 genes, top 200 per gene)
- **~19,155 LOEUF scores** -- gene constraint from gnomAD

### 4. Start the backend

```bash
cd backend
uvicorn app.main:app --reload --port 8000
```

The API is now running at `http://localhost:8000`. Verify with:

```bash
curl http://localhost:8000/health
# {"status":"ok","app":"HumanProof"}
```

### 5. Set up and start the frontend

```bash
cd frontend
npm install
npm run dev
```

The app is now running at `http://localhost:3000`.

> **Note:** The Next.js dev server proxies all `/api/*` requests to the backend at `localhost:8000` via rewrites configured in `next.config.ts`. No CORS issues to worry about in development.

---

## Pipeline Architecture

HumanProof implements a 6-module safety assessment pipeline:

```
Biologic Sequence
       |
       v
+------------------------------------------+
| Module 1: Structure Prediction           |  <-- Phase 3
| Predict 3D structure, screen against     |
| ~5K human surface/secreted proteins      |
+------------------------------------------+
       |
       v
+------------------------------------------+
| Module 2: Off-Target Ranking             |  <-- Phase 3
| Filter to actionable off-target          |
| candidates via confidence thresholds     |
+------------------------------------------+
       |
       v
+------------------------------------------+
| Module 3: Cell-Type Expression Map       |  <-- Active (MVP)
| Map targets to cell types where          |
| expressed; binding x expression =        |
| vulnerability                            |
+------------------------------------------+
       |
       v
+------------------------------------------+
| Module 4: pLOF Association Lookup        |  <-- Active (MVP)
| Retrieve loss-of-function phenotypes     |
| as predicted consequences of             |
| pharmacologic inhibition                 |
+------------------------------------------+
       |
       v
+------------------------------------------+
| Module 5: Gene Network Propagation       |  <-- Phase 2
| Propagate perturbations through          |
| cell-type gene regulatory networks       |
+------------------------------------------+
       |
       v
+------------------------------------------+
| Module 6: Integrated Risk Scoring        |  <-- Phase 2
| Aggregate into Human Safety Risk Score   |
| (HSRS) for regulatory submissions        |
+------------------------------------------+
```

**Current MVP (Phase 1)** implements Modules 1 (mock predictor), 3, and 4 with demo data.

---

## Features

### Binding Prediction Pipeline

Submit biologic sequences (antibody heavy/light chains, nanobody VHH domains, or peptides) for proteome-wide binding prediction. The system:

1. Validates the amino acid sequence
2. Runs async binding prediction against all database targets
3. Ranks results by predicted binding affinity (Kd, binding score)
4. Enriches each hit with expression breadth, pLOF count, and risk classification
5. Flags critical safety signals (targets with high constraint scores)

### Gene Safety Cards

For any gene, view an integrated safety card with:

- **Risk Gauge** -- semicircular gauge showing overall risk class (low/moderate/high/critical) based on LOEUF score
- **Expression Heatmap** -- D3.js visualization of expression across tissues and cell types
- **PheWAS Plot** -- Manhattan-style phenome-wide association scatter plot
- **pLOF Association Table** -- sortable table of all loss-of-function phenotype associations
- **SHAP Interpretability** -- per-feature SHAP values from the DR+PU safety model explaining each gene's risk score

### Executive Safety Dashboard

Overview of all protein-coding genes with risk distribution breakdown and a sortable gene table.

### Expression Atlas Explorer

Browse expression patterns for any gene across tissues, with heatmap visualization and data tables.

### pLOF Association Explorer

Explore phenome-wide associations for any gene with interactive PheWAS plots.

---

## Tech Stack

| Layer         | Technology                                                  |
|---------------|-------------------------------------------------------------|
| **Backend**   | Python 3.11+, FastAPI, SQLAlchemy 2.0 (async), Pydantic v2 |
| **Database**  | SQLite via aiosqlite (swappable to PostgreSQL)              |
| **Frontend**  | Next.js 16 (App Router), React 19, TypeScript 5            |
| **Styling**   | Tailwind CSS v4, shadcn/ui components, Radix UI primitives  |
| **Charts**    | D3.js v7 (custom scientific visualizations)                 |
| **Icons**     | Lucide React                                                |

---

## Project Structure

```
HumanProof/
├── data/
│   ├── cellxgene/                          # CZ CellxGene census aggregations (real)
│   │   ├── celltype_metadata.csv
│   │   ├── gene_metadata.csv
│   │   ├── celltype_log1p_mean_expression.csv.gz
│   │   ├── celltype_mean_expression.csv.gz
│   │   └── celltype_fraction_expressing.csv.gz
│   ├── genebass_pLoF_filtered.pkl          # Genebass UK Biobank pLoF results (real)
│   ├── LOEUF_scores.csv.gz                 # LOEUF constraint scores (real)
│   ├── opentargets_target_prioritisation.parquet  # Open Targets Platform 25.12
│   ├── safety_model_dr.py                  # DR+PU safety model (5-step pipeline)
│   ├── export_shap_dr.py                   # Export per-gene SHAP for all 17,745 genes
│   ├── export_shap_all_genes.py            # Expanded feature matrix builder
│   ├── plot_dr_diagnostics.py              # Diagnostic figures for DR model
│   └── safety_model_output/                # Model artefacts (gitignored)
│       └── dr/                             # DR+PU model outputs
│           ├── model_final.json            # Final DR XGBRegressor (Step 5)
│           ├── model_outcome.json          # Outcome model (Step 3)
│           ├── predictions.csv             # DR safety scores for all 17,745 genes
│           ├── gene_shap_dr.json           # Per-gene SHAP values (~455 MB)
│           └── figures/                    # Diagnostic plots
│               ├── 01_propensity_overlap.png   # P(S=1|X) distributions
│               ├── 02_score_distribution.png   # DR score by safety label
│               ├── 03_pseudo_outcome.png        # AIPW pseudo-outcomes
│               └── 04_feature_importance.png   # Top 20 features by mean |SHAP|
│
├── backend/
│   ├── pyproject.toml              # Python dependencies
│   ├── load_real_data.py           # Real data loader (CellxGene + LOEUF + Genebass)
│   ├── generate_demo_data.py       # Legacy synthetic data generator
│   └── app/
│       ├── main.py                 # FastAPI app + lifespan
│       ├── config.py               # Settings (env: HUMANPROOF_*)
│       ├── core/
│       │   └── database.py         # Async SQLAlchemy engine, session, Base
│       ├── models/
│       │   ├── expression.py       # ExpressionSummary table
│       │   ├── plof.py             # PLOFAssociation + GeneDosageSensitivity
│       │   └── job.py              # PredictionJob + BindingPrediction
│       ├── schemas/
│       │   ├── expression.py       # Expression response models
│       │   ├── plof.py             # pLOF response models
│       │   ├── safety.py           # SafetyCard, GeneSearchResult
│       │   ├── sequence.py         # SequenceSubmission + validation
│       │   └── job.py              # Job status, results, pipeline output
│       ├── services/
│       │   ├── sequence_validator.py  # Input validation logic
│       │   └── binding_service.py     # Abstract predictor + mock impl
│       └── api/v1/
│           ├── router.py           # Aggregates all route modules
│           ├── expression.py       # GET /expression/*
│           ├── plof.py             # GET /plof/*
│           ├── targets.py          # GET /targets/* (search, safety card)
│           └── jobs.py             # POST+GET /pipeline/* (binding pred)
│
├── frontend/
│   ├── package.json
│   ├── next.config.ts              # API proxy rewrites
│   └── src/
│       ├── types/api.ts            # All TypeScript interfaces
│       ├── lib/
│       │   ├── api-client.ts       # Typed fetch wrapper + all endpoints
│       │   └── hooks/use-debounce.ts
│       ├── components/
│       │   ├── layout/
│       │   │   ├── Sidebar.tsx     # Collapsible navigation
│       │   │   └── Header.tsx      # Dynamic page title + badges
│       │   ├── visualizations/
│       │   │   ├── ExpressionHeatmap.tsx  # D3 horizontal bar heatmap
│       │   │   ├── PheWASPlot.tsx         # D3 Manhattan scatter
│       │   │   ├── RiskGauge.tsx          # SVG semicircular gauge
│       │   │   └── BindingChart.tsx       # D3 binding affinity bars
│       │   ├── pipeline/
│       │   │   ├── SequenceInput.tsx      # Multi-type sequence form
│       │   │   └── JobProgress.tsx        # 3-step progress stepper
│       │   └── ui/                 # 13 shadcn/ui primitives
│       └── app/
│           ├── layout.tsx          # Root layout (sidebar + header)
│           ├── page.tsx            # Home: gene search
│           ├── dashboard/          # Executive safety dashboard
│           ├── pipeline/           # Binding prediction submission
│           ├── pipeline/[jobId]/   # Prediction results + viz
│           ├── targets/[geneId]/   # Gene safety card
│           ├── explore/expression/ # Expression atlas explorer
│           ├── explore/plof/       # pLOF association explorer
│           ├── methodology/        # DR+PU pipeline documentation
│           └── about/              # Platform description
│
└── data/                           # Generated demo JSON (gitignored)
```

---

## API Reference

All endpoints are prefixed with `/api/v1`.

### Expression Module

| Method | Path                         | Description                          |
|--------|------------------------------|--------------------------------------|
| GET    | `/expression/{gene}`         | Expression profile for a gene        |
| GET    | `/expression/heatmap/data`   | Heatmap matrix for multiple genes    |
| GET    | `/expression/tissues/list`   | All tissues with their cell types    |

### pLOF Module

| Method | Path                        | Description                          |
|--------|-----------------------------|--------------------------------------|
| GET    | `/plof/{gene}`              | pLOF association profile for a gene  |
| GET    | `/plof/phewas/data`         | PheWAS data for multiple genes       |
| GET    | `/plof/dosage/{gene}`       | Dosage sensitivity (pLI, LOEUF)      |

### Targets / Safety Module

| Method | Path                              | Description                      |
|--------|-----------------------------------|----------------------------------|
| GET    | `/targets/search?q=`             | Autocomplete gene search         |
| GET    | `/targets/{gene}/safety-card`     | Integrated safety card           |

### Pipeline / Binding Prediction Module

| Method | Path                                  | Description                          |
|--------|---------------------------------------|--------------------------------------|
| POST   | `/pipeline/validate`                  | Validate a biologic sequence         |
| POST   | `/pipeline/jobs`                      | Submit prediction job (async)        |
| GET    | `/pipeline/jobs`                      | List recent prediction jobs          |
| GET    | `/pipeline/jobs/{id}`                 | Poll job status                      |
| GET    | `/pipeline/jobs/{id}/binding`         | Raw binding predictions              |
| GET    | `/pipeline/jobs/{id}/pipeline`        | Enriched results + safety data       |

### Health Check

| Method | Path       | Description          |
|--------|------------|----------------------|
| GET    | `/health`  | Application health   |

---

## Environment Variables

All backend settings use the `HUMANPROOF_` prefix:

| Variable                   | Default                                  | Description              |
|----------------------------|------------------------------------------|--------------------------|
| `HUMANPROOF_APP_NAME`      | `HumanProof`                             | Application name         |
| `HUMANPROOF_DATABASE_URL`  | `sqlite+aiosqlite:///./humanproof.db`    | Database connection URL  |
| `HUMANPROOF_CORS_ORIGINS`  | `["http://localhost:3000"]`              | Allowed CORS origins     |
| `HUMANPROOF_DATA_DIR`      | `<project>/data`                         | Data directory path      |

Frontend environment (`.env.local`):

| Variable                | Default                    | Description              |
|-------------------------|----------------------------|--------------------------|
| `NEXT_PUBLIC_API_URL`   | `http://localhost:8000`    | Backend API URL          |

---

## Real Data Sources

The database is populated from three real biological data sources via `load_real_data.py`:

### Expression: CZ CellxGene Census
- **Source**: CZ CellxGene single-cell census, aggregated per cell type across all normal human tissues
- **Metric**: log1p mean expression (log1p of mean raw UMI counts per cell type) and fraction of cells expressing each gene
- **Coverage**: 60 biologically curated cell types across 14 tissues (Brain, Heart, Liver, Lung, Kidney, Intestine, Bone marrow, Lymph node, Spleen, Skin, Pancreas, Muscle, Adrenal, Breast)
- **Records**: ~1.06M (17,745 genes × 60 cell types)
- **Files**: `data/cellxgene/celltype_log1p_mean_expression.csv.gz`, `celltype_fraction_expressing.csv.gz`

### pLoF Associations: Genebass UK Biobank
- **Source**: Genebass exome-wide burden test results from UK Biobank
- **Content**: Phenome-wide association results for protein-truncating variants (pLoF), covering 19 phenotype categories (metabolic, cardiovascular, neurological, renal, hepatic, hematologic, respiratory, musculoskeletal, cancer, and more)
- **Records**: ~3.6M associations across 18,092 genes (top 200 hits per gene)
- **File**: `data/genebass_pLoF_filtered.pkl`

### Gene Constraint: LOEUF (gnomAD)
- **Source**: Loss-of-Function Observed/Expected Upper bound Fraction from gnomAD, distributed via [grr.iossifovlab.com](https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/)
- **Metric**: LOEUF score (0.03–2.0); lower = more intolerant to loss-of-function
- **Coverage**: ~19,155 protein-coding genes
- **Risk classification** (gnomAD-aligned thresholds):
  - `critical`: LOEUF < 0.35
  - `high`: LOEUF 0.35–0.70
  - `moderate`: LOEUF 0.70–1.00
  - `low`: LOEUF ≥ 1.00
- **Note**: pLI scores are not included from this source; pLI fields are set to 0.0
- **File**: `data/LOEUF_scores.csv.gz`

### Target Prioritisation: Open Targets Platform 25.12
- **Source**: [Open Targets Platform](https://platform.opentargets.org/) release 25.12, `target_prioritisation` dataset
- **Coverage**: 78,725 targets × 17 curated features
- **Key features**:

  | Feature | Description |
  |---------|-------------|
  | `hasSafetyEvent` | Known drug safety event (weak label for ML) |
  | `mouseKOScore` | Severity of knockout phenotype in mice |
  | `geneticConstraint` | Tolerance to genetic variation |
  | `isInMembrane` / `isSecreted` | Protein localisation |
  | `hasPocket` / `hasLigand` / `hasSmallMoleculeBinder` | Druggability indicators |
  | `paralogMaxIdentityPercentage` | Functional redundancy via paralogues |
  | `mouseOrthologMaxIdentityPercentage` | Cross-species target conservation |
  | `tissueSpecificity` / `tissueDistribution` | Expression breadth |
  | `maxClinicalTrialPhase` | Highest drug development stage reached |

- **File**: `data/opentargets_target_prioritisation.parquet` (Snappy-compressed Parquet)

### ML Safety Model — DR+PU (Doubly-Robust + Positive-Unlabeled)
- **Script**: `data/safety_model_dr.py`
- **Algorithm**: 5-step doubly-robust pipeline with XGBoost and SHAP interpretability
- **Features** (191 total, all discovered dynamically — no hardcoded filters):
  - Open Targets platform features (14)
  - LOEUF constraint score (1)
  - pLoF Genebass aggregates: 4 global + 38 per-category (19 categories × min p + max |β|) (42)
  - Cell-type expression: 60 cell types × mean + % expressing (120)
  - Organ-level max expression: 14 organs (14)

The pipeline addresses two structural challenges in drug safety prediction:

1. **Selection bias** — genes in clinical trials are a non-random, druggability-enriched subset. Corrected via inverse propensity weighting (IPW).
2. **Positive-unlabeled (PU) structure** — undrugged genes are unlabeled, not confirmed safe. Corrected via Elkan–Noto prior estimation.

**Five-step pipeline** (training population: 17,745 DB genes, 1,506 drugged, 409 positives):
1. **Propensity model** P(S=1|X) — 5-fold CV XGBoost, AUROC ~0.91; generates stabilized IPW weights
2. **PU prior** π_p ≈ 9% — Elkan–Noto method estimates fraction of all genes with true safety liability
3. **Outcome model** P(Y=1|X) — trained on drugged genes with IPW+PU weights; OOF AUROC 0.691, AUPRC 0.426
4. **AIPW pseudo-outcomes** Ỹ — doubly-robust correction combining both models
5. **Final DR model** m_DR(X) — XGBRegressor on all 17,745 genes; deployed predictor for the UI

To run:
```bash
# Requires conda Python with xgboost, shap, scikit-learn, pyarrow
~/miniconda3/bin/python data/safety_model_dr.py

# Export per-gene SHAP for all 17,745 genes (~455 MB output)
~/miniconda3/bin/python data/export_shap_dr.py

# Generate diagnostic figures
~/miniconda3/bin/python data/plot_dr_diagnostics.py
```

**Diagnostic figures** (saved to `data/safety_model_output/dr/figures/`):

| File | Description |
|------|-------------|
| `01_propensity_overlap.png` | P(S=1\|X) distribution for drugged vs undrugged genes |
| `02_score_distribution.png` | DR safety score by label (positive / drugged-safe / undrugged) |
| `03_pseudo_outcome.png` | AIPW pseudo-outcomes Ỹ + DR vs naive m̂ scatter |
| `04_feature_importance.png` | Top 20 features by mean \|SHAP\| (DR model) |

---

## Binding Prediction Architecture

The binding prediction module uses a pluggable predictor pattern:

```
BindingPredictor (abstract)
    |
    +-- MockBindingPredictor   (demo: deterministic, sequence-hash seeded)
    |
    +-- Boltz2Predictor        (future: GPU, pip install boltz)
```

**MockBindingPredictor** generates realistic predictions by:
- Using the sequence MD5 hash as a deterministic seed
- Calculating base affinity from sequence properties (hydrophobicity, charge)
- Boosting scores for known druggable targets (EGFR, ERBB2, PDCD1, etc.)
- Computing thermodynamic properties (delta-G, Kd) from binding scores
- Annotating known binding domains and interaction types

To swap to Boltz-2, implement the `BindingPredictor` interface and update `get_predictor()` in `binding_service.py`.

---

## Development

### Backend development

```bash
cd backend
source .venv/bin/activate
uvicorn app.main:app --reload --port 8000
```

The `--reload` flag enables hot-reloading on file changes. API docs are at `http://localhost:8000/docs` (Swagger) and `http://localhost:8000/redoc` (ReDoc).

### Frontend development

```bash
cd frontend
npm run dev
```

Next.js runs on port 3000 with hot module replacement. The proxy to the backend is configured in `next.config.ts`.

### Running tests

```bash
cd backend
pip install -e ".[dev]"
pytest
```

### Building for production

```bash
cd frontend
npm run build
npm run start
```

---

## Roadmap

- [x] **Real expression data**: CZ CellxGene census per-cell-type aggregations (60 cell types, 14 tissues)
- [x] **Real pLoF data**: Genebass UK Biobank exome-wide burden test results
- [x] **Real constraint data**: LOEUF scores from gnomAD
- [x] **Open Targets integration**: target prioritisation features (78K targets, 17 features)
- [x] **ML safety model**: DR+PU doubly-robust model across all ~17,745 protein-coding genes
- [x] **Safety score in UI**: DR safety scores surfaced in gene safety cards and dashboard
- [x] **SHAP interpretability**: Per-feature SHAP explanations in gene safety cards (all 191 features)
- [x] **Methodology page**: `/methodology` page documenting the 5-step DR+PU pipeline
- [ ] **pLI scores**: Add gnomAD pLI alongside LOEUF (full constraint metrics)
- [ ] **Phase 2**: Gene network propagation (Module 5), integrated HSRS scoring (Module 6)
- [ ] **Phase 3**: Boltz-2 structure prediction (Module 1), off-target ranking (Module 2)
- [ ] **PostgreSQL**: Swap SQLite for production database
- [ ] **Authentication**: User accounts and saved analyses
- [ ] **PDF export**: Regulatory-ready safety reports
- [ ] **Batch mode**: Multi-sequence comparison analysis

---

## License

This project is proprietary. All rights reserved.
