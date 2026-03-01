# Contributing to HumanProof

This guide covers the project architecture, conventions, and workflows for contributing to HumanProof.

---

## Table of Contents

- [Architecture Overview](#architecture-overview)
- [Real Data Sources](#real-data-sources)
- [Backend Architecture](#backend-architecture)
- [Frontend Architecture](#frontend-architecture)
- [Data Flow](#data-flow)
- [Adding a New Module](#adding-a-new-module)
- [Code Conventions](#code-conventions)
- [Common Tasks](#common-tasks)
- [Python Environments](#python-environments)

---

## Architecture Overview

HumanProof is a monorepo with two independently deployable services:

```
                    +-----------+
                    |  Browser  |
                    +-----+-----+
                          |
                          | HTTP (port 3000)
                          v
                  +-------+--------+
                  |   Next.js 16   |
                  |  (App Router)  |
                  |   TypeScript   |
                  +-------+--------+
                          |
                          | /api/* rewrite proxy
                          | (next.config.ts)
                          v
                  +-------+--------+
                  |   FastAPI      |
                  |  (async)       |
                  |  Python 3.11+  |
                  +-------+--------+
                          |
                          | SQLAlchemy 2.0 async
                          v
                  +-------+--------+
                  |    SQLite      |
                  |  (aiosqlite)   |
                  +----------------+
```

### Why this architecture?

- **Next.js App Router** gives us server components, file-based routing, and zero-config TypeScript
- **FastAPI** provides async request handling, automatic OpenAPI docs, and Pydantic validation
- **SQLAlchemy 2.0 async** allows non-blocking database operations
- **SQLite** keeps the MVP zero-infrastructure; the async engine is swappable to PostgreSQL
- **Next.js rewrites** proxy `/api/*` to the backend, eliminating CORS issues and keeping a single origin for the browser

---

## Real Data Sources

HumanProof uses real biological data across all three evidence modules. The database is built by `backend/load_real_data.py`.

### Expression: CZ CellxGene Census

- **What**: Aggregated single-cell RNA-seq expression from the CZ CellxGene census
- **Metrics stored**:
  - `mean_expression` — log1p of mean UMI count per cell type (from `celltype_log1p_mean_expression.csv.gz`)
  - `pct_expressed` — fraction of cells with count > 0 (from `celltype_fraction_expressing.csv.gz`)
  - `n_cells` — total cells of that type across all tissues (from `celltype_metadata.csv`)
- **Cell type selection**: 60 cell types across 14 tissues, curated via `CELLXGENE_TISSUE_MAP` in `load_real_data.py`
- **Data directory**: `data/cellxgene/`

The expression matrices are wide-format gzip-compressed CSVs: rows = cell types, columns = gene symbols. The loader reads all columns and loads all ~17,745 genes covered by the DB.

### pLoF: Genebass UK Biobank

- **What**: Exome-wide burden test results from the UK Biobank, covering all protein-coding genes
- **Metrics stored**: combined p-value, burden test p-value, SKAT p-value, effect size (beta), SE, direction
- **File**: `data/genebass_pLoF_filtered.pkl` (pandas DataFrame)
- **Phenotype classification**: rule-based keyword mapping in `classify_phenotype()` (19 categories including metabolic, cardiovascular, neurological, renal, hepatic, hematologic, respiratory, musculoskeletal, cancer, procedural, medication, gastrointestinal, and more)

### Gene Constraint: LOEUF

- **What**: Loss-of-Function Observed/Expected Upper bound Fraction from gnomAD
- **Source**: [https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/](https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/)
- **File**: `data/LOEUF_scores.csv.gz`
- **Risk classification thresholds** (LOEUF-based, gnomAD-aligned):

  | Risk class | LOEUF range | Interpretation |
  |------------|-------------|----------------|
  | critical   | < 0.35      | Highly intolerant to LOF |
  | high       | 0.35–0.70   | Moderately intolerant |
  | moderate   | 0.70–1.00   | Somewhat tolerant |
  | low        | ≥ 1.00      | Tolerant |

- **Note**: pLI scores are not available from this source and are stored as `0.0`. The `classify_severity()` function in `app/api/v1/plof.py` uses LOEUF directly for critical classification.

### Target Prioritisation: Open Targets Platform 25.12

- **What**: Curated gene-level features from the Open Targets Platform covering druggability, genetic constraint, mouse phenotypes, and drug development history
- **Source**: `http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/target_prioritisation/`
- **File**: `data/opentargets_target_prioritisation.parquet` (Snappy Parquet, ~908 KB)
- **Coverage**: 78,725 human gene targets × 17 features
- **Key fields**:

  | Field | Type | Notes |
  |-------|------|-------|
  | `targetId` | str | Ensembl gene ID (join key) |
  | `hasSafetyEvent` | float | `-1.0` = confirmed safety event; `NaN` = not recorded |
  | `maxClinicalTrialPhase` | float | 0.25/0.50/0.75/1.0 = Phase I–IV; `NaN` = not drugged |
  | `mouseKOScore` | float | Aggregate phenotypic severity from mouse knockout |
  | `geneticConstraint` | float | Target-level genetic constraint (OT-normalised) |
  | `isInMembrane` / `isSecreted` | float | Protein localisation |
  | `hasPocket` / `hasLigand` / `hasSmallMoleculeBinder` | float | Druggability flags |
  | `paralogMaxIdentityPercentage` | float | Functional redundancy |
  | `mouseOrthologMaxIdentityPercentage` | float | Conservation in mouse |
  | `tissueSpecificity` / `tissueDistribution` | float | Expression breadth |

- **Encoding**: binary fields use `1.0` = true, `0.0` = false, `NaN` = not assessed. `hasSafetyEvent` uses `-1.0` = true (Open Targets convention).

### ML Safety Model — DR+PU (Doubly-Robust + Positive-Unlabeled)

- **Script**: `data/safety_model_dr.py`
- **Purpose**: Predict genome-wide target safety liability via a 5-step doubly-robust pipeline that simultaneously corrects for selection bias (IPW) and positive-unlabeled label noise (Elkan–Noto)
- **Algorithm**: XGBoost (`xgboost>=3.2`) with SHAP feature attribution (`shap>=0.47`)
- **Feature groups** (191 total — all discovered dynamically from df.columns, no hardcoded filters):
  1. Open Targets platform features (14)
  2. LOEUF constraint score (1)
  3. pLoF Genebass global aggregates: min p-value, n significant, max |β|, n phenotypes (4)
  4. pLoF per-category: 19 phenotype categories × min p + max |β| (38)
  5. Cell-type expression: 60 cell types × mean expression + % expressing (120)
  6. Organ-level max expression: 14 organs (14)

**Five-step DR+PU pipeline** (training population: 17,745 DB genes, 1,506 drugged, 409 positives):

| Step | Name | Description | Key output |
|------|------|-------------|------------|
| 1 | Propensity model | 5-fold CV XGBoost; predicts P(S=1\|X) | Stabilized IPW weights; AUROC ~0.91 |
| 2 | PU prior | Elkan–Noto method estimates π_p = P(Y=1) | π_p ≈ 9% |
| 3 | Outcome model | IPW+PU weighted XGBoost on drugged genes; 5-fold OOF | OOF AUROC 0.691, AUPRC 0.426 |
| 4 | AIPW pseudo-outcomes | Ỹ = m̂(X) + (S/π̂)(Y − m̂_oof) | Doubly-robust bias correction |
| 5 | Final DR model | XGBRegressor on all 17,745 genes with pseudo-outcome targets | Deployed predictor m_DR(X) ∈ [0,1] |

**Outputs** (in `data/safety_model_output/dr/`):
- `model_final.json` — final DR XGBRegressor (Step 5)
- `model_outcome.json` — outcome model (Step 3)
- `predictions.csv` — DR safety scores for all 17,745 genes
- `gene_shap_dr.json` — per-gene SHAP values (~455 MB; generated by `export_shap_dr.py`)

---

## Backend Architecture

### Directory Layout

```
backend/
├── pyproject.toml                  # Dependencies and project metadata
├── load_real_data.py               # Real data loader (CellxGene + LOEUF + Genebass)
├── generate_demo_data.py           # Legacy synthetic data generator (kept for reference)
│
└── app/
    ├── main.py                     # FastAPI app factory + lifespan events
    ├── config.py                   # Pydantic Settings (env: HUMANPROOF_*)
    │
    ├── core/
    │   └── database.py             # Engine, session factory, Base class
    │
    ├── models/                     # SQLAlchemy ORM models (database tables)
    │   ├── __init__.py             # Re-exports all models (required for metadata)
    │   ├── expression.py           # ExpressionSummary
    │   ├── plof.py                 # PLOFAssociation, GeneDosageSensitivity
    │   └── job.py                  # PredictionJob, BindingPrediction
    │
    ├── schemas/                    # Pydantic v2 request/response models
    │   ├── expression.py           # CellTypeExpression, GeneExpressionProfile, ...
    │   ├── plof.py                 # PLOFAssociationSchema, GenePLOFProfile, ...
    │   ├── safety.py               # SafetyCard, GeneSearchResult
    │   ├── sequence.py             # SequenceSubmission, SequenceValidationResult
    │   └── job.py                  # PredictionJobStatus, PipelineResult, ...
    │
    ├── services/                   # Business logic (not tied to HTTP)
    │   ├── sequence_validator.py   # Sequence validation (length, AA chars, type)
    │   └── binding_service.py      # Abstract predictor, mock impl, job runner
    │
    └── api/v1/                     # HTTP route handlers
        ├── router.py               # Aggregates all sub-routers
        ├── expression.py           # /api/v1/expression/*
        ├── plof.py                 # /api/v1/plof/*
        ├── targets.py              # /api/v1/targets/*
        └── jobs.py                 # /api/v1/pipeline/*
```

### Key Design Decisions

#### Models vs. Schemas

- **Models** (`app/models/`) are SQLAlchemy ORM classes that map to database tables. They use `Mapped[]` type annotations (SQLAlchemy 2.0 style).
- **Schemas** (`app/schemas/`) are Pydantic v2 models for API request/response validation. They are separate from ORM models to allow different shapes for input, output, and database storage.

#### The Service Layer

Business logic lives in `app/services/`, not in route handlers. Route handlers (`app/api/v1/`) only handle HTTP concerns (parsing parameters, returning responses, error codes). This makes the logic testable without HTTP.

**Example: Binding Prediction Flow**

```
POST /api/v1/pipeline/jobs
    |
    v
jobs.py (route handler)
    |-- validate_submission()          <-- services/sequence_validator.py
    |-- submit_prediction_job()        <-- services/binding_service.py
    |-- background_tasks.add_task()    <-- FastAPI BackgroundTasks
    |
    v
run_prediction() runs asynchronously:
    |-- get all gene targets from DB
    |-- predictor.predict()            <-- MockBindingPredictor or Boltz2Predictor
    |-- filter significant binders
    |-- store BindingPrediction records
    |-- mark job as completed
```

#### Abstract Predictor Pattern

The binding service uses an abstract base class:

```python
class BindingPredictor(ABC):
    @abstractmethod
    async def predict(self, sequence, sequence_type, gene_targets) -> list[dict]: ...

    @property
    @abstractmethod
    def name(self) -> str: ...

class MockBindingPredictor(BindingPredictor):
    # Deterministic mock using sequence hash as seed
    ...

# Future:
class Boltz2Predictor(BindingPredictor):
    # Real GPU-accelerated structure prediction
    ...
```

To swap predictors, modify `get_predictor()` in `binding_service.py`.

#### Async Database Sessions

All database operations use `async with` sessions:

```python
async def get_db() -> AsyncSession:
    async with async_session() as session:
        yield session
```

Route handlers receive the session via FastAPI's `Depends(get_db)`.

#### Background Tasks

Long-running predictions use FastAPI's `BackgroundTasks` (not Celery/Redis). The task creates its own database session to avoid lifecycle issues:

```python
async def run_in_background(job_id: str):
    async with async_session() as session:
        await run_prediction(session, job_id)

background_tasks.add_task(run_in_background, job.id)
```

---

## Frontend Architecture

### Directory Layout

```
frontend/src/
├── types/
│   └── api.ts                      # All TypeScript interfaces (mirrors backend schemas)
│
├── lib/
│   ├── api-client.ts               # Typed fetch wrapper + all API methods
│   ├── utils.ts                    # cn() utility for className merging
│   └── hooks/
│       └── use-debounce.ts         # Debounce hook for search inputs
│
├── components/
│   ├── layout/
│   │   ├── Sidebar.tsx             # Collapsible nav with active-route highlighting
│   │   └── Header.tsx              # Dynamic page title + demo badge
│   │
│   ├── visualizations/             # D3.js scientific charts
│   │   ├── ExpressionHeatmap.tsx   # Horizontal bars, Viridis color, expression + pct
│   │   ├── PheWASPlot.tsx          # Manhattan scatter, 12 categories, jittered
│   │   ├── RiskGauge.tsx           # SVG semicircular gauge, animated needle
│   │   └── BindingChart.tsx        # Ranked bars, risk-class colors, confidence dots
│   │
│   ├── pipeline/                   # Binding prediction UI components
│   │   ├── SequenceInput.tsx       # Multi-type form (antibody/nanobody/peptide)
│   │   └── JobProgress.tsx         # 3-step progress stepper with polling
│   │
│   └── ui/                         # shadcn/ui primitives (13 components)
│       ├── button.tsx
│       ├── card.tsx
│       ├── input.tsx
│       ├── badge.tsx
│       ├── table.tsx
│       ├── tabs.tsx
│       ├── tooltip.tsx
│       ├── separator.tsx
│       ├── dialog.tsx
│       ├── popover.tsx
│       ├── scroll-area.tsx
│       ├── skeleton.tsx
│       └── command.tsx
│
└── app/                            # Next.js App Router pages
    ├── layout.tsx                  # Root layout: sidebar + header + main
    ├── page.tsx                    # /          Gene search home
    ├── pipeline/
    │   ├── page.tsx                # /pipeline        Sequence submission
    │   └── [jobId]/
    │       └── page.tsx            # /pipeline/:id    Results + visualizations
    ├── targets/
    │   └── [geneId]/
    │       └── page.tsx            # /targets/:gene   Safety card
    ├── dashboard/
    │   └── page.tsx                # /dashboard       Executive overview
    ├── explore/
    │   ├── expression/
    │   │   └── page.tsx            # /explore/expression  Expression atlas
    │   └── plof/
    │       └── page.tsx            # /explore/plof        pLOF explorer
    ├── methodology/
    │   └── page.tsx                # /methodology     DR+PU pipeline documentation
    └── about/
        └── page.tsx                # /about           Platform description
```

### Key Design Decisions

#### App Router (Next.js 16)

All pages use the App Router with `"use client"` directives for interactive components. Pages are organized by feature:

- `/pipeline/*` -- Binding prediction workflow
- `/targets/*` -- Gene safety cards
- `/explore/*` -- Data exploration modules
- `/dashboard` -- Executive overview

Dynamic routes use `[param]` folder naming (e.g., `[jobId]`, `[geneId]`).

#### API Client Pattern

All API calls go through a single typed client (`lib/api-client.ts`):

```typescript
const API_BASE = "";  // Empty = same-origin (proxied by Next.js)

async function fetchApi<T>(path: string, options?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, options);
  if (!res.ok) throw new Error(`API error ${res.status}`);
  return res.json();
}

export const api = {
  searchGenes: (q) => fetchApi<GeneSearchResult[]>(`/api/v1/targets/search?q=${q}`),
  getSafetyCard: (gene) => fetchApi<SafetyCard>(`/api/v1/targets/${gene}/safety-card`),
  submitPredictionJob: (sub) => postApi<PredictionJobStatus>("/api/v1/pipeline/jobs", sub),
  // ... all other endpoints
};
```

`API_BASE` is empty because Next.js rewrites `/api/*` to the backend -- so all fetches are same-origin.

#### D3.js Visualization Pattern

All D3 charts follow the same pattern:

```typescript
export function MyChart({ data }: Props) {
  const svgRef = useRef<SVGSVGElement>(null);

  useEffect(() => {
    if (!svgRef.current || data.length === 0) return;
    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();           // Clear previous render

    // ... D3 drawing code ...

    // Tooltip: append to body, remove on cleanup
    const tooltip = d3.select("body").append("div")...;
    return () => { tooltip.remove(); };
  }, [data]);                              // Re-render on data change

  return <svg ref={svgRef} className="w-full" />;
}
```

Key conventions:
- Use `viewBox` for responsive SVGs
- Tooltips are appended to `<body>` and cleaned up in the effect's return
- Use `useRef` for the SVG element, never manipulate React-managed DOM

#### Component Organization

Components are organized by domain, not by type:

- `components/layout/` -- App shell (sidebar, header)
- `components/visualizations/` -- D3.js charts (reusable across pages)
- `components/pipeline/` -- Binding prediction-specific UI
- `components/ui/` -- Generic primitives (from shadcn/ui)

#### Type Safety

`types/api.ts` mirrors every backend Pydantic schema as a TypeScript interface. When adding a new backend schema, always add the corresponding TypeScript type.

---

## Data Flow

### Gene Search Flow

```
User types in search box
    |
    v
useDebounce(300ms)
    |
    v
api.searchGenes(query)
    |
    v
GET /api/v1/targets/search?q=EGF
    |
    v
SQL: SELECT DISTINCT gene_symbol FROM expression_summary WHERE gene_symbol LIKE 'EGF%'
    |
    v
GeneSearchResult[] returned to frontend
    |
    v
User clicks gene -> navigate to /targets/EGFR
    |
    v
api.getSafetyCard("EGFR")
    |
    v
GET /api/v1/targets/EGFR/safety-card
    |
    v
Backend joins expression + pLOF + dosage tables
    |
    v
SafetyCard { risk_class, expression_summary, plof_summary, dosage_sensitivity }
    |
    v
Render RiskGauge + ExpressionHeatmap + PheWASPlot
```

### Binding Prediction Flow

```
User fills sequence form -> clicks "Run Binding Prediction"
    |
    v
api.submitPredictionJob(submission)
    |
    v
POST /api/v1/pipeline/jobs
    |
    v
1. Validate sequence (sequence_validator.py)
2. Create PredictionJob record (status: "pending")
3. Add background task -> redirect to /pipeline/{job_id}
    |
    v
Frontend polls: api.getJobStatus(jobId) every 2 seconds
    |
    v
Background task runs:
    1. Mark job "running"
    2. Load all gene targets from DB
    3. predictor.predict(sequence, targets)
    4. Filter binders with score > 0.3
    5. Store BindingPrediction records
    6. Mark job "completed"
    |
    v
Poll detects status == "completed"
    |
    v
api.getPipelineResults(jobId)
    |
    v
GET /api/v1/pipeline/jobs/{id}/pipeline
    |
    v
Backend enriches each binding hit:
    - Count expression tissues
    - Get top tissue
    - Count pLOF associations
    - Get top phenotype
    - Get risk class from dosage sensitivity
    |
    v
PipelineResult { binding_profiles: BindingProfileSummary[] }
    |
    v
Render BindingChart + BindingTable
    |
    v
User clicks gene -> navigate to /targets/{gene} safety card
```

---

## Adding a New Module

### Backend: Adding a new data module

1. **Create model** in `app/models/your_module.py`:

```python
from sqlalchemy import String, Float
from sqlalchemy.orm import Mapped, mapped_column
from app.core.database import Base

class YourModel(Base):
    __tablename__ = "your_table"
    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    gene_symbol: Mapped[str] = mapped_column(String(32), index=True)
    # ... your columns
```

2. **Register model** in `app/models/__init__.py`:

```python
from app.models.your_module import YourModel
__all__ = [..., "YourModel"]
```

3. **Create schema** in `app/schemas/your_module.py`:

```python
from pydantic import BaseModel

class YourSchema(BaseModel):
    gene_symbol: str
    # ... your fields
    class Config:
        from_attributes = True
```

4. **Create routes** in `app/api/v1/your_module.py`:

```python
from fastapi import APIRouter, Depends
from app.core.database import get_db

router = APIRouter(prefix="/your-module", tags=["Your Module"])

@router.get("/{gene}")
async def get_data(gene: str, db = Depends(get_db)):
    # ... query and return
```

5. **Register router** in `app/api/v1/router.py`:

```python
from app.api.v1.your_module import router as your_router
api_router.include_router(your_router)
```

6. **Add demo data** in `generate_demo_data.py` and re-run it.

### Frontend: Adding a new page

1. **Add TypeScript types** in `types/api.ts`
2. **Add API methods** in `lib/api-client.ts`
3. **Create page** at `app/your-module/page.tsx`
4. **Add to sidebar** in `components/layout/Sidebar.tsx` (NAV_ITEMS array)
5. **Add page title** in `components/layout/Header.tsx` (PAGE_TITLES map)

### Adding a new visualization

1. Create component in `components/visualizations/YourChart.tsx`
2. Follow the D3.js pattern (useRef + useEffect + cleanup)
3. Use `viewBox` for responsive sizing
4. Clean up tooltips in effect return function
5. Accept data as props, re-render on data change

---

## Code Conventions

### Backend

- **Python 3.11+** -- use modern syntax (`str | None`, `list[str]`)
- **Async everywhere** -- all DB operations use `await`
- **Pydantic v2** -- use `model_validate()` not `from_orm()`
- **SQLAlchemy 2.0** -- use `Mapped[]` annotations, `select()` not `.query()`
- **No raw SQL** -- always use the ORM query builder
- **Route handler rule**: thin handlers, business logic in services

### Frontend

- **TypeScript strict** -- no `any` types
- **"use client"** -- all interactive components must declare this
- **Tailwind classes** -- use `cn()` utility for conditional classes
- **Import aliases** -- use `@/` prefix (maps to `src/`)
- **D3 in useEffect** -- never mix D3 DOM manipulation with React rendering
- **API client** -- all fetch calls go through `lib/api-client.ts`, never raw `fetch`

### Naming

| Item            | Convention              | Example                        |
|-----------------|-------------------------|--------------------------------|
| Python files    | `snake_case.py`         | `binding_service.py`           |
| Python classes  | `PascalCase`            | `BindingPredictor`             |
| Python funcs    | `snake_case`            | `run_prediction()`             |
| TS components   | `PascalCase.tsx`        | `BindingChart.tsx`             |
| TS interfaces   | `PascalCase`            | `PredictionJobStatus`          |
| API routes      | `kebab-case`            | `/pipeline/jobs`               |
| DB tables       | `snake_case`            | `binding_predictions`          |
| CSS classes     | Tailwind utilities      | `className="flex gap-3"`       |

---

## Common Tasks

### Rebuild the database from real data

```bash
cd backend
python load_real_data.py
```

This deletes the existing `humanproof.db` and rebuilds it from the three real data sources:

| Source | File | Contents |
|--------|------|----------|
| **CZ CellxGene** | `data/cellxgene/celltype_log1p_mean_expression.csv.gz` | log1p mean expression per cell type per gene |
| **CZ CellxGene** | `data/cellxgene/celltype_fraction_expressing.csv.gz` | Fraction of cells expressing each gene |
| **Genebass** | `data/genebass_pLoF_filtered.pkl` | UK Biobank pLoF burden-test associations |
| **gnomAD/LOEUF** | `data/LOEUF_scores.csv.gz` | LOEUF gene constraint scores |

### Run the safety prediction model (DR+PU)

Requires the conda Python environment (see [Python Environments](#python-environments)):

```bash
# Step 1–5: Train the full DR+PU pipeline (~10 min)
~/miniconda3/bin/python data/safety_model_dr.py

# Export per-gene SHAP for all 17,745 genes (~455 MB output, ~20 min)
~/miniconda3/bin/python data/export_shap_dr.py

# Generate 4 diagnostic figures
~/miniconda3/bin/python data/plot_dr_diagnostics.py
```

`safety_model_dr.py` executes the 5-step pipeline:
1. Loads `data/opentargets_target_prioritisation.parquet`
2. Loads LOEUF, pLoF (all 19 categories), and expression (all 60 cell types) from `backend/humanproof.db`
3. Joins all features into a 17,745 × 191 matrix (features discovered dynamically)
4. Runs Steps 1–5 (propensity → PU prior → outcome → AIPW → DR regressor)
5. Writes outputs to `data/safety_model_output/dr/`

`export_shap_dr.py` loads the final DR model and exports per-gene SHAP for all DB genes:
- Reads `model_final.json` from `data/safety_model_output/dr/`
- Computes SHAP with batch TreeExplainer
- Writes `gene_shap_dr.json` (~455 MB)

Diagnostic figures (saved to `data/safety_model_output/dr/figures/`):

| Figure | Content |
|--------|---------|
| `01_propensity_overlap.png` | P(S=1\|X) density for drugged vs undrugged genes |
| `02_score_distribution.png` | DR score histograms by label (log scale) |
| `03_pseudo_outcome.png` | AIPW Ỹ distribution + DR vs naive m̂ scatter |
| `04_feature_importance.png` | Top 20 features by mean \|SHAP\| (color-coded by group) |

To re-download the Open Targets parquet:

```bash
curl -L "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/target_prioritisation/part-00000-cad91f39-c3ab-4d9a-9a62-407309b45590-c000.snappy.parquet" \
     -o data/opentargets_target_prioritisation.parquet
```

#### Re-downloading LOEUF scores

```bash
curl -L "https://grr.iossifovlab.com/gene_properties/gene_scores/LOEUF/LOEUF_scores.csv.gz" \
     -o data/LOEUF_scores.csv.gz
```

#### Dynamic feature engineering

All feature columns are discovered dynamically — there are no hardcoded cell-type or pLoF-category lists. The key functions are:

- `get_feature_cols(df)` in `safety_model.py` — returns all `expr_ct_*`, `expr_organ_*`, and `plof_*` columns present in `df`, plus the static OT and LOEUF features
- `_ct_col(cell_type)` — sanitizes a cell-type label to a valid column name: `"expr_ct_" + re.sub(r"[^a-z0-9]+", "_", ct.lower()).strip("_")`
- `get_readable_label(feat)` in `export_gene_shap.py` — generates human-readable labels for any feature algorithmically (no static dict needed for expression/pLoF features)

### Add a new shadcn/ui component

```bash
cd frontend
npx shadcn@latest add <component-name>
```

### View API documentation

With the backend running:
- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

### Reset the database

```bash
cd backend
rm humanproof.db
python generate_demo_data.py
# Restart the backend
```

### Check TypeScript types

```bash
cd frontend
npx tsc --noEmit
```

### Production build check

```bash
cd frontend
npm run build
```

This runs TypeScript checking and generates optimized output. All pages must compile without errors.

---

## Python Environments

HumanProof uses two Python runtimes for different purposes:

| Runtime | Path | Use |
|---------|------|-----|
| **System Python 3.14** (Homebrew) | `/opt/homebrew/opt/python@3.14/bin/python3` | Backend FastAPI server (`backend/load_real_data.py`, `uvicorn`) |
| **Conda Python 3.13** (Miniconda) | `~/miniconda3/bin/python` | Data science scripts (`data/safety_model_dr.py`, `export_shap_dr.py`, `plot_dr_diagnostics.py`) |

The conda environment has the full ML stack installed:

```
pyarrow >= 23      # Parquet I/O
xgboost >= 3.2     # Gradient boosted trees
scikit-learn >= 1.8
shap               # SHAP feature attribution
lightgbm >= 4.6    # Alternative to XGBoost (available but not used by default)
pandas >= 2        # Data manipulation
```

To install missing ML packages into conda:

```bash
~/miniconda3/bin/pip install xgboost lightgbm shap scikit-learn pyarrow
```

The backend venv (created in step 1 of Quick Start) is used only for the FastAPI server and does not need the ML stack.
