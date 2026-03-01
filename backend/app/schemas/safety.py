from pydantic import BaseModel
from typing import Literal

from app.schemas.expression import CellTypeExpression
from app.schemas.plof import PLOFAssociationSchema, GeneDosageSensitivitySchema


class ShapFeature(BaseModel):
    name: str
    label: str
    group: Literal["ot", "genetics", "expression"]
    shap_value: float
    feature_value: float | None


class GeneSHAP(BaseModel):
    gene_symbol: str
    ensembl_id: str
    safety_score: float
    base_value: float
    model: Literal["A", "B", "DR"]
    safety_label: str
    features: list[ShapFeature]


class ExpressionSummaryForCard(BaseModel):
    n_tissues: int
    n_cell_types: int
    n_expressing_tissues: int    # tissues with ≥1 cell type at pct_expressed > 10%
    n_expressing_cell_types: int # cell types with pct_expressed > 10%
    top_tissue: str
    top_cell_type: str
    max_expression: float
    expression_breadth: float
    top_entries: list[CellTypeExpression]


class ScoreDistribution(BaseModel):
    bins: list[float]           # 21 edges defining 20 bins [0, 0.05, ..., 1.0]
    drugged_safety: list[int]   # 20 counts per bin: drugged + safety event
    drugged_no_safety: list[int]
    undrugged: list[int]


class PLOFSummaryForCard(BaseModel):
    n_associations: int
    n_significant: int   # p < 5e-8 (genome-wide significant)
    n_suggestive: int    # 5e-8 ≤ p < 1e-5
    top_phenotype: str | None
    max_severity: str
    organ_systems_affected: list[str]
    top_entries: list[PLOFAssociationSchema]


class SafetyCard(BaseModel):
    gene_symbol: str
    ensembl_id: str
    risk_class: Literal["low", "moderate", "high", "critical"]
    humanproof_score: float | None
    humanproof_model: Literal["A", "B", "DR"] | None
    is_drugged: bool
    clinical_phase: float | None
    has_safety_event: bool | None
    expression_summary: ExpressionSummaryForCard
    plof_summary: PLOFSummaryForCard
    dosage_sensitivity: GeneDosageSensitivitySchema | None


class GeneSearchResult(BaseModel):
    gene_symbol: str
    ensembl_id: str
    n_tissues: int
    n_plof_associations: int
    risk_class: str
