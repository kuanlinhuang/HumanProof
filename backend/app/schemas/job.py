"""Pydantic schemas for prediction job responses."""

from pydantic import BaseModel
from datetime import datetime


class BindingPredictionSchema(BaseModel):
    """A single binding prediction result."""

    gene_symbol: str
    ensembl_id: str
    binding_score: float
    confidence: float
    binding_site: str | None = None
    interaction_type: str
    delta_g: float
    kd_nm: float
    rank: int

    class Config:
        from_attributes = True


class PredictionJobStatus(BaseModel):
    """Response schema for job status queries."""

    job_id: str
    status: str
    sequence_type: str
    sequence_name: str
    created_at: datetime
    completed_at: datetime | None = None
    n_targets_found: int = 0
    error_message: str | None = None
    predictor_used: str = "mock"

    class Config:
        from_attributes = True


class PredictionJobResult(BaseModel):
    """Full job result with binding predictions."""

    job_id: str
    status: str
    sequence_type: str
    sequence_name: str
    created_at: datetime
    completed_at: datetime | None = None
    n_targets_found: int = 0
    predictor_used: str = "mock"
    predictions: list[BindingPredictionSchema] = []


class PredictionJobList(BaseModel):
    """List of prediction jobs."""

    jobs: list[PredictionJobStatus]
    total: int


class BindingProfileSummary(BaseModel):
    """Summary binding profile linking to safety data."""

    gene_symbol: str
    ensembl_id: str
    binding_score: float
    confidence: float
    kd_nm: float
    interaction_type: str
    risk_class: str | None = None
    n_tissues_expressed: int = 0
    n_plof_associations: int = 0
    top_tissue: str | None = None
    top_phenotype: str | None = None


class PipelineResult(BaseModel):
    """Complete pipeline result linking binding to safety."""

    job_id: str
    sequence_name: str
    sequence_type: str
    n_targets: int
    binding_profiles: list[BindingProfileSummary]
    completed_at: datetime | None = None
