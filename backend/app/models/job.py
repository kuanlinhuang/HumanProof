"""Models for binding prediction jobs and results."""

from sqlalchemy import String, Float, Integer, Text, DateTime, Index
from sqlalchemy.orm import Mapped, mapped_column
from datetime import datetime

from app.core.database import Base


class PredictionJob(Base):
    """Tracks async binding prediction jobs."""

    __tablename__ = "prediction_jobs"

    id: Mapped[str] = mapped_column(String(36), primary_key=True)  # UUID
    status: Mapped[str] = mapped_column(String(16), default="pending")  # pending, running, completed, failed
    sequence_type: Mapped[str] = mapped_column(String(16))  # antibody, nanobody, peptide
    sequence_name: Mapped[str] = mapped_column(String(128), default="")
    sequence: Mapped[str] = mapped_column(Text)
    heavy_chain: Mapped[str | None] = mapped_column(Text, nullable=True)
    light_chain: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime, nullable=True)
    n_targets_found: Mapped[int] = mapped_column(Integer, default=0)
    error_message: Mapped[str | None] = mapped_column(Text, nullable=True)
    predictor_used: Mapped[str] = mapped_column(String(32), default="mock")  # mock, boltz2

    __table_args__ = (
        Index("idx_job_status", "status"),
        Index("idx_job_created", "created_at"),
    )


class BindingPrediction(Base):
    """Individual binding prediction results for a job."""

    __tablename__ = "binding_predictions"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    job_id: Mapped[str] = mapped_column(String(36), index=True)
    gene_symbol: Mapped[str] = mapped_column(String(32), index=True)
    ensembl_id: Mapped[str] = mapped_column(String(32))
    binding_score: Mapped[float] = mapped_column(Float)  # 0-1, higher = stronger binding
    confidence: Mapped[float] = mapped_column(Float)  # 0-1, prediction confidence
    binding_site: Mapped[str | None] = mapped_column(String(128), nullable=True)  # e.g., "Domain: ECD (residues 1-621)"
    interaction_type: Mapped[str] = mapped_column(String(32))  # competitive, allosteric, orthosteric
    delta_g: Mapped[float] = mapped_column(Float)  # kcal/mol binding free energy
    kd_nm: Mapped[float] = mapped_column(Float)  # dissociation constant in nM
    rank: Mapped[int] = mapped_column(Integer)  # rank by binding_score

    __table_args__ = (
        Index("idx_binding_job_gene", "job_id", "gene_symbol"),
        Index("idx_binding_score", "job_id", "binding_score"),
    )
