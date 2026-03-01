from sqlalchemy import String, Float, Integer, Index
from sqlalchemy.orm import Mapped, mapped_column

from app.core.database import Base


class PLOFAssociation(Base):
    __tablename__ = "plof_associations"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    gene_symbol: Mapped[str] = mapped_column(String(32), index=True)
    phenotype: Mapped[str] = mapped_column(String(256))
    phenotype_category: Mapped[str] = mapped_column(String(64))
    organ_system: Mapped[str] = mapped_column(String(64))
    p_value: Mapped[float] = mapped_column(Float)
    p_value_burden: Mapped[float | None] = mapped_column(Float, nullable=True)
    p_value_skat: Mapped[float | None] = mapped_column(Float, nullable=True)
    beta: Mapped[float] = mapped_column(Float)
    se: Mapped[float] = mapped_column(Float)
    n_carriers: Mapped[int | None] = mapped_column(Integer, nullable=True)
    direction: Mapped[str] = mapped_column(String(16))  # "loss" or "gain"

    __table_args__ = (
        Index("idx_plof_gene_pheno", "gene_symbol", "phenotype"),
    )


class GeneDosageSensitivity(Base):
    __tablename__ = "gene_dosage_sensitivity"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    gene_symbol: Mapped[str] = mapped_column(String(32), unique=True, index=True)
    ensembl_id: Mapped[str] = mapped_column(String(32))
    pli_score: Mapped[float] = mapped_column(Float)
    loeuf_score: Mapped[float] = mapped_column(Float)
    mis_z_score: Mapped[float] = mapped_column(Float)
    risk_class: Mapped[str] = mapped_column(String(16))  # low, moderate, high, critical
