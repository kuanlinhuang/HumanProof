from sqlalchemy import String, Float, Integer, Index
from sqlalchemy.orm import Mapped, mapped_column

from app.core.database import Base


class ExpressionSummary(Base):
    __tablename__ = "expression_summary"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    gene_symbol: Mapped[str] = mapped_column(String(32), index=True)
    ensembl_id: Mapped[str] = mapped_column(String(32))
    cell_type: Mapped[str] = mapped_column(String(128))
    tissue: Mapped[str] = mapped_column(String(64), index=True)
    organ: Mapped[str] = mapped_column(String(64))
    mean_expression: Mapped[float] = mapped_column(Float)
    pct_expressed: Mapped[float] = mapped_column(Float)
    n_cells: Mapped[int] = mapped_column(Integer)

    __table_args__ = (
        Index("idx_expression_gene_cell", "gene_symbol", "cell_type"),
    )
