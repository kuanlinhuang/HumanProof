from pydantic import BaseModel
from typing import Literal


class PLOFAssociationSchema(BaseModel):
    phenotype: str
    phenotype_category: str
    organ_system: str
    p_value: float
    p_value_burden: float | None = None
    p_value_skat: float | None = None
    beta: float
    se: float
    n_carriers: int | None = None
    direction: str


class GenePLOFProfile(BaseModel):
    gene_symbol: str
    n_associations: int
    associations: list[PLOFAssociationSchema]
    top_phenotype: str | None
    max_severity: Literal["benign", "moderate", "severe", "critical"]


class PheWASResponse(BaseModel):
    genes: list[str]
    associations: list[dict]  # [{gene, phenotype, category, p_value, beta, ...}]


class GeneDosageSensitivitySchema(BaseModel):
    gene_symbol: str
    pli_score: float
    loeuf_score: float
    mis_z_score: float
    risk_class: str
