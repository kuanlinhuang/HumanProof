from pydantic import BaseModel


class CellTypeExpression(BaseModel):
    cell_type: str
    tissue: str
    organ: str
    mean_expression: float
    pct_expressed: float
    n_cells: int


class GeneExpressionProfile(BaseModel):
    gene_symbol: str
    ensembl_id: str
    cell_types: list[CellTypeExpression]


class ExpressionHeatmapResponse(BaseModel):
    genes: list[str]
    cell_types: list[str]
    tissues: list[str]
    matrix: list[list[float]]  # genes × cell_types


class TissueInfo(BaseModel):
    tissue: str
    organ: str
    cell_types: list[str]


class TissuesResponse(BaseModel):
    tissues: list[TissueInfo]
