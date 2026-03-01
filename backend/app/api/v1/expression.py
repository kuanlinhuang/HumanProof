from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import select, distinct, func
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.database import get_db
from app.models.expression import ExpressionSummary
from app.schemas.expression import (
    GeneExpressionProfile,
    CellTypeExpression,
    ExpressionHeatmapResponse,
    TissuesResponse,
    TissueInfo,
)

router = APIRouter(prefix="/expression", tags=["expression"])


@router.get("/{gene_symbol}", response_model=GeneExpressionProfile)
async def get_gene_expression(gene_symbol: str, db: AsyncSession = Depends(get_db)):
    gene_upper = gene_symbol.upper()
    result = await db.execute(
        select(ExpressionSummary)
        .where(ExpressionSummary.gene_symbol == gene_upper)
        .order_by(ExpressionSummary.mean_expression.desc())
    )
    rows = result.scalars().all()
    if not rows:
        raise HTTPException(status_code=404, detail=f"Gene {gene_upper} not found")

    return GeneExpressionProfile(
        gene_symbol=gene_upper,
        ensembl_id=rows[0].ensembl_id,
        cell_types=[
            CellTypeExpression(
                cell_type=r.cell_type,
                tissue=r.tissue,
                organ=r.organ,
                mean_expression=r.mean_expression,
                pct_expressed=r.pct_expressed,
                n_cells=r.n_cells,
            )
            for r in rows
        ],
    )


@router.get("/heatmap/data", response_model=ExpressionHeatmapResponse)
async def get_expression_heatmap(
    genes: str = Query(..., description="Comma-separated gene symbols"),
    tissues: str | None = Query(None, description="Comma-separated tissues to filter"),
    db: AsyncSession = Depends(get_db),
):
    gene_list = [g.strip().upper() for g in genes.split(",")]

    query = select(ExpressionSummary).where(
        ExpressionSummary.gene_symbol.in_(gene_list)
    )
    if tissues:
        tissue_list = [t.strip() for t in tissues.split(",")]
        query = query.where(ExpressionSummary.tissue.in_(tissue_list))

    result = await db.execute(query)
    rows = result.scalars().all()

    all_cell_types = sorted(set(r.cell_type for r in rows))
    all_tissues = sorted(set(r.tissue for r in rows))
    found_genes = sorted(set(r.gene_symbol for r in rows))

    # Build lookup
    lookup = {}
    for r in rows:
        lookup[(r.gene_symbol, r.cell_type)] = r.mean_expression

    matrix = [
        [lookup.get((gene, ct), 0.0) for ct in all_cell_types]
        for gene in found_genes
    ]

    return ExpressionHeatmapResponse(
        genes=found_genes,
        cell_types=all_cell_types,
        tissues=all_tissues,
        matrix=matrix,
    )


@router.get("/tissues/list", response_model=TissuesResponse)
async def get_tissues(db: AsyncSession = Depends(get_db)):
    result = await db.execute(
        select(
            ExpressionSummary.tissue,
            ExpressionSummary.organ,
            ExpressionSummary.cell_type,
        ).distinct()
    )
    rows = result.all()

    tissue_map: dict[str, TissueInfo] = {}
    for tissue, organ, cell_type in rows:
        if tissue not in tissue_map:
            tissue_map[tissue] = TissueInfo(tissue=tissue, organ=organ, cell_types=[])
        if cell_type not in tissue_map[tissue].cell_types:
            tissue_map[tissue].cell_types.append(cell_type)

    return TissuesResponse(tissues=sorted(tissue_map.values(), key=lambda t: t.tissue))
