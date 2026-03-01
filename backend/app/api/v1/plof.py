from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.database import get_db
from app.models.plof import PLOFAssociation, GeneDosageSensitivity
from app.schemas.plof import (
    GenePLOFProfile,
    PLOFAssociationSchema,
    PheWASResponse,
    GeneDosageSensitivitySchema,
)

router = APIRouter(prefix="/plof", tags=["plof"])


def classify_severity(associations: list, dosage: GeneDosageSensitivity | None) -> str:
    sig_count = sum(1 for a in associations if a.p_value < 5e-8)
    has_severe_pheno = any(
        a.phenotype_category in ("death", "cancer", "cardiovascular", "neurological")
        and a.p_value < 5e-8
        for a in associations
    )
    if dosage and dosage.loeuf_score < 0.2 and has_severe_pheno:
        return "critical"
    if has_severe_pheno or (dosage and dosage.pli_score > 0.9):
        return "severe"
    if sig_count > 3 or (dosage and dosage.pli_score > 0.5):
        return "moderate"
    return "benign"


@router.get("/{gene_symbol}", response_model=GenePLOFProfile)
async def get_gene_plof(gene_symbol: str, db: AsyncSession = Depends(get_db)):
    gene_upper = gene_symbol.upper()

    result = await db.execute(
        select(PLOFAssociation)
        .where(PLOFAssociation.gene_symbol == gene_upper)
        .order_by(PLOFAssociation.p_value)
    )
    associations = result.scalars().all()
    if not associations:
        raise HTTPException(status_code=404, detail=f"No pLOF data for {gene_upper}")

    dosage_result = await db.execute(
        select(GeneDosageSensitivity).where(
            GeneDosageSensitivity.gene_symbol == gene_upper
        )
    )
    dosage = dosage_result.scalar_one_or_none()

    severity = classify_severity(associations, dosage)

    return GenePLOFProfile(
        gene_symbol=gene_upper,
        n_associations=len(associations),
        associations=[
            PLOFAssociationSchema(
                phenotype=a.phenotype,
                phenotype_category=a.phenotype_category,
                organ_system=a.organ_system,
                p_value=a.p_value,
                p_value_burden=a.p_value_burden,
                p_value_skat=a.p_value_skat,
                beta=a.beta,
                se=a.se,
                n_carriers=a.n_carriers,
                direction=a.direction,
            )
            for a in associations
        ],
        top_phenotype=associations[0].phenotype if associations else None,
        max_severity=severity,
    )


@router.get("/phewas/data", response_model=PheWASResponse)
async def get_phewas(
    genes: str = Query(..., description="Comma-separated gene symbols"),
    p_threshold: float = Query(0.05, description="P-value threshold"),
    db: AsyncSession = Depends(get_db),
):
    gene_list = [g.strip().upper() for g in genes.split(",")]

    result = await db.execute(
        select(PLOFAssociation)
        .where(
            PLOFAssociation.gene_symbol.in_(gene_list),
            PLOFAssociation.p_value <= p_threshold,
        )
        .order_by(PLOFAssociation.p_value)
    )
    rows = result.scalars().all()
    found_genes = sorted(set(r.gene_symbol for r in rows))

    return PheWASResponse(
        genes=found_genes,
        associations=[
            {
                "gene": r.gene_symbol,
                "phenotype": r.phenotype,
                "category": r.phenotype_category,
                "organ_system": r.organ_system,
                "p_value": r.p_value,
                "p_value_burden": r.p_value_burden,
                "p_value_skat": r.p_value_skat,
                "beta": r.beta,
                "se": r.se,
                "n_carriers": r.n_carriers,
                "direction": r.direction,
            }
            for r in rows
        ],
    )


@router.get("/dosage/{gene_symbol}", response_model=GeneDosageSensitivitySchema)
async def get_dosage_sensitivity(
    gene_symbol: str, db: AsyncSession = Depends(get_db)
):
    gene_upper = gene_symbol.upper()
    result = await db.execute(
        select(GeneDosageSensitivity).where(
            GeneDosageSensitivity.gene_symbol == gene_upper
        )
    )
    dosage = result.scalar_one_or_none()
    if not dosage:
        raise HTTPException(
            status_code=404, detail=f"No dosage data for {gene_upper}"
        )

    return GeneDosageSensitivitySchema(
        gene_symbol=dosage.gene_symbol,
        pli_score=dosage.pli_score,
        loeuf_score=dosage.loeuf_score,
        mis_z_score=dosage.mis_z_score,
        risk_class=dosage.risk_class,
    )
