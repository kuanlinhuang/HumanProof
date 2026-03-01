"""Binding prediction service with pluggable predictor backends.

Architecture:
    BindingPredictor (abstract) <-- MockBindingPredictor (demo)
                                <-- Boltz2Predictor (future, GPU required)

The mock predictor generates realistic-looking binding predictions based on
sequence features (length, composition) matched against the gene database.
When Boltz-2 is available, swap by setting HUMANPROOF_PREDICTOR=boltz2.
"""

import asyncio
import hashlib
import math
import random
import uuid
from abc import ABC, abstractmethod
from datetime import datetime

from sqlalchemy import select, func
from sqlalchemy.ext.asyncio import AsyncSession

from app.models.expression import ExpressionSummary
from app.models.plof import GeneDosageSensitivity
from app.models.job import PredictionJob, BindingPrediction
from app.schemas.sequence import SequenceSubmission
from app.schemas.job import (
    BindingPredictionSchema,
    PredictionJobResult,
    BindingProfileSummary,
    PipelineResult,
)


# ──────────────────────────── ABSTRACT PREDICTOR ──────────────────────────── #

class BindingPredictor(ABC):
    """Abstract interface for binding prediction backends."""

    @abstractmethod
    async def predict(
        self,
        sequence: str,
        sequence_type: str,
        gene_targets: list[dict],
    ) -> list[dict]:
        """Run binding prediction. Returns list of binding result dicts."""
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """Name of the predictor backend."""
        ...


# ──────────────────────────── MOCK PREDICTOR ──────────────────────────────── #

# Realistic target domains for known druggable genes
GENE_BINDING_DOMAINS = {
    "EGFR": {"domain": "ECD Domain III", "residues": "310-514", "type": "orthosteric"},
    "ERBB2": {"domain": "ECD Domain IV", "residues": "510-630", "type": "orthosteric"},
    "PDCD1": {"domain": "IgV-like N-terminal", "residues": "35-145", "type": "orthosteric"},
    "CD274": {"domain": "IgV domain", "residues": "18-134", "type": "orthosteric"},
    "TNF": {"domain": "Trimer interface", "residues": "80-170", "type": "competitive"},
    "VEGFA": {"domain": "Receptor-binding site", "residues": "1-110", "type": "competitive"},
    "IL6": {"domain": "Site I/III", "residues": "42-190", "type": "allosteric"},
    "CTLA4": {"domain": "MYPPPY loop", "residues": "97-107", "type": "orthosteric"},
    "IL1B": {"domain": "Beta-trefoil", "residues": "1-152", "type": "competitive"},
    "TGFB1": {"domain": "Growth factor domain", "residues": "279-390", "type": "competitive"},
    "KIT": {"domain": "ECD D4-D5", "residues": "309-519", "type": "orthosteric"},
    "MET": {"domain": "SEMA domain", "residues": "25-515", "type": "orthosteric"},
    "FGFR1": {"domain": "D2-D3 linker", "residues": "147-365", "type": "allosteric"},
    "ALK": {"domain": "ECD", "residues": "1-1038", "type": "orthosteric"},
}


class MockBindingPredictor(BindingPredictor):
    """
    Demo predictor generating realistic binding profiles.

    Uses sequence hash as a deterministic seed so the same sequence always
    produces the same results. Scoring heuristics:
      - Sequence length and composition influence base affinity
      - Known druggable targets get higher scores
      - Adds realistic noise for delta_G and Kd calculations
    """

    @property
    def name(self) -> str:
        return "mock"

    async def predict(
        self,
        sequence: str,
        sequence_type: str,
        gene_targets: list[dict],
    ) -> list[dict]:
        # Simulate computation time (1-3 seconds)
        await asyncio.sleep(random.uniform(1.0, 3.0))

        # Use sequence hash for deterministic results
        seq_hash = int(hashlib.md5(sequence.encode()).hexdigest(), 16)
        rng = random.Random(seq_hash)

        # Sequence features that influence binding
        seq_len = len(sequence.replace(":", ""))
        hydrophobic_frac = sum(1 for aa in sequence if aa in "AILMFWVP") / max(1, seq_len)
        charged_frac = sum(1 for aa in sequence if aa in "DEKRH") / max(1, seq_len)

        # Base affinity modifier from sequence properties
        if sequence_type == "antibody":
            base_affinity = 0.55 + hydrophobic_frac * 0.3  # antibodies bind well
        elif sequence_type == "nanobody":
            base_affinity = 0.45 + hydrophobic_frac * 0.35  # smaller but can be potent
        else:
            base_affinity = 0.3 + charged_frac * 0.3  # peptides vary more

        results = []
        for target in gene_targets:
            gene = target["gene_symbol"]
            ensembl_id = target["ensembl_id"]

            # Known druggable targets get a boost
            domain_info = GENE_BINDING_DOMAINS.get(gene, None)
            if domain_info:
                target_boost = rng.uniform(0.15, 0.35)
            else:
                target_boost = 0.0

            # Gene-specific random variation (seeded)
            gene_hash = int(hashlib.md5(f"{sequence}{gene}".encode()).hexdigest(), 16)
            gene_rng = random.Random(gene_hash)

            raw_score = base_affinity + target_boost + gene_rng.gauss(0, 0.15)
            binding_score = max(0.01, min(0.99, raw_score))

            # Confidence correlates with binding score (strong predictions are more confident)
            confidence = max(0.3, min(0.99, binding_score + gene_rng.gauss(0.1, 0.1)))

            # Thermodynamic properties from binding score
            # ΔG = RT ln(Kd), approximate
            delta_g = -8.0 * binding_score + gene_rng.gauss(0, 0.5)  # kcal/mol
            kd_nm = 10 ** (3 - 4 * binding_score + gene_rng.gauss(0, 0.3))  # nM scale

            # Binding site info
            if domain_info:
                binding_site = f"{domain_info['domain']} (residues {domain_info['residues']})"
                interaction_type = domain_info["type"]
            else:
                binding_site = f"Predicted surface patch"
                interaction_type = gene_rng.choice(["orthosteric", "allosteric", "competitive"])

            results.append({
                "gene_symbol": gene,
                "ensembl_id": ensembl_id,
                "binding_score": round(binding_score, 4),
                "confidence": round(confidence, 4),
                "binding_site": binding_site,
                "interaction_type": interaction_type,
                "delta_g": round(delta_g, 3),
                "kd_nm": round(max(0.01, kd_nm), 3),
            })

        # Sort by binding score descending and assign ranks
        results.sort(key=lambda x: x["binding_score"], reverse=True)
        for i, r in enumerate(results):
            r["rank"] = i + 1

        return results


# ──────────────────────────── SERVICE LAYER ───────────────────────────────── #

def get_predictor() -> BindingPredictor:
    """Factory: return the configured predictor backend."""
    # Future: check settings.predictor_backend == "boltz2" → Boltz2Predictor
    return MockBindingPredictor()


async def submit_prediction_job(
    db: AsyncSession,
    submission: SequenceSubmission,
) -> PredictionJob:
    """Create a new prediction job record."""
    job = PredictionJob(
        id=str(uuid.uuid4()),
        status="pending",
        sequence_type=submission.sequence_type,
        sequence_name=submission.sequence_name,
        sequence=submission.get_primary_sequence(),
        heavy_chain=submission.heavy_chain,
        light_chain=submission.light_chain,
        predictor_used=get_predictor().name,
    )
    db.add(job)
    await db.commit()
    await db.refresh(job)
    return job


async def run_prediction(db: AsyncSession, job_id: str) -> None:
    """Execute binding prediction for a job (runs as background task)."""
    # Load job
    result = await db.execute(select(PredictionJob).where(PredictionJob.id == job_id))
    job = result.scalar_one_or_none()
    if not job:
        return

    try:
        # Mark as running
        job.status = "running"
        await db.commit()

        # Get all gene targets from the database
        gene_query = await db.execute(
            select(
                GeneDosageSensitivity.gene_symbol,
                GeneDosageSensitivity.ensembl_id,
            ).distinct()
        )
        gene_targets = [
            {"gene_symbol": row.gene_symbol, "ensembl_id": row.ensembl_id}
            for row in gene_query.all()
        ]

        if not gene_targets:
            job.status = "failed"
            job.error_message = "No gene targets found in database."
            await db.commit()
            return

        # Run prediction
        predictor = get_predictor()
        predictions = await predictor.predict(
            sequence=job.sequence,
            sequence_type=job.sequence_type,
            gene_targets=gene_targets,
        )

        # Filter to meaningful binders (score > 0.3)
        significant = [p for p in predictions if p["binding_score"] > 0.3]

        # Re-rank after filtering
        for i, pred in enumerate(significant):
            pred["rank"] = i + 1

        # Store predictions
        for pred in significant:
            binding = BindingPrediction(
                job_id=job_id,
                gene_symbol=pred["gene_symbol"],
                ensembl_id=pred["ensembl_id"],
                binding_score=pred["binding_score"],
                confidence=pred["confidence"],
                binding_site=pred["binding_site"],
                interaction_type=pred["interaction_type"],
                delta_g=pred["delta_g"],
                kd_nm=pred["kd_nm"],
                rank=pred["rank"],
            )
            db.add(binding)

        job.status = "completed"
        job.completed_at = datetime.utcnow()
        job.n_targets_found = len(significant)
        await db.commit()

    except Exception as e:
        job.status = "failed"
        job.error_message = str(e)
        await db.commit()


async def get_pipeline_result(db: AsyncSession, job_id: str) -> PipelineResult | None:
    """
    Get integrated pipeline results: binding predictions enriched with
    safety data (expression + pLOF) for each predicted target.
    """
    # Load job
    result = await db.execute(select(PredictionJob).where(PredictionJob.id == job_id))
    job = result.scalar_one_or_none()
    if not job or job.status != "completed":
        return None

    # Load binding predictions
    pred_result = await db.execute(
        select(BindingPrediction)
        .where(BindingPrediction.job_id == job_id)
        .order_by(BindingPrediction.rank)
    )
    predictions = pred_result.scalars().all()

    # Enrich each prediction with safety data
    profiles = []
    for pred in predictions:
        gene = pred.gene_symbol

        # Get expression breadth
        expr_count = await db.execute(
            select(func.count(func.distinct(ExpressionSummary.tissue)))
            .where(ExpressionSummary.gene_symbol == gene)
        )
        n_tissues = expr_count.scalar() or 0

        # Get top tissue
        top_tissue_q = await db.execute(
            select(ExpressionSummary.tissue, func.max(ExpressionSummary.mean_expression))
            .where(ExpressionSummary.gene_symbol == gene)
            .group_by(ExpressionSummary.tissue)
            .order_by(func.max(ExpressionSummary.mean_expression).desc())
            .limit(1)
        )
        top_tissue_row = top_tissue_q.first()
        top_tissue = top_tissue_row[0] if top_tissue_row else None

        # Get pLOF count
        from app.models.plof import PLOFAssociation
        plof_count = await db.execute(
            select(func.count(PLOFAssociation.id))
            .where(PLOFAssociation.gene_symbol == gene)
        )
        n_plof = plof_count.scalar() or 0

        # Get top phenotype
        top_pheno_q = await db.execute(
            select(PLOFAssociation.phenotype)
            .where(PLOFAssociation.gene_symbol == gene)
            .order_by(PLOFAssociation.p_value.asc())
            .limit(1)
        )
        top_pheno_row = top_pheno_q.first()
        top_phenotype = top_pheno_row[0] if top_pheno_row else None

        # Get risk class
        dosage_q = await db.execute(
            select(GeneDosageSensitivity.risk_class)
            .where(GeneDosageSensitivity.gene_symbol == gene)
        )
        dosage_row = dosage_q.first()
        risk_class = dosage_row[0] if dosage_row else None

        profiles.append(BindingProfileSummary(
            gene_symbol=pred.gene_symbol,
            ensembl_id=pred.ensembl_id,
            binding_score=pred.binding_score,
            confidence=pred.confidence,
            kd_nm=pred.kd_nm,
            interaction_type=pred.interaction_type,
            risk_class=risk_class,
            n_tissues_expressed=n_tissues,
            n_plof_associations=n_plof,
            top_tissue=top_tissue,
            top_phenotype=top_phenotype,
        ))

    return PipelineResult(
        job_id=job.id,
        sequence_name=job.sequence_name,
        sequence_type=job.sequence_type,
        n_targets=len(profiles),
        binding_profiles=profiles,
        completed_at=job.completed_at,
    )
