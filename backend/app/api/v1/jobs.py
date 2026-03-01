"""API routes for binding prediction jobs."""

from fastapi import APIRouter, Depends, HTTPException, BackgroundTasks
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.database import get_db, async_session
from app.models.job import PredictionJob, BindingPrediction
from app.schemas.sequence import SequenceSubmission, SequenceValidationResult
from app.schemas.job import (
    PredictionJobStatus,
    PredictionJobResult,
    PredictionJobList,
    BindingPredictionSchema,
    PipelineResult,
)
from app.services.sequence_validator import validate_submission
from app.services.binding_service import (
    submit_prediction_job,
    run_prediction,
    get_pipeline_result,
)

router = APIRouter(prefix="/pipeline", tags=["Pipeline"])


# ─── Validation endpoint ─────────────────────────────────────────────── #

@router.post("/validate", response_model=SequenceValidationResult)
async def validate_sequence(submission: SequenceSubmission):
    """Validate a biologic sequence before submitting a prediction job."""
    return validate_submission(submission)


# ─── Job submission ───────────────────────────────────────────────────── #

@router.post("/jobs", response_model=PredictionJobStatus, status_code=201)
async def create_prediction_job(
    submission: SequenceSubmission,
    background_tasks: BackgroundTasks,
    db: AsyncSession = Depends(get_db),
):
    """
    Submit a biologic sequence for binding prediction.

    The prediction runs asynchronously. Poll GET /pipeline/jobs/{job_id}
    for status updates.
    """
    # Validate first
    validation = validate_submission(submission)
    if not validation.valid:
        raise HTTPException(status_code=422, detail=validation.errors)

    # Create job
    job = await submit_prediction_job(db, submission)

    # Run prediction in background
    async def run_in_background(job_id: str):
        async with async_session() as session:
            await run_prediction(session, job_id)

    background_tasks.add_task(run_in_background, job.id)

    return PredictionJobStatus(
        job_id=job.id,
        status=job.status,
        sequence_type=job.sequence_type,
        sequence_name=job.sequence_name,
        created_at=job.created_at,
        predictor_used=job.predictor_used,
    )


# ─── Job status ───────────────────────────────────────────────────────── #

@router.get("/jobs/{job_id}", response_model=PredictionJobStatus)
async def get_job_status(job_id: str, db: AsyncSession = Depends(get_db)):
    """Get the status of a prediction job."""
    result = await db.execute(
        select(PredictionJob).where(PredictionJob.id == job_id)
    )
    job = result.scalar_one_or_none()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    return PredictionJobStatus(
        job_id=job.id,
        status=job.status,
        sequence_type=job.sequence_type,
        sequence_name=job.sequence_name,
        created_at=job.created_at,
        completed_at=job.completed_at,
        n_targets_found=job.n_targets_found,
        error_message=job.error_message,
        predictor_used=job.predictor_used,
    )


# ─── Job results ──────────────────────────────────────────────────────── #

@router.get("/jobs/{job_id}/binding", response_model=PredictionJobResult)
async def get_binding_results(job_id: str, db: AsyncSession = Depends(get_db)):
    """Get binding prediction results for a completed job."""
    result = await db.execute(
        select(PredictionJob).where(PredictionJob.id == job_id)
    )
    job = result.scalar_one_or_none()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job is not completed. Current status: {job.status}",
        )

    # Load predictions
    pred_result = await db.execute(
        select(BindingPrediction)
        .where(BindingPrediction.job_id == job_id)
        .order_by(BindingPrediction.rank)
    )
    predictions = pred_result.scalars().all()

    return PredictionJobResult(
        job_id=job.id,
        status=job.status,
        sequence_type=job.sequence_type,
        sequence_name=job.sequence_name,
        created_at=job.created_at,
        completed_at=job.completed_at,
        n_targets_found=job.n_targets_found,
        predictor_used=job.predictor_used,
        predictions=[
            BindingPredictionSchema.model_validate(p) for p in predictions
        ],
    )


# ─── Integrated pipeline results ──────────────────────────────────────── #

@router.get("/jobs/{job_id}/pipeline", response_model=PipelineResult)
async def get_pipeline_results(job_id: str, db: AsyncSession = Depends(get_db)):
    """
    Get integrated pipeline results: binding predictions enriched with
    expression and pLOF safety data for each target gene.
    """
    pipeline = await get_pipeline_result(db, job_id)
    if not pipeline:
        raise HTTPException(
            status_code=404,
            detail="Pipeline results not available. Job may not be completed.",
        )
    return pipeline


# ─── List recent jobs ─────────────────────────────────────────────────── #

@router.get("/jobs", response_model=PredictionJobList)
async def list_jobs(
    limit: int = 20,
    offset: int = 0,
    db: AsyncSession = Depends(get_db),
):
    """List recent prediction jobs."""
    # Count total
    from sqlalchemy import func
    count_result = await db.execute(select(func.count(PredictionJob.id)))
    total = count_result.scalar() or 0

    # Fetch jobs
    result = await db.execute(
        select(PredictionJob)
        .order_by(PredictionJob.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    jobs = result.scalars().all()

    return PredictionJobList(
        jobs=[
            PredictionJobStatus(
                job_id=j.id,
                status=j.status,
                sequence_type=j.sequence_type,
                sequence_name=j.sequence_name,
                created_at=j.created_at,
                completed_at=j.completed_at,
                n_targets_found=j.n_targets_found,
                error_message=j.error_message,
                predictor_used=j.predictor_used,
            )
            for j in jobs
        ],
        total=total,
    )
