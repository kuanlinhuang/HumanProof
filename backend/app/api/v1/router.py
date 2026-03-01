from fastapi import APIRouter

from app.api.v1.expression import router as expression_router
from app.api.v1.plof import router as plof_router
from app.api.v1.targets import router as targets_router
from app.api.v1.jobs import router as jobs_router

api_router = APIRouter(prefix="/api/v1")
api_router.include_router(expression_router)
api_router.include_router(plof_router)
api_router.include_router(targets_router)
api_router.include_router(jobs_router)
