from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy import text

from app.config import settings
from app.core.database import engine, init_db
from app.api.v1.router import api_router


@asynccontextmanager
async def lifespan(app: FastAPI):
    await init_db()
    yield


app = FastAPI(
    title=settings.app_name,
    description="In Silico Drug Target Safety Platform",
    version="0.1.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(api_router)


@app.get("/")
async def root():
    """Landing route so hitting the base URL returns useful info instead of a
    bare 404 ``{"detail":"Not Found"}``."""
    return {
        "app": settings.app_name,
        "status": "ok",
        "message": f"{settings.app_name} API is running.",
        "docs": "/docs",
        "health": "/health",
        "api_base": "/api/v1",
    }


@app.get("/health")
async def health_check():
    # Liveness stays 200 as long as the process is up (so a DB blip doesn't
    # fail Railway's healthcheck), but we report DB reachability so operators
    # can tell at a glance whether Postgres is actually connected. Only the
    # exception type is surfaced — never the message — to avoid leaking the
    # connection string / credentials.
    database = "ok"
    try:
        async with engine.connect() as conn:
            await conn.execute(text("SELECT 1"))
    except Exception as exc:  # noqa: BLE001 - reported for ops debugging
        database = f"error: {type(exc).__name__}"
    return {"status": "ok", "app": settings.app_name, "database": database}
