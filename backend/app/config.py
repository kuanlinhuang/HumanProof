import os
from pathlib import Path

from pydantic import field_validator
from pydantic_settings import BaseSettings


def normalize_database_url(url: str) -> str:
    """Coerce a Postgres URL into the async driver form SQLAlchemy needs.

    Managed hosts (Railway, Heroku, Render, …) inject a database URL using the
    ``postgres://`` or ``postgresql://`` scheme. SQLAlchemy's async engine
    requires an explicit async driver, i.e. ``postgresql+asyncpg://``. Passing
    the raw scheme makes ``create_async_engine`` fail at import time, which
    takes the whole service down. Normalize it here so pointing the app at a
    managed Postgres "just works".
    """
    if url.startswith("postgres://"):
        url = "postgresql://" + url[len("postgres://") :]
    if url.startswith("postgresql://"):
        url = "postgresql+asyncpg://" + url[len("postgresql://") :]
    return url


class Settings(BaseSettings):
    app_name: str = "HumanProof"
    database_url: str = "sqlite+aiosqlite:///./humanproof.db"
    cors_origins: list[str] = ["http://localhost:3000"]
    data_dir: Path = Path(__file__).parent.parent.parent / "data"

    model_config = {"env_prefix": "HUMANPROOF_"}

    @field_validator("database_url", mode="before")
    @classmethod
    def _coerce_database_url(cls, value):
        return normalize_database_url(value) if isinstance(value, str) else value


def _resolve_settings() -> Settings:
    # Railway/Heroku-style platforms inject a bare ``DATABASE_URL`` (not the
    # app-prefixed ``HUMANPROOF_DATABASE_URL``). Honor it when the prefixed
    # variable is not set, so simply attaching a managed Postgres connects the
    # app instead of silently falling back to the local SQLite default.
    if "HUMANPROOF_DATABASE_URL" not in os.environ and os.environ.get("DATABASE_URL"):
        return Settings(database_url=os.environ["DATABASE_URL"])
    return Settings()


settings = _resolve_settings()
