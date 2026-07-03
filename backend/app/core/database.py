import asyncio
import logging

from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.orm import DeclarativeBase

from app.config import settings

logger = logging.getLogger(__name__)

# Bound connection establishment so an unreachable / blackholed / DNS-stalled
# Postgres fails fast instead of hanging startup (``init_db``) and ``/health``.
# asyncpg's ``timeout`` covers DNS resolution + TCP connect; it's asyncpg-only,
# so scope it to that driver (SQLite would misinterpret it as a busy timeout).
CONNECT_TIMEOUT_SECONDS = 10.0
_connect_args = (
    {"timeout": CONNECT_TIMEOUT_SECONDS} if "asyncpg" in settings.database_url else {}
)

engine = create_async_engine(
    settings.database_url,
    echo=False,
    # Recycle connections dropped while the DB was offline (e.g. a Postgres
    # restart) instead of handing out dead ones on the next request.
    pool_pre_ping=True,
    connect_args=_connect_args,
)
async_session = async_sessionmaker(engine, class_=AsyncSession, expire_on_commit=False)


class Base(DeclarativeBase):
    pass


async def get_db() -> AsyncSession:
    async with async_session() as session:
        yield session


async def init_db(retries: int = 5, base_delay: float = 1.0) -> bool:
    """Create tables, retrying to tolerate a database that is still starting up.

    Managed platforms (e.g. Railway) frequently boot the app before Postgres is
    ready, and the database can restart independently. Rather than crash the
    whole service on a transient connection error, retry with exponential
    backoff. Returns True on success. On persistent failure we log and return
    False so the process still starts and ``/health`` can report the database
    as unreachable instead of the container 502-ing with no signal.
    """
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            async with engine.begin() as conn:
                await conn.run_sync(Base.metadata.create_all)
            return True
        except Exception as exc:  # noqa: BLE001 - retried / reported below
            error_text = str(exc).lower()
            if "already exists" in error_text and "sqlite" in settings.database_url:
                return True
            last_exc = exc
            if attempt < retries:
                delay = base_delay * (2 ** (attempt - 1))
                logger.warning(
                    "Database not ready (attempt %d/%d): %s. Retrying in %.1fs.",
                    attempt,
                    retries,
                    type(exc).__name__,
                    delay,
                )
                await asyncio.sleep(delay)

    logger.error(
        "Database initialization failed after %d attempts: %s. "
        "Starting anyway — check /health and the database service.",
        retries,
        type(last_exc).__name__ if last_exc else "unknown",
    )
    return False
