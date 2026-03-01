from pydantic_settings import BaseSettings
from pathlib import Path


class Settings(BaseSettings):
    app_name: str = "HumanProof"
    database_url: str = "sqlite+aiosqlite:///./humanproof.db"
    cors_origins: list[str] = ["http://localhost:3000"]
    data_dir: Path = Path(__file__).parent.parent.parent / "data"

    model_config = {"env_prefix": "HUMANPROOF_"}


settings = Settings()
