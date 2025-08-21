"""Settings for the tropic api."""

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """
    Settings for the tropic website.

    The default way to modify these is to set environment variables with the
    "TROPIC_" prefix. E.g. TROPIC_DATABASE_URL="mongodb://localhost:27017"
    """

    DATABASE_URL: str = Field(
        "mongodb://localhost:27017",
        description="URL for the MongoDB database.",
    )
    DATABASE_NAME: str = Field("tropic", description="Name of the MongoDB database.")
    API_HOST: str = Field("127.0.0.1", description="Host for the API server.")
    API_PORT: int = Field(8000, description="Port for the API server.")
    API_WORKERS: int = Field(1, description="Number of workers for the API server.")
    model_config = SettingsConfigDict(env_prefix="TROPIC_")
