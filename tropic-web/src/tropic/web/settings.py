"""Settings for the tropic website."""

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """
    Settings for the tropic website.

    The default way to modify these is to set environment variables with the
    "TROPIC_" prefix. E.g. TROPIC_API_ENDPOINT="http://localhost:8000"
    """

    API_ENDPOINT: str = Field(
        "http://localhost:8000",
        description="URL for the tropic API. This is used to fetch data for the website.",
    )
    REQUEST_TIMEOUT: int = Field(
        5, description="Number of seconds before API requests timeout."
    )
    model_config = SettingsConfigDict(env_prefix="TROPIC_")
