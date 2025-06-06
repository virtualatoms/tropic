"""Settings for the roppy website."""

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """
    Settings for the roppy website.

    The default way to modify these is to modify ~/.roppy.env or through
    through environment variables by using the "ROPPY" prefix.
    E.g. ROPPY_API_ENDPOINT="http://localhost:8000"
    """

    # database settings
    API_ENDPOINT: str = Field(
        "http://localhost:8000",
        description="URL for the roppy API. This is used to fetch data for the website.",
    )

    SEARCH_NUM_ROWS: int = Field(
        10, description="Number of rows to return in the search table."
    )

    model_config = SettingsConfigDict(env_prefix="roppy_", env_file="~/.roppy.env")
