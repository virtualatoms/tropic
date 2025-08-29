"""Tool to rebuild the TROPIC database from CSV files."""

import argparse
import asyncio
import logging
import pickle

from beanie import WriteRules, init_beanie
from motor.motor_asyncio import AsyncIOMotorClient

from tropic.api import SETTINGS
from tropic.api.documents import MonomerDocument, ReactionDocument

logging.basicConfig(level=logging.INFO)


async def load_db(filename: str) -> None:
    """Rebuild the TROPIC database."""
    logging.info(f"Connecting to database: {SETTINGS.DATABASE_URL}")  # noqa: G004

    client = AsyncIOMotorClient(
        SETTINGS.DATABASE_URL,
        username=SETTINGS.DATABASE_USERNAME,
        password=SETTINGS.DATABASE_PASSWORD,
        authSource=SETTINGS.DATABASE_AUTH_SOURCE,
    )
    await init_beanie(
        database=client[SETTINGS.DATABASE_NAME],
        document_models=[ReactionDocument, MonomerDocument],
    )

    with open(filename, "rb") as file:  # noqa:PTH123,ASYNC230
        reactions = pickle.load(file)  # noqa:S301

    logging.info("Uploading reactions")
    for reaction in reactions:
        await reaction.save(link_rule=WriteRules.WRITE)
    logging.info("Finished uploading")


def main() -> None:
    """Entry point for the build script."""
    parser = argparse.ArgumentParser(description="Rebuild the TROPIC database.")
    parser.add_argument(
        "-f",
        "--filename",
        default="data/db.pkl",
        type=str,
        help="Path to the pickle file to load.",
    )
    args = parser.parse_args()
    asyncio.run(load_db(args.filename))
