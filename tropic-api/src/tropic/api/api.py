"""Central starting point for the Tropic API."""

from typing import Any

import uvicorn
from beanie import init_beanie
from fastapi import FastAPI, HTTPException
from fastapi_filter import FilterDepends
from fastapi_pagination import add_pagination
from fastapi_pagination.ext.beanie import paginate
from motor.motor_asyncio import AsyncIOMotorClient

from tropic.api import SETTINGS
from tropic.api.documents import MonomerSummaryDocument
from tropic.api.filters import MonomerSummaryFilter
from tropic.api.paginate import BigPage

app = FastAPI(title="Polymerization API")
add_pagination(app)


@app.on_event("startup")
async def startup_event() -> None:
    """Initialize the Beanie ODM with the MongoDB client."""
    client = AsyncIOMotorClient(SETTINGS.DATABASE_URL)
    await init_beanie(
        database=client[SETTINGS.DATABASE_NAME],
        document_models=[MonomerSummaryDocument],
    )


@app.get("/monomers", response_model=BigPage[MonomerSummaryDocument])
async def get_monomers(
    monomer_filter: MonomerSummaryFilter = FilterDepends(MonomerSummaryFilter),  # noqa: B008
) -> Any:
    """Retrieve a paginated list of monomers with optional filtering."""
    query = monomer_filter.filter(MonomerSummaryDocument.find({}))
    query = monomer_filter.sort(query)
    query = query.find(fetch_links=False)
    return await paginate(query)


@app.get("/monomers/{monomer_id}")
async def get_monomer(monomer_id: str) -> MonomerSummaryDocument:
    """Retrieve a specific monomer by its ID."""
    document = await MonomerSummaryDocument.find_one({"monomer_id": monomer_id})
    if not document:
        raise HTTPException(status_code=404, detail="Monomer not found")
    return document


def main() -> None:
    """Run the FastAPI application."""
    uvicorn.run(app, host=SETTINGS.API_HOST, port=SETTINGS.API_PORT)
