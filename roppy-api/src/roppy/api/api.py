from fastapi import FastAPI, HTTPException
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from typing import Any
import uvicorn
from roppy.api.documents import MonomerSummaryDocument, PolymerisationSummaryDocument
from roppy.api.filters import MonomerSummaryFilter
from roppy.core.models import MonomerSummary, MonomerSummaryBrief
from fastapi_filter import FilterDepends
from fastapi_pagination import Page, add_pagination
from fastapi_pagination.ext.beanie import paginate


app = FastAPI(title="Polymerization API")
add_pagination(app)

DATABASE_URL = "mongodb://localhost:27017"
DATABASE_NAME = "roppy"


@app.on_event("startup")
async def startup_event():
    client = AsyncIOMotorClient(DATABASE_URL)
    await init_beanie(
        database=client[DATABASE_NAME],
        document_models=[MonomerSummaryDocument, PolymerisationSummaryDocument],
    )


@app.get("/monomers", response_model=Page[MonomerSummaryBrief])
async def get_monomers(
    monomer_filter: MonomerSummaryFilter = FilterDepends(MonomerSummaryFilter),
) -> Any:
    query = monomer_filter.filter(MonomerSummaryDocument.find({}))
    query = monomer_filter.sort(query)
    query = query.find(fetch_links=False)
    return await paginate(query.project(MonomerSummaryBrief))


@app.get("/monomers/{monomer_id}", response_model=MonomerSummary)
async def get_monomer(monomer_id: str):
    document = await MonomerSummaryDocument.find_one(
        {"monomer_id": monomer_id}, fetch_links=True
    )
    if not document:
        raise HTTPException(status_code=404, detail="Monomer not found")
    return document


def main():
    uvicorn.run(app, host="0.0.0.0", port=8000)
