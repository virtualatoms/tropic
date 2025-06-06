from typing import Any

import uvicorn
from beanie import init_beanie
from fastapi import FastAPI, HTTPException
from fastapi_filter import FilterDepends
from fastapi_pagination import add_pagination
from fastapi_pagination.ext.beanie import paginate
from motor.motor_asyncio import AsyncIOMotorClient

from roppy.api import DATABASE_NAME, DATABASE_URL
from roppy.api.documents import MonomerSummaryDocument
from roppy.api.filters import MonomerSummaryFilter
from roppy.api.paginate import BigPage
from roppy.core.models import MonomerSummary

app = FastAPI(title="Polymerization API")
add_pagination(app)


@app.on_event("startup")
async def startup_event():
    client = AsyncIOMotorClient(DATABASE_URL)
    await init_beanie(
        database=client[DATABASE_NAME],
        document_models=[MonomerSummaryDocument],
    )


@app.get("/monomers", response_model=BigPage[MonomerSummary])
async def get_monomers(
    monomer_filter: MonomerSummaryFilter = FilterDepends(MonomerSummaryFilter),  # noqa: B008
) -> Any:
    query = monomer_filter.filter(MonomerSummaryDocument.find({}))
    query = monomer_filter.sort(query)
    query = query.find(fetch_links=False)
    return await paginate(query)


@app.get("/monomers/{monomer_id}", response_model=MonomerSummary)
async def get_monomer(monomer_id: str):
    document = await MonomerSummaryDocument.find_one({"monomer_id": monomer_id})
    if not document:
        raise HTTPException(status_code=404, detail="Monomer not found")
    return document


def main():
    uvicorn.run(app, host="0.0.0.0", port=8000)
