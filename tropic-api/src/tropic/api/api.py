"""Central starting point for the Tropic API."""

import argparse

import uvicorn
from beanie import init_beanie
from fastapi import FastAPI, HTTPException
from fastapi_filter import FilterDepends
from fastapi_pagination import add_pagination
from fastapi_pagination.ext.beanie import paginate
from motor.motor_asyncio import AsyncIOMotorClient

from tropic.api import SETTINGS
from tropic.api.documents import MonomerDocument, ReactionDocument
from tropic.api.filters import MonomerSummariesFilter, ReactionFilter
from tropic.api.paginate import BigPage

app = FastAPI(title="TROPIC API")
add_pagination(app)


@app.on_event("startup")
async def startup_event() -> None:
    """Initialize the Beanie ODM with the MongoDB client."""
    client = AsyncIOMotorClient(SETTINGS.DATABASE_URL)
    await init_beanie(
        database=client[SETTINGS.DATABASE_NAME],
        document_models=[ReactionDocument, MonomerDocument],
    )


@app.get("/reactions")
async def get_reactions(
    reaction_filter: ReactionFilter = FilterDepends(ReactionFilter),  # noqa: B008
) -> list[ReactionDocument]:
    """Retrieve a paginated list of reactions with optional filtering."""
    query = reaction_filter.filter(ReactionDocument.find({}))
    query = query.find(fetch_links=True)
    return await query.to_list()


@app.get("/reactions/{reaction_id}")
async def get_reaction(reaction_id: str) -> ReactionDocument:
    """Retrieve a specific reaction by its ID."""
    document = await ReactionDocument.find_one(
        {"reaction_id": reaction_id},
        fetch_links=True,
    )
    if not document:
        raise HTTPException(status_code=404, detail="Reaction not found")
    return document


@app.get("/monomer-summaries", response_model=BigPage, include_in_schema=False)
async def get_monomers(
    reaction_filter: MonomerSummariesFilter = FilterDepends(  # noqa: B008
        MonomerSummariesFilter,
    ),
    has_comp: bool | None = None,
    has_exp: bool | None = None,
) -> list[dict]:
    """Private endpoint to retrieve monomer summaries."""
    pipeline = [
        {
            "$group": {
                "_id": "$monomer.smiles",
                "monomer": {
                    "$first": {
                        "monomer_id": "$monomer.monomer_id",
                        "smiles": "$monomer.smiles",
                        "ring_size": "$monomer.ring_size",
                    },
                },
                "has_exp": {"$max": "$parameters.is_experimental"},
                "has_comp": {"$max": {"$not": "$parameters.is_experimental"}},
            },
        },
        {"$sort": {"monomer.monomer_id": 1}},
        {"$project": {"_id": 0}},
    ]
    match = []
    if has_comp is not None:
        match.append({"has_comp": has_comp})
    if has_exp is not None:
        match.append({"has_exp": has_exp})
    if match:
        pipeline += [{"$match": {"$and": match}}]

    query = reaction_filter.filter(ReactionDocument.find({}))
    return await paginate(query.aggregate(pipeline))


@app.get("/monomer-summaries/{monomer_id}", include_in_schema=False)
async def get_monomer(monomer_id: str) -> dict:
    """Private endpoint to retrieve a specific monomer by its ID."""
    pipeline = [
        {
            "$group": {
                "_id": "$monomer.smiles",
                "has_exp": {"$max": "$parameters.is_experimental"},
                "has_comp": {"$max": {"$not": "$parameters.is_experimental"}},
                "monomer": {
                    "$first": {
                        "monomer_id": "$monomer.monomer_id",
                        "smiles": "$monomer.smiles",
                        "inchi": "$monomer.inchi",
                        "molecular_weight": "$monomer.molecular_weight",
                        "functional_groups": "$monomer.functional_groups",
                        "iupac_name": "$monomer.iupac_name",
                        "pubchem_cid": "$monomer.pubchem_cid",
                        "ring_size": "$monomer.ring_size",
                        "xyz": "$monomer.xyz",
                    },
                },
                "data": {
                    "$push": {
                        "type": "$type",
                        "reaction_id": "$reaction_id",
                        "is_experimental": "$parameters.is_experimental",
                        "state_summary": "$parameters.state_summary",
                        "initial_monomer_conc": "$parameters.initial_monomer_conc",
                        "bulk_monomer_conc": "$parameters.bulk_monomer_conc",
                        "medium": "$parameters.medium",
                        "repeating_units": "$product.repeating_units",
                        "delta_h": "$thermo.delta_h",
                        "delta_s": "$thermo.delta_s",
                        "ceiling_temperature": "$thermo.ceiling_temperature",
                        "method": "$parameters.method",
                        "year": "$metadata.year",
                        "doi": "$metadata.doi",
                        "formatted_reference": "$metadata.formatted_reference",
                    },
                },
            },
        },
    ]
    query = ReactionDocument.find(
        {"monomer.monomer_id": monomer_id},
        fetch_links=True,
    )
    documents = await query.aggregate(pipeline).to_list()
    if not documents:
        raise HTTPException(status_code=404, detail="Monomer not found")
    return documents[0]


def main() -> None:
    """Run the FastAPI application."""
    parser = argparse.ArgumentParser(description="Run the Tropic API.")
    parser.add_argument(
        "--host",
        type=str,
        default=SETTINGS.API_HOST,
        help="Host for the API server.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=SETTINGS.API_PORT,
        help="Port for the API server.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=SETTINGS.API_WORKERS,
        help="Number of workers for the API server.",
    )
    args = parser.parse_args()

    uvicorn.run(
        app,
        host=args.host,
        port=args.port,
        workers=args.workers,
    )
