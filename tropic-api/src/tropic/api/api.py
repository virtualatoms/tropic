"""Central starting point for the Tropic API."""

import logging

import uvicorn
from beanie import init_beanie
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi_filter import FilterDepends
from motor.motor_asyncio import AsyncIOMotorClient

from tropic.api import SETTINGS
from tropic.api.documents import MonomerDocument, ReactionDocument
from tropic.api.filters import MonomerSummariesFilter, ReactionFilter
from tropic.core.models import Reaction

logging.basicConfig(level=logging.INFO)

ORIGINS = [
    "http://localhost",
    "http://localhost:3000",
    f"http://{SETTINGS.API_HOST}:80",
]

logging.info("Allowed origins:")
for origin in ORIGINS:
    logging.info(f" - {origin}")  # noqa: G004
logging.info(f"Using MongoDB database: {SETTINGS.DATABASE_URL}")  # noqa: G004

app = FastAPI(title="TROPIC API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
async def startup_event() -> None:
    """Initialize the Beanie ODM with the MongoDB client."""
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


@app.get("/reactions")
async def get_reactions(
    reaction_filter: ReactionFilter = FilterDepends(ReactionFilter),  # noqa: B008
    include_svg: bool | None = False,
) -> list[Reaction]:
    """Retrieve a list of reactions with optional filtering."""
    query = reaction_filter.filter(ReactionDocument.find({}))
    query = query.find(fetch_links=True)
    documents = await query.to_list()

    if not include_svg:
        for doc in documents:
            doc.monomer.svg = None

    return documents


@app.get("/reactions/{reaction_id}")
async def get_reaction(
    reaction_id: str,
    include_svg: bool | None = False,
) -> Reaction:
    """Retrieve a specific reaction by its ID."""
    document = await ReactionDocument.find_one(
        {"reaction_id": reaction_id},
        fetch_links=True,
    )

    if not include_svg and document:
        document.monomer.svg = None

    if not document:
        raise HTTPException(status_code=404, detail="Reaction not found")
    return document


@app.get("/monomer-summaries", include_in_schema=False)
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
                        "svg": "$monomer.svg",
                    },
                },
                "data": {
                    "$push": {
                        "is_experimental": "$parameters.is_experimental",
                        "delta_h": "$thermo.delta_h",
                        "delta_s": "$thermo.delta_s",
                        "ceiling_temperature": "$thermo.ceiling_temperature",
                        "doi": "$metadata.doi",
                        "formatted_reference": "$metadata.formatted_reference",
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
    return await query.aggregate(pipeline).to_list()


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
                        "svg": "$monomer.svg",
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
    uvicorn.run(
        app,
        host=SETTINGS.API_HOST,
        port=SETTINGS.API_PORT,
        workers=SETTINGS.API_WORKERS,
    )
