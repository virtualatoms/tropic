from fastapi import FastAPI, HTTPException, Query
from beanie import Document, init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from typing import Optional, List, Any
import uvicorn
from pydantic import BaseModel

# Import our previous models
from models import (
    Monomer,
    Product,
    Initiator,
    Thermo,
    Conditions,
    Parameters,
    PolymerizationModel,
    SOLVENT_ALIASES,
)


# Create Beanie Document Model
class PolymerizationDocument(Document, PolymerizationModel):
    class Settings:
        name = "polymerizations"


# Initialize FastAPI app
app = FastAPI(title="Polymerization API")

# Database connection settings
DATABASE_URL = "mongodb://localhost:27017"
DATABASE_NAME = "polymer_db"


# Database initialization
@app.on_event("startup")
async def startup_event():
    client = AsyncIOMotorClient(DATABASE_URL)
    await init_beanie(
        database=client[DATABASE_NAME],
        document_models=[PolymerizationDocument, MonomerOverview],
    )


# Helper function to build query filters
def build_query_filter(
    monomer_smiles: Optional[str] = None,
    monomer_mw_min: Optional[float] = None,
    monomer_mw_max: Optional[float] = None,
    monomer_functional_group: Optional[str] = None,
    product_repeating_unit: Optional[str] = None,
    product_mw_min: Optional[float] = None,
    product_mw_max: Optional[float] = None,
    initiator_smiles: Optional[str] = None,
    initiator_functional_group: Optional[str] = None,
    delta_h_min: Optional[float] = None,
    delta_h_max: Optional[float] = None,
    delta_s_min: Optional[float] = None,
    delta_s_max: Optional[float] = None,
    temperature_min: Optional[float] = None,
    temperature_max: Optional[float] = None,
    solvent: Optional[str] = None,
    monomer_state: Optional[str] = None,
    polymer_state: Optional[str] = None,
    pressure_min: Optional[float] = None,
    pressure_max: Optional[float] = None,
    method: Optional[str] = None,
    functional: Optional[str] = None,
    basis_set: Optional[str] = None,
    experimental: Optional[bool] = None,
    doi: Optional[str] = None,
) -> dict:
    query = {}

    # Monomer filters
    if monomer_smiles:
        query["monomer.smiles"] = {"$regex": monomer_smiles, "$options": "i"}
    if monomer_mw_min is not None or monomer_mw_max is not None:
        query["monomer.molecular_weight"] = {}
        if monomer_mw_min is not None:
            query["monomer.molecular_weight"]["$gte"] = monomer_mw_min
        if monomer_mw_max is not None:
            query["monomer.molecular_weight"]["$lte"] = monomer_mw_max
    if monomer_functional_group:
        query["monomer.functional_group"] = {
            "$regex": monomer_functional_group,
            "$options": "i",
        }

    # Product filters
    if product_repeating_unit:
        query["product.repeating_unit"] = {
            "$regex": product_repeating_unit,
            "$options": "i",
        }
    if product_mw_min is not None or product_mw_max is not None:
        query["product.molecular_weight"] = {}
        if product_mw_min is not None:
            query["product.molecular_weight"]["$gte"] = product_mw_min
        if product_mw_max is not None:
            query["product.molecular_weight"]["$lte"] = product_mw_max

    # Initiator filters
    if initiator_smiles:
        query["initiator.smiles"] = {"$regex": initiator_smiles, "$options": "i"}
    if initiator_functional_group:
        query["initiator.functional_group"] = {
            "$regex": initiator_functional_group,
            "$options": "i",
        }

    # Thermo filters
    if delta_h_min is not None or delta_h_max is not None:
        query["thermo.delta_h"] = {}
        if delta_h_min is not None:
            query["thermo.delta_h"]["$gte"] = delta_h_min
        if delta_h_max is not None:
            query["thermo.delta_h"]["$lte"] = delta_h_max
    if delta_s_min is not None or delta_s_max is not None:
        query["thermo.delta_s"] = {}
        if delta_s_min is not None:
            query["thermo.delta_s"]["$gte"] = delta_s_min
        if delta_s_max is not None:
            query["thermo.delta_s"]["$lte"] = delta_s_max

    # Conditions filters
    if temperature_min is not None or temperature_max is not None:
        query["conditions.temperature"] = {}
        if temperature_min is not None:
            query["conditions.temperature"]["$gte"] = temperature_min
        if temperature_max is not None:
            query["conditions.temperature"]["$lte"] = temperature_max
    if solvent:
        canonical_solvent = SOLVENT_ALIASES.get(solvent.lower(), solvent)
        query["conditions.solvent"] = canonical_solvent
    if monomer_state:
        query["conditions.monomer_state"] = monomer_state
    if polymer_state:
        query["conditions.polymer_state"] = polymer_state
    if pressure_min is not None or pressure_max is not None:
        query["conditions.pressure"] = {}
        if pressure_min is not None:
            query["conditions.pressure"]["$gte"] = pressure_min
        if pressure_max is not None:
            query["conditions.pressure"]["$lte"] = pressure_max

    # Parameters filters
    if method:
        query["parameters.method"] = method
    if functional:
        query["parameters.functional"] = {"$regex": functional, "$options": "i"}
    if basis_set:
        query["parameters.basis_set"] = {"$regex": basis_set, "$options": "i"}

    # Other filters
    if experimental is not None:
        query["experimental"] = experimental
    if doi:
        query["doi"] = {"$regex": doi, "$options": "i"}

    return query


# API endpoints
@app.post("/polymerization/", response_model=PolymerizationModel)
async def create_polymerization(polymerization: PolymerizationModel):
    document = PolymerizationDocument(**polymerization.dict())
    await document.insert()
    return document


@app.get("/polymerization/{id}", response_model=PolymerizationModel)
async def get_polymerization(id: str):
    document = await PolymerizationDocument.get(id)
    if not document:
        raise HTTPException(status_code=404, detail="Polymerization not found")
    return document


@app.get("/polymerizations/", response_model=List[PolymerizationModel])
async def search_polymerizations(
    # Monomer filters
    monomer_smiles: Optional[str] = Query(
        None, description="Filter by monomer SMILES pattern"
    ),
    monomer_mw_min: Optional[float] = Query(
        None, description="Minimum monomer molecular weight"
    ),
    monomer_mw_max: Optional[float] = Query(
        None, description="Maximum monomer molecular weight"
    ),
    monomer_functional_group: Optional[str] = Query(
        None, description="Filter by monomer functional group"
    ),
    # Product filters
    product_repeating_unit: Optional[str] = Query(
        None, description="Filter by product repeating unit"
    ),
    product_mw_min: Optional[float] = Query(
        None, description="Minimum product molecular weight"
    ),
    product_mw_max: Optional[float] = Query(
        None, description="Maximum product molecular weight"
    ),
    # Initiator filters
    initiator_smiles: Optional[str] = Query(
        None, description="Filter by initiator SMILES pattern"
    ),
    initiator_functional_group: Optional[str] = Query(
        None, description="Filter by initiator functional group"
    ),
    # Thermo filters
    delta_h_min: Optional[float] = Query(None, description="Minimum delta H value"),
    delta_h_max: Optional[float] = Query(None, description="Maximum delta H value"),
    delta_s_min: Optional[float] = Query(None, description="Minimum delta S value"),
    delta_s_max: Optional[float] = Query(None, description="Maximum delta S value"),
    # Conditions filters
    temperature_min: Optional[float] = Query(None, description="Minimum temperature"),
    temperature_max: Optional[float] = Query(None, description="Maximum temperature"),
    solvent: Optional[str] = Query(None, description="Filter by solvent"),
    monomer_state: Optional[str] = Query(None, description="Filter by monomer state"),
    polymer_state: Optional[str] = Query(None, description="Filter by polymer state"),
    pressure_min: Optional[float] = Query(None, description="Minimum pressure"),
    pressure_max: Optional[float] = Query(None, description="Maximum pressure"),
    # Parameters filters
    method: Optional[str] = Query(None, description="Filter by computational method"),
    functional: Optional[str] = Query(None, description="Filter by DFT functional"),
    basis_set: Optional[str] = Query(None, description="Filter by basis set"),
    # Other filters
    experimental: Optional[bool] = Query(
        None, description="Filter by experimental flag"
    ),
    doi: Optional[str] = Query(None, description="Filter by DOI"),
    # Pagination
    skip: int = Query(0, description="Number of records to skip"),
    limit: int = Query(100, description="Maximum number of records to return"),
):
    query = build_query_filter(
        monomer_smiles=monomer_smiles,
        monomer_mw_min=monomer_mw_min,
        monomer_mw_max=monomer_mw_max,
        monomer_functional_group=monomer_functional_group,
        product_repeating_unit=product_repeating_unit,
        product_mw_min=product_mw_min,
        product_mw_max=product_mw_max,
        initiator_smiles=initiator_smiles,
        initiator_functional_group=initiator_functional_group,
        delta_h_min=delta_h_min,
        delta_h_max=delta_h_max,
        delta_s_min=delta_s_min,
        delta_s_max=delta_s_max,
        temperature_min=temperature_min,
        temperature_max=temperature_max,
        solvent=solvent,
        monomer_state=monomer_state,
        polymer_state=polymer_state,
        pressure_min=pressure_min,
        pressure_max=pressure_max,
        method=method,
        functional=functional,
        basis_set=basis_set,
        experimental=experimental,
        doi=doi,
    )

    documents = (
        await PolymerizationDocument.find(query).skip(skip).limit(limit).to_list()
    )
    return documents


@app.put("/polymerization/{id}", response_model=PolymerizationModel)
async def update_polymerization(id: str, polymerization: PolymerizationModel):
    document = await PolymerizationDocument.get(id)
    if not document:
        raise HTTPException(status_code=404, detail="Polymerization not found")

    await document.update({"$set": polymerization.dict(exclude_unset=True)})
    return document


@app.delete("/polymerization/{id}")
async def delete_polymerization(id: str):
    document = await PolymerizationDocument.get(id)
    if not document:
        raise HTTPException(status_code=404, detail="Polymerization not found")

    await document.delete()
    return {"message": "Polymerization deleted"}


app.post("/monomer-overview/", response_model=MonomerOverview)


async def create_monomer_overview(monomer: MonomerOverview):
    document = MonomerOverview(**monomer.dict())
    await document.insert()
    return document


@app.get("/monomer-overview/{monomer_id}", response_model=MonomerOverview)
async def get_monomer_overview(monomer_id: str):
    document = await MonomerOverview.find_one({"monomer_id": monomer_id})
    if not document:
        raise HTTPException(status_code=404, detail="Monomer overview not found")
    return document


@app.get("/monomer-overviews/", response_model=List[MonomerOverview])
async def search_monomer_overviews(
    monomer_id: Optional[str] = Query(None, description="Filter by monomer ID"),
    smiles: Optional[str] = Query(None, description="Filter by SMILES pattern"),
    ring_size: Optional[Union[int, Literal["many"]]] = Query(
        None, description="Filter by ring size"
    ),
    min_delta_h: Optional[float] = Query(None, description="Minimum average delta H"),
    max_delta_h: Optional[float] = Query(None, description="Maximum average delta H"),
    min_delta_s: Optional[float] = Query(None, description="Minimum average delta S"),
    max_delta_s: Optional[float] = Query(None, description="Maximum average delta S"),
    skip: int = Query(0, description="Number of records to skip"),
    limit: int = Query(100, description="Maximum number of records to return"),
):
    query = {}

    if monomer_id:
        query["monomer_id"] = monomer_id
    if smiles:
        query["smiles"] = {"$regex": smiles, "$options": "i"}
    if ring_size is not None:
        query["ring_size"] = ring_size

    # Handle range queries for thermodynamic properties
    if min_delta_h is not None or max_delta_h is not None:
        query["average_delta_h"] = {}
        if min_delta_h is not None:
            query["average_delta_h"]["$gte"] = min_delta_h
        if max_delta_h is not None:
            query["average_delta_h"]["$lte"] = max_delta_h

    if min_delta_s is not None or max_delta_s is not None:
        query["average_delta_s"] = {}
        if min_delta_s is not None:
            query["average_delta_s"]["$gte"] = min_delta_s
        if max_delta_s is not None:
            query["average_delta_s"]["$lte"] = max_delta_s

    documents = await MonomerOverview.find(query).skip(skip).limit(limit).to_list()
    return documents


@app.put("/monomer-overview/{monomer_id}", response_model=MonomerOverview)
async def update_monomer_overview(monomer_id: str, monomer: MonomerOverview):
    document = await MonomerOverview.find_one({"monomer_id": monomer_id})
    if not document:
        raise HTTPException(status_code=404, detail="Monomer overview not found")

    await document.update({"$set": monomer.dict(exclude_unset=True)})
    return document


@app.delete("/monomer-overview/{monomer_id}")
async def delete_monomer_overview(monomer_id: str):
    document = await MonomerOverview.find_one({"monomer_id": monomer_id})
    if not document:
        raise HTTPException(status_code=404, detail="Monomer overview not found")

    await document.delete()
    return {"message": "Monomer overview deleted"}


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
