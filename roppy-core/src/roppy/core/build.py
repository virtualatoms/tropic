"""
1. Parse input csvs to generate Polymerisation documents
2. Add Polymerisation objects to the db
"""

import asyncio
from beanie import init_beanie
from csv import DictReader
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.core.models import Polymerisation, MonomerSummary

CLIENT_URL = "mongodb://localhost:27017"

# WARN: must ensure bijective mapping between smiles strings and chemical structures


def format_polymerisation_data(data: dict[str, str]) -> dict[str, dict[str, str]]:
    return {
        "monomer": {"smiles": data["monomer_smiles"]},
        "initiator": {"smiles": data["initiator_smiles"]},
        "product": {
            "number_of_units": data["product_number_of_units"],
            "repeating_unit": data["product_repeating_unit"],
        },
        "parameters": {
            "temperature": data["temperature"],
            "pressure": data["pressure"],
            "solvent": data["solvent"],
            "solvent_conc": data["solvent_conc"],
            "solvent_model": data["solvent_model"],
            "monomer_state": data["monomer_state"],
            "polymer_state": data["polymer_state"],
            "method": data["comp_method"],
            "functional": data["comp_functional"],
            "basis_set": data["comp_basis_set"],
            "dispersion": data["comp_dispersion"],
            "forcefield": data["comp_forcefield"],
        },
        "thermo": {
            "delta_h": data["delta_h"],
            "delta_s": data["delta_s"],
            "ceiling_temperature": data["ceiling_temperature"],
        },
        "metadata": {
            "comment": data["polymerisation_comment"],
            "doi": data["polymerisation_doi"],
            "url": data["polymerisation_link"],
        },
    }


# TODO: create defined directory for data csvs
async def parse_polymerisations():
    with open("data.csv") as fstream:
        reader = DictReader(fstream)
        polys = [
            Polymerisation.model_validate_strings(format_polymerisation_data(data))
            for data in reader
        ]
        await Polymerisation.insert_many(polys)


async def find_polymerisations():

    # async for poly in Polymerisation.find({})
    result = await Polymerisation.find({}).first_or_none()
    print(result)


async def rebuild_db():

    client = AsyncIOMotorClient(CLIENT_URL)
    database = client["roppy"]
    polymerisations = database["polymerisations"]

    await init_beanie(
        database=database, document_models=[Polymerisation, MonomerSummary]
    )

    await parse_polymerisations()
    await find_polymerisations()


if __name__ == "__main__":

    asyncio.run(rebuild_db())
