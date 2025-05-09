"""
1. Parse input csvs to generate Polymerisation documents
2. Add Polymerisation objects to the db
"""

import asyncio
from beanie import init_beanie
from csv import DictReader
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.core.models import Polymerisation, Monomer, MonomerSummary

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
async def parse_data() -> None:

    with open("data.csv") as fstream:
        reader = DictReader(fstream)
        polys = [
            Polymerisation.model_validate_strings(format_polymerisation_data(data))
            for data in reader
        ]

    for i, poly in enumerate(polys):
        poly.display_id = i

    await Polymerisation.insert_many(polys)


async def clear_database():
    await Polymerisation.find_all().delete()
    await Monomer.find_all().delete()
    await MonomerSummary.find_all().delete()


async def create_monomer_summaries():

    # insert monomers into db
    for poly in await Polymerisation.find_all().to_list():
        if not await Monomer.find(
            Monomer.smiles == poly.monomer.smiles
        ).first_or_none():
            await Monomer.insert_one(poly.monomer)
            await poly.save()
        else:
            monomer = await Monomer.find_one(Monomer.smiles == poly.monomer.smiles)
            if isinstance(monomer, Monomer):
                poly.monomer = monomer
            await poly.save()

    for i, monomer in enumerate(await Monomer.find_all().to_list()):
        monomer_summary = MonomerSummary.model_validate(
            {
                "display_id": i,
                "monomer": monomer,
                "polymerisations": await Polymerisation.find(
                    Polymerisation.monomer.smiles == monomer.smiles
                ).to_list(),
            }
        )

        await MonomerSummary.insert_one(monomer_summary)


async def find_polymerisations():

    # result = await Polymerisation.find({}).first_or_none()
    result = await Polymerisation.find_all().to_list()
    print(result)
    print(len(result))


async def find_monomers():

    # result = await Monomer.find({}).first_or_none()
    result = await Monomer.find_all().to_list()
    print(result)
    print(len(result))


async def find_monomer_summaries():

    # result = await Monomer.find({}).first_or_none()
    result = await MonomerSummary.find_all().to_list()
    print(result[-1])
    print(len(result))
    for poly in result[-1].polymerisations:
        print(poly.thermo.delta_h)


async def rebuild_db():

    client = AsyncIOMotorClient(CLIENT_URL)
    database = client["roppy"]
    poly_collection = database["polymerisations"]
    monomer_collection = database["monomers"]
    monomer_summary_collection = database["monomer_summaries"]

    await init_beanie(
        database=database, document_models=[Polymerisation, Monomer, MonomerSummary]
    )

    await clear_database()
    await parse_data()
    await create_monomer_summaries()
    # await find_polymerisations()
    # await find_monomers()
    # await find_monomer_summaries()


if __name__ == "__main__":

    asyncio.run(rebuild_db())
