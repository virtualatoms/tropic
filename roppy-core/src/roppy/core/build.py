import os
import asyncio
from beanie import init_beanie
from csv import DictReader
from motor.motor_asyncio import AsyncIOMotorClient
from typing import Any, Type
from roppy.core.models import (
    Molecule,
    Monomer,
    Initiator,
    Polymerisation,
    MonomerSummary,
    InitiatorSummary,
    DataRow,
)

CLIENT_URL = "mongodb://localhost:27017"

# WARN: must ensure bijective mapping between smiles strings and chemical structures
# NOTE: achieved using RdKit Smiles standardisation?


def format_polymerisation_data(data: dict[str, str]) -> dict[str, Any]:
    return {
        "type": data["polymerisation_type"],
        "monomer": {"smiles": data["monomer_smiles"]},
        "initiator": {"smiles": data["initiator_smiles"]},
        "product": {
            "smiles": data["polymer_smiles"],
            "length": data["polymer_length"],
            "dispersity": data["dispersity"],
            "deg_of_poly": data["degree_of_polymerisation"],
            "n_avg_molar_mass": data["number_average_molar_mass"],
            "m_avg_molar_mass": data["mass_average_molar_mass"],
        },
        "parameters": {
            "is_experimental": data["is_experimental"],
            "temperature": data["temperature"],
            "pressure": data["pressure"],
            "solvent": data["solvent"],
            "solvent_conc": data["solvent_conc"],
            "monomer_conc": data["monomer_conc"],
            "monomer_state": data["monomer_state"],
            "polymer_state": data["polymer_state"],
            "solvent_model": data["solvent_model"],
            "method": data["method"],
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
            "date": data["date"],
            "comment": data["comment"],
            "doi": data["doi"],
            "url": data["url"],
        },
    }


async def clear_database():
    await Polymerisation.find_all().delete()
    await Monomer.find_all().delete()
    await MonomerSummary.find_all().delete()
    await Initiator.find_all().delete()
    await InitiatorSummary.find_all().delete()


async def parse_data() -> None:

    polys = []
    for input_file in [
        file
        for file in os.listdir("./data")
        if (file.startswith("input") and file.endswith("csv"))
    ]:
        with open(f"./data/{input_file}") as fstream:
            for data in DictReader(fstream):
                polys.append(
                    Polymerisation.model_validate_strings(
                        format_polymerisation_data(data)
                    )
                )

    for i, poly in enumerate(polys):
        poly.polymerisation_id = i

    await Polymerisation.insert_many(polys)


# THIS SHOULD WORK
async def create_monomer_summaries():

    # insert unique monomers into db as collection
    for poly in await Polymerisation.find_all().to_list():
        monomer = await Monomer.find(
            Monomer.smiles == poly.monomer.smiles
        ).first_or_none()
        if not monomer:
            await Monomer.insert_one(poly.monomer)
        else:
            # reassign polymerisation monomer to existing monome
            poly.monomer = monomer
            await poly.save()

    # generate monomer summaries for each unique monomer
    for i, monomer in enumerate(await Monomer.find_all().to_list()):
        # polymerisation(s) guaranteed to exist due to previous step
        polymerisations = await Polymerisation.find(
            Polymerisation.monomer.smiles == monomer.smiles
        ).to_list()
        data = [
            DataRow.model_validate(
                {
                    "type": poly.type,
                    "is_experimental": poly.parameters.is_experimental,
                    "monomer_state": poly.parameters.monomer_state,
                    "polymer_state": poly.parameters.polymer_state,
                    "monomer_conc": poly.parameters.monomer_conc,
                    "solvent": poly.parameters.solvent,
                    "solvent_conc": poly.parameters.solvent_conc,
                    "delta_h": poly.thermo.delta_h,
                    "delta_s": poly.thermo.delta_s,
                    "ceiling_temperature": poly.thermo.ceiling_temperature,
                    "date": poly.metadata.date,
                }
            )
            for poly in polymerisations
        ]

        monomer_summary = MonomerSummary.model_validate(
            {
                "monomer_id": i,
                "monomer": monomer,
                "polymerisations": polymerisations,
                "data": data,
            }
        )

        await MonomerSummary.insert_one(monomer_summary)


async def create_initiator_summaries():

    # insert initiators into db as collection
    for poly in await Polymerisation.find_all().to_list():
        if poly.initiator is not None:
            initiator = await Initiator.find(
                Initiator.smiles == poly.initiator.smiles
            ).first_or_none()
            if not initiator:
                await Initiator.insert_one(poly.initiator)
            else:
                poly.initiator = initiator
                await poly.save()

    # insert initiator summaries into db as collection
    for i, initiator in enumerate(await Initiator.find_all().to_list()):
        polymerisations = (
            await Polymerisation.find(Polymerisation.initiator is not None)
            .find(
                Polymerisation.initiator.smiles == initiator.smiles,
            )
            .to_list()
        )
        initiator_summary = InitiatorSummary.model_validate(
            {
                "initiator_id": i,
                "initiator": initiator,
                "polymerisations": polymerisations,
            }
        )

        await InitiatorSummary.insert_one(initiator_summary)


# UNIFIED FUNCTION TO INSERT MONOMERS AND INITIATORS
async def create_summaries(MolType: Type[Molecule]) -> None:

    field_name = MolType.__name__.lower()

    # insert all molecules into db as collection
    for poly in await Polymerisation.find_all().to_list():
        field = getattr(poly, field_name)
        if field.smiles:
            molecule = await MolType.find(
                MolType.smiles == field.smiles,
            ).first_or_none()
            if not molecule:
                await MolType.insert_one(field)
            else:
                setattr(poly, field_name, field)
                await poly.save()

    # insert molecule summaries into db as collection
    for i, molecule in enumerate(await MolType.find_all().to_list()):
        polymerisations = await Polymerisation.find(
            Polymerisation.monomer.smiles == molecule.smiles
        ).to_list()
        summary = MonomerSummary.model_validate(
            {
                "monomer_id": i,
                "monomer": monomer,
                "polymerisations": polymerisations,
            }
        )

        await MonomerSummary.insert_one(monomer_summary)


async def find_polymerisations():

    # result = await Polymerisation.find({}).first_or_none()
    result = await Polymerisation.find_all().to_list()
    # print(result)
    print(len(result))


async def find_monomers():

    # result = await Monomer.find({}).first_or_none()
    result = await Monomer.find_all().to_list()
    # print(result)
    print(len(result))


async def find_monomer_summaries():

    # result = await Monomer.find({}).first_or_none()
    result = await MonomerSummary.find_all().to_list()
    print(len(result))
    print(result[-1])


async def rebuild_db():

    client = AsyncIOMotorClient(CLIENT_URL)
    database = client["roppy"]

    await init_beanie(
        database=database,
        document_models=[
            Polymerisation,
            Monomer,
            MonomerSummary,
            Initiator,
            InitiatorSummary,
        ],
    )

    await clear_database()
    await parse_data()
    await create_monomer_summaries()
    # await create_initiator_summaries()
    # await find_polymerisations()
    # await find_monomers()
    # await find_monomer_summaries()


if __name__ == "__main__":

    asyncio.run(rebuild_db())
