import asyncio
from beanie import init_beanie
from csv import DictReader
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.core.models import Polymerisation, Monomer, MonomerSummary

CLIENT_URL = "mongodb://localhost:27017"

# WARN: must ensure bijective mapping between smiles strings and chemical structures
# NOTE: achieved using RdKit Smiles standardisation?


def format_polymerisation_data(data: dict[str, str]) -> dict[str, dict[str, str]]:
    return {
        "monomer": {"smiles": data["monomer_smiles"]},
        "initiator": {"smiles": data["initiator_smiles"]},
        "product": {
            "repeating_unit": data["polymer_repeating_unit"],
            "dispersity": data["polymer_dispersity"],
            "deg_of_poly": data["polymer_degree_of_polymerisation"],
            "n_avg_molar_mass": data["polymer_number_average_molar_mass"],
            "m_avg_molar_mass": data["polymer_mass_average_molar_mass"],
        },
        "parameters": {
            "temperature": data["temperature"],
            "pressure": data["pressure"],
            "solvent": data["solvent"],
            "solvent_conc": data["solvent_conc"],
            "monomer_conc": data["monomer_conc"],
            "monomer_state": data["monomer_state"],
            "polymer_state": data["polymer_state"],
            "solvent_model": data["solvent_model"],
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


async def clear_database():
    await Polymerisation.find_all().delete()
    await Monomer.find_all().delete()
    await MonomerSummary.find_all().delete()


async def parse_data() -> None:

    with open("./data/data.csv") as fstream:
        polys = [
            Polymerisation.model_validate_strings(format_polymerisation_data(data))
            for data in DictReader(fstream)
        ]

    for i, poly in enumerate(polys):
        poly.display_id = i

    await Polymerisation.insert_many(polys)


async def create_monomer_summaries():

    # insert monomers into db as collection
    for poly in await Polymerisation.find_all().to_list():
        monomer = await Monomer.find(
            Monomer.smiles == poly.monomer.smiles
        ).first_or_none()
        if not monomer:
            await Monomer.insert_one(poly.monomer)
        else:
            poly.monomer = monomer
            await poly.save()

    # insert monomer summaries into db as collection
    for i, monomer in enumerate(await Monomer.find_all().to_list()):
        polymerisations = await Polymerisation.find(
            Polymerisation.monomer.smiles == monomer.smiles
        ).to_list()
        monomer_summary = MonomerSummary.model_validate(
            {
                "display_id": i,
                "monomer": monomer,
                # "polymerisations": polymerisations,
                "exp_polymerisations": [
                    poly for poly in polymerisations if poly.is_experimental
                ],
                "comp_polymerisations": [
                    poly for poly in polymerisations if not poly.is_experimental
                ],
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
        database=database, document_models=[Polymerisation, Monomer, MonomerSummary]
    )

    await clear_database()
    await parse_data()
    await create_monomer_summaries()
    await find_polymerisations()
    # await find_monomers()
    await find_monomer_summaries()


if __name__ == "__main__":

    asyncio.run(rebuild_db())
