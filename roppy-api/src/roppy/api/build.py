import asyncio
from collections import defaultdict
from csv import DictReader
from pathlib import Path
from typing import Any

import requests_cache
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient

from roppy.api import DATABASE_NAME, DATABASE_URL
from roppy.api.documents import MonomerSummaryDocument, PolymerisationDocument
from roppy.core.models import DataRow

session = requests_cache.CachedSession("doi_cache")


def format_polymerisation_data(data: dict[str, str]) -> dict[str, Any]:
    formatted_reference = get_formatted_reference(data.get("doi"))
    return {
        "type": data["polymerisation_type"],
        "monomer": {"smiles": data["monomer_smiles"]},
        "initiator": {"smiles": data["initiator_smiles"]},
        "product": {
            "smiles": data["polymer_smiles"],
            "repeating_units": data["repeating_units"],
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
            "initial_monomer_conc": data["initial_monomer_conc"],
            "bulk_monomer_conc": data["bulk_monomer_conc"],
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
            "formatted_reference": formatted_reference,
        },
    }


def get_formatted_reference(doi: str) -> str:
    """Formats the reference using the citation API."""
    if not doi:
        return None

    response = session.get(
        f"https://citation.doi.org/format?doi={doi}&style=royal-society-of-chemistry&lang=en-US"
    )
    if response.status_code != 200:
        return None
    return response.text.strip()[2:-1]


async def clear_database():
    await PolymerisationDocument.find_all().delete()
    await MonomerSummaryDocument.find_all().delete()


async def parse_data() -> None:
    polys = []
    for input_file in Path("data").glob("input*.csv"):
        with open(input_file) as fstream:
            for data in DictReader(fstream):
                polys.append(PolymerisationDocument(**format_polymerisation_data(data)))

    for i, poly in enumerate(polys):
        poly.polymerisation_id = f"poly-{i + 1}"

    await PolymerisationDocument.insert_many(polys)


async def create_monomer_summaries():
    summaries = defaultdict(list)
    monomers = {}
    polymerisations = await PolymerisationDocument.find_all().to_list()
    for poly in polymerisations:
        row = DataRow(
            type=poly.type,
            is_experimental=poly.parameters.is_experimental,
            monomer_state=poly.parameters.monomer_state,
            polymer_state=poly.parameters.polymer_state,
            initial_monomer_conc=poly.parameters.initial_monomer_conc,
            bulk_monomer_conc=poly.parameters.bulk_monomer_conc,
            solvent=poly.parameters.solvent,
            delta_h=poly.thermo.delta_h,
            delta_s=poly.thermo.delta_s,
            ceiling_temperature=poly.thermo.ceiling_temperature,
            date=poly.metadata.date,
            doi=poly.metadata.doi,
            formatted_reference=poly.metadata.formatted_reference,
        )
        summaries[poly.monomer.smiles].append(row)
        monomers[poly.monomer.smiles] = poly.monomer

    monomer_summaries = []
    for i, (monomer_smiles, monomer) in enumerate(monomers.items()):
        monomer_summary = MonomerSummaryDocument(
            monomer_id=f"monomer-{i + 1}",
            monomer=monomer,
            data=summaries[monomer_smiles],
        )
        monomer_summaries.append(monomer_summary)

    await MonomerSummaryDocument.insert_many(monomer_summaries)


async def rebuild_db():
    client = AsyncIOMotorClient(DATABASE_URL)
    await init_beanie(
        database=client[DATABASE_NAME],
        document_models=[MonomerSummaryDocument, PolymerisationDocument],
    )
    await clear_database()
    await parse_data()
    await create_monomer_summaries()


def main():
    asyncio.run(rebuild_db())
