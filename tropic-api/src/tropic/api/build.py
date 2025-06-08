from __future__ import annotations

import argparse
import asyncio
from collections import defaultdict
from csv import DictReader
from pathlib import Path
from typing import Any

import requests
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from requests_cache import DO_NOT_CACHE, NEVER_EXPIRE, install_cache

from tropic.api import DATABASE_NAME, DATABASE_URL
from tropic.api.documents import MonomerSummaryDocument, PolymerisationDocument
from tropic.core.models import DataRow

install_cache(
    "doi_cache",
    cache_control=True,
    urls_expire_after={
        "*.doi.org": NEVER_EXPIRE,
        "*.nih.gov": NEVER_EXPIRE,
        "*": DO_NOT_CACHE,
    },
)


def format_polymerisation_data(data: dict[str, str]) -> dict[str, Any]:

    pubchem_data = get_iupac_name_cid(data["monomer_smiles"])
    if pubchem_data is None:
        iupac_name, cid = None, None
    else:
        iupac_name, cid = pubchem_data

    return {
        "type": data["polymerisation_type"],
        "monomer": {
            "smiles": data["monomer_smiles"],
            "iupac_name": iupac_name,
            "pubchem_cid": cid,
            # "common_name": get_common_name(data["monomer_smiles"]),
        },
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
            "delta_g": data["delta_g"],
            "ceiling_temperature": data["ceiling_temperature"],
        },
        "metadata": {
            "year": data["date"],
            "comment": data["comment"],
            "doi": data["doi"],
            "url": data["url"],
            "formatted_reference": get_formatted_reference(data["doi"]),
        },
    }


def get_formatted_reference(doi: str) -> str | None:
    """Formats the reference using the citation API."""
    if not doi:
        return None

    response = requests.get(
        f"https://citation.doi.org/format?doi={doi}&style=royal-society-of-chemistry&lang=en-US",
        timeout=5,
    )
    if response.status_code != 200:
        return None

    return response.text.strip()[2:-1]


def get_iupac_name_cid(smiles: str) -> tuple[str, int] | None:
    """Fetches IUPAC name for a given SMILES string."""
    if not smiles:
        return None

    response = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON",
        timeout=5,
    )
    if response.status_code != 200:
        return None

    data = response.json()
    properties = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    return properties.get("IUPACName"), properties.get("CID")


# Common names are a bit wacky from pubchem, so we are not using them for now.
# def get_common_name(smiles: str) -> str:
#     """Fetches common name for a given SMILES string."""
#     if not smiles:
#         return None

#     response = requests.get(
#         f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/synonyms/JSON"
#     )
#     if response.status_code != 200:
#         return None

#     data = response.json()
#     properties = data.get("InformationList", {}).get("Information", [{}])[0]
#     print(properties.get("Synonym", [None]))
#     return properties.get("Synonym", [None])[0]


async def clear_database(summaries_only: bool = False):
    if not summaries_only:
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
            polymerisation_id=poly.polymerisation_id,
            is_experimental=poly.parameters.is_experimental,
            state_summary=poly.parameters.state_summary,
            initial_monomer_conc=poly.parameters.initial_monomer_conc,
            bulk_monomer_conc=poly.parameters.bulk_monomer_conc,
            solvent=poly.parameters.solvent,
            delta_h=poly.thermo.delta_h,
            delta_s=poly.thermo.delta_s,
            ceiling_temperature=poly.thermo.ceiling_temperature,
            year=poly.metadata.year,
            doi=poly.metadata.doi,
            repeating_units=poly.product.repeating_units,
            method=poly.parameters.method,
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


async def rebuild_db(summaries_only: bool = False):
    client = AsyncIOMotorClient(DATABASE_URL)
    await init_beanie(
        database=client[DATABASE_NAME],
        document_models=[MonomerSummaryDocument, PolymerisationDocument],
    )
    await clear_database(summaries_only=summaries_only)

    if not summaries_only:
        await parse_data()

    await create_monomer_summaries()


def main():
    parser = argparse.ArgumentParser(description="Rebuild the TROPIC database.")
    parser.add_argument(
        "-s",
        "--summaries-only",
        action="store_true",
        help="Rebuild the monomer summaries only.",
    )
    args = parser.parse_args()
    asyncio.run(rebuild_db(args.summaries_only))
