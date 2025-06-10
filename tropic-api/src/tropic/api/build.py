"""Tool to rebuild the TROPIC database from CSV files."""

import argparse
import asyncio
from pathlib import Path
from typing import Any

import requests
from beanie import init_beanie, WriteRules
from motor.motor_asyncio import AsyncIOMotorClient
from requests_cache import DO_NOT_CACHE, NEVER_EXPIRE, install_cache
from tqdm import tqdm

from tropic.api import SETTINGS
from tropic.api.documents import PolymerisationDocument, MonomerDocument

install_cache(
    "doi_cache",
    cache_control=True,
    urls_expire_after={
        "*.doi.org": NEVER_EXPIRE,
        "*.nih.gov": NEVER_EXPIRE,
        "*": DO_NOT_CACHE,
    },
)


def get_polymerisation_document(
    data: dict[str, str], monomers: dict[str, MonomerDocument], polymerisation_id: int
) -> dict[str, Any]:
    """Format the polymerisation data into a structured dictionary."""
    if data["monomer_smiles"] not in monomers:
        # iupac_name, cid = get_iupac_name_cid(data["monomer_smiles"])
        iupac_name = None
        cid = None
        monomers[data["monomer_smiles"]] = MonomerDocument(
            smiles=data["monomer_smiles"],
            iupac_name=iupac_name,
            pubchem_cid=cid,
            monomer_id=f"monomer-{len(monomers) + 1}",
        )
    return PolymerisationDocument(
        **{
            "polymerisation_id": f"poly-{polymerisation_id}",
            "type": data["polymerisation_type"],
            "monomer": monomers[data["monomer_smiles"]],  # use the existing monomer
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
                "initiator": data["initiator_smiles"],
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
    )


def get_formatted_reference(doi: str) -> str | None:
    """Format the reference using the citation API."""
    return None
    if not doi:
        return None

    response = requests.get(
        f"https://citation.doi.org/format?doi={doi}&style=royal-society-of-chemistry&lang=en-US",
        timeout=5,
    )
    if response.status_code != requests.codes.OK:
        return None

    return response.text.strip()[2:-1]


def get_iupac_name_cid(smiles: str) -> tuple[str, int] | tuple[None, None]:
    """Fetch IUPAC name for a given SMILES string."""
    if not smiles:
        return None, None

    response = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON",
        timeout=5,
    )
    if response.status_code != requests.codes.OK:
        return None, None

    data = response.json()
    properties = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    return properties.get("IUPACName"), properties.get("CID")

async def create_polymerisations_monomers(monomers: dict[str, MonomerDocument]) -> None:
    """Parse the input CSV files and create PolymerisationDocument instances."""
    import numpy as np
    import pandas as pd

    all_files = Path("data").glob("input*.csv")
    data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
    data = data.replace({np.nan: None})
    for polymerisation_id, row in tqdm(data.iterrows(), total=data.shape[0]):
        poly = get_polymerisation_document(row, monomers, polymerisation_id + 1)
        await poly.save(link_rule=WriteRules.WRITE)
        


async def build_db(skip_monomers: bool = False) -> None:
    """Rebuild the TROPIC database."""
    client = AsyncIOMotorClient(SETTINGS.DATABASE_URL)
    await init_beanie(
        database=client[SETTINGS.DATABASE_NAME],
        document_models=[PolymerisationDocument, MonomerDocument],
    )

    if not skip_monomers:
        await MonomerDocument.find_all().delete()
    await PolymerisationDocument.find_all().delete()

    monomers = {
        monomer.smiles: monomer for monomer 
        in await MonomerDocument.find_all().to_list()
    }

    await create_polymerisations_monomers(monomers)



def main() -> None:
    """Entry point for the build script."""
    parser = argparse.ArgumentParser(description="Rebuild the TROPIC database.")
    parser.add_argument(
        "-s",
        "--skip-monomers",
        action="store_true",
        help="Skip rebuilding the monomers.",
    )
    args = parser.parse_args()
    asyncio.run(build_db(args.skip_monomers))
