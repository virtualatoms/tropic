"""Tool to rebuild the TROPIC database from CSV files."""

import argparse
import asyncio
import base64
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from beanie import WriteRules, init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles
from requests_cache import DO_NOT_CACHE, NEVER_EXPIRE, install_cache
from scipy.constants import physical_constants
from scipy.stats import linregress
from tqdm import tqdm

from tropic.api import SETTINGS
from tropic.api.documents import MonomerDocument, ReactionDocument

logging.basicConfig(level=logging.INFO)


RDLogger.DisableLog("rdApp.*")

install_cache(
    "data/doi_cache",
    cache_control=True,
    urls_expire_after={
        "*.doi.org": NEVER_EXPIRE,
        "*.nih.gov": NEVER_EXPIRE,
        "*": DO_NOT_CACHE,
    },
)


def get_reaction_document(
    data: dict[str, str],
    monomers: dict[str, MonomerDocument],
    reaction_id: int,
) -> ReactionDocument:
    """Format the reaction data into a structured dictionary."""
    monomer_smiles = StandardizeSmiles(data["monomer_smiles"])
    if monomer_smiles not in monomers:
        iupac_name, cid = get_iupac_name_cid(data["monomer_smiles"])
        monomers[monomer_smiles] = MonomerDocument(
            smiles=monomer_smiles,
            iupac_name=iupac_name,
            pubchem_cid=cid,
            monomer_id=f"monomer-{len(monomers) + 1}",
            svg=smiles_to_image(monomer_smiles),
        )

    vanthoff_data = None
    if data["vanthoff_file"]:
        vanthoff_file = f"data/vanthoff/{data['vanthoff_file']}"
        if not Path(vanthoff_file).exists():
            msg = f"Van't Hoff file not found: {data['vanthoff_file']}"
            logging.warning(msg)
        else:
            vanthoff_data = load_vanthoff(vanthoff_file)

    extrapolation_data = None
    if data["extrapolation_file"]:
        extrap_file = f"data/extrapolation/{data['extrapolation_file']}"
        if not Path(extrap_file).exists():
            msg = f"Extrapolation file not found: {data['extrapolation_file']}"
            logging.warning(msg)
        else:
            extrapolation_data = load_extrapolation(extrap_file)

    return ReactionDocument(
        reaction_id=f"reaction-{reaction_id}",
        type=data["polymerisation_type"],
        monomer=monomers[monomer_smiles],
        product={
            "smiles": data["polymer_smiles"],
            "repeating_units": data["repeating_units"],
            "dispersity": data["dispersity"],
            "deg_of_poly": data["degree_of_polymerisation"],
            "n_avg_molar_mass": data["number_average_molar_mass"],
            "m_avg_molar_mass": data["mass_average_molar_mass"],
        },
        parameters={
            "is_experimental": data["is_experimental"],
            "temperature": data["temperature"],
            "pressure": data["pressure"],
            "solvent": data["solvent"],
            "cosolvent": data["cosolvent"],
            "solvent_cosolvent_ratio": data["solvent/cosolvent"],
            "initiator": data["initiator_smiles"],
            "initial_monomer_conc": data["initial_monomer_conc"],
            "bulk_monomer_conc": data["bulk_monomer_conc"],
            "medium": data["medium"],
            "monomer_state": data["monomer_state"],
            "polymer_state": data["polymer_state"],
            "solvent_model": data["solvent_model"],
            "method": data["method"],
            "functional": data["comp_functional"],
            "basis_set": data["comp_basis_set"],
            "dispersion": data["comp_dispersion"],
            "forcefield": data["comp_forcefield"],
        },
        thermo={
            "delta_h": data["delta_h"],
            "delta_h_std": data["delta_h_std"],
            "delta_s": data["delta_s"],
            "delta_s_std": data["delta_s_std"],
            "delta_g": data["delta_g"],
            "ceiling_temperature": data["reported_ceiling_temperature"],
            "vanthoff": vanthoff_data,
            "extrapolation": extrapolation_data,
        },
        metadata={
            "year": data["date"],
            "comment": data["comment"],
            "doi": data["doi"],
            "url": data["url"],
            "formatted_reference": get_formatted_reference(data["doi"]),
            "flag": data["flag"],
        },
    )


def get_formatted_reference(doi: str) -> str | None:
    """Format the reference using the citation API."""
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


def smiles_to_image(smiles: str, size: tuple[int, int] = (150, 100)) -> str:
    """Convert a SMILES string to a base64-encoded SVG image."""
    drawer = Draw.MolDraw2DSVG(*size)
    dopts = drawer.drawOptions()
    dopts.bondLineWidth = 1.0  # default is 2.
    dopts.clearBackground = False  # default is True
    mol = Chem.MolFromSmiles(smiles)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return base64.b64encode(svg.encode("utf-8")).decode("utf-8")


def load_vanthoff(filename: str) -> dict[str, list[float]]:
    """Load Vanthoff plot data."""
    gas_constant = physical_constants["molar gas constant"][0]

    df_vanthoff = pd.read_excel(filename)
    temperature = df_vanthoff["T (K)"]
    equilibrium_concentration = df_vanthoff["Measured [M]eq"]

    inverse_temperature = 1 / temperature
    r_ln_equilibrium_concentration = gas_constant * np.log(equilibrium_concentration)

    return {
        "temperature": temperature.tolist(),
        "inverse_temperature": inverse_temperature.tolist(),
        "equilibrium_concentration": equilibrium_concentration.tolist(),
        "r_ln_equilibrium_concentration": r_ln_equilibrium_concentration.tolist(),
    }


def load_extrapolation(filename: str) -> dict[str, list[float]]:
    """Load a computational extrapolation file."""
    df_extrap = pd.read_csv(filename)

    x = 1 / df_extrap["repeating_units"]
    y = df_extrap["delta_h"]
    slope, intercept, _, _, std_err = linregress(x, y)
    return {
        "repeating_units": df_extrap["repeating_units"].tolist(),
        "inverse_repeating_units": (1 / df_extrap["repeating_units"]).tolist(),
        "delta_h": df_extrap["delta_h"].tolist(),
        "slope": slope,
        "intercept": intercept,
        "std_err": np.sqrt(std_err),
    }


async def create_reactions_monomers(monomers: dict[str, MonomerDocument]) -> None:
    """Parse the input CSV files and create ReactionDocument instances."""
    all_files = Path("data").glob("input*.xlsx")
    data = pd.concat((pd.read_excel(f) for f in all_files), ignore_index=True)
    data = data.replace({np.nan: None})
    for reaction_id, row in tqdm(data.iterrows(), total=data.shape[0]):
        reaction = get_reaction_document(row, monomers, reaction_id + 1)
        await reaction.save(link_rule=WriteRules.WRITE)


async def build_db(skip_monomers: bool = False) -> None:
    """Rebuild the TROPIC database."""
    logging.info(f"Connecting to database: {SETTINGS.DATABASE_URL}")  # noqa: G004

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

    if not skip_monomers:
        await MonomerDocument.find_all().delete()
    await ReactionDocument.find_all().delete()

    monomers = {
        monomer.smiles: monomer
        for monomer in await MonomerDocument.find_all().to_list()
    }

    await create_reactions_monomers(monomers)


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
