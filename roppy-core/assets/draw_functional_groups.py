import asyncio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.core.models import (
    Monomer,
    Initiator,
    Polymerisation,
    MonomerSummary,
    InitiatorSummary,
)

CLIENT_URL = "mongodb://localhost:27017"


def get_experimental(df: pd.Series) -> str:

    if df.any():
        if df.all():
            return "experimental"
        else:
            return "both"
    else:
        return "computational"


def get_func_group(func_groups: str) -> str:

    if (
        ("O=C([O][R])[N]([R])[R]" in func_groups)
        or ("O=[C]([R])[N]([R])[R]" in func_groups)
        or ("O=[C]1[R][O]C[N]1[R]" in func_groups)
        or ("O=C1OC[N]1[R]L" in func_groups)
    ):
        return "Lm"
    elif ("S=C[S][R]" in func_groups) or ("O=C([S][R])[S][R]" in func_groups):
        return "dtL"
    elif (
        ("O=[C]([R])[S][R]" in func_groups)
        or ("O=[C]([R])SC[O][R]" in func_groups)
        or ("O=C([O][R])[S][R]" in func_groups)
    ):
        return "tL"
    elif ("O=C([O][R])[O][R]" in func_groups) or ("O=C1OCO1" in func_groups):
        return "CC"
    elif (
        ("O=[C]([R])[O][R]" in func_groups)
        or ("C=CC(=O)[O][R]" in func_groups)
        or ("C=CO[C](=O)[R]" in func_groups)
        or ("O=[C](O)[R]" in func_groups)
    ):
        return "L"
    else:
        return "other"


async def draw_dataset():

    polys = await Polymerisation.find_all().to_list()
    data = {
        "poly_id": [poly.polymerisation_id for poly in polys],
        "monomer_smiles": [poly.monomer.smiles for poly in polys],
        "func_group": [poly.monomer.functional_groups for poly in polys],
        "is_experimental": [poly.parameters.is_experimental for poly in polys],
        "type": [poly.type for poly in polys],
    }
    df = pd.DataFrame(data)
    # print(df["func_group"].to_string())

    plt.style.use("seaborn-v0_8")

    monomer_grouping = df.groupby("monomer_smiles", group_keys=True)
    df = monomer_grouping.first()
    df["experimental"] = monomer_grouping["is_experimental"].apply(get_experimental)
    df["func_group_pretty"] = df["func_group"].apply(get_func_group)
    print(df[["func_group", "func_group_pretty"]].to_string())

    # func_group_type = df.groupby("monomer_smiles", group_keys=True).first()
    # print(func_group_type.to_string)

    # print(experimental_type.to_string())

    # plt.show()
    plt.savefig("assets/experimental_dataset.svg")


async def draw():

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

    await draw_dataset()


if __name__ == "__main__":

    asyncio.run(draw())
