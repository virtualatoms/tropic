import asyncio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from tropic.api.documents import PolymerisationDocument, MonomerSummaryDocument

CLIENT_URL = "mongodb://localhost:27017"


def get_experimental(df: pd.Series) -> str:
    if df.any():
        if df.all():
            return "experimental"
        else:
            return "both"
    else:
        return "computational"


def get_func_group(func_groups: list[str]) -> str:
    if (
        ("O=C([O][R])[N]([R])[R]" in func_groups)
        or ("O=[C]([R])[N]([R])[R]" in func_groups)
        or ("O=[C]1[R][O]C[N]1[R]" in func_groups)
        or ("O=C1OC[N]1[R]" in func_groups)
    ):
        return "Lm"
    if (
        ("S=C([O][R])[S][R]" in func_groups)
        or ("S=C([O][R])[O][R]" in func_groups)
        or ("S=C[O][R]" in func_groups)
    ):
        return "tnL"
    elif (
        ("S=C[S][R]" in func_groups)
        or ("O=C([S][R])[S][R]" in func_groups)
        or ("S=C([S][R])[S][R]" in func_groups)
    ):
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
    polys = await PolymerisationDocument.find_all().to_list()
    data = {
        "poly_id": [poly.polymerisation_id for poly in polys],
        "monomer_smiles": [poly.monomer.smiles for poly in polys],
        "func_group": [poly.monomer.functional_groups for poly in polys],
        "is_experimental": [poly.parameters.is_experimental for poly in polys],
    }
    df = pd.DataFrame(data)

    plt.style.use("seaborn-v0_8-paper")
    _, ax = plt.subplots(subplot_kw={"aspect": "equal"})

    # Group dataframe by monomer_smiles and add experimental/func_group columns
    monomer_grouping = df.groupby("monomer_smiles", group_keys=True)
    df = monomer_grouping.first()
    df["experimental"] = monomer_grouping["is_experimental"].apply(get_experimental)
    df["func_group_pretty"] = df["func_group"].apply(get_func_group)

    experimental_count = df.groupby(["experimental"]).count()
    func_group_count = (
        df.groupby(["experimental", "func_group_pretty"]).count().reset_index()
    )
    func_group_count = func_group_count.groupby("experimental").apply(
        lambda x: [x["func_group_pretty"].to_list(), x["poly_id"].to_list()]
    )
    func_group_labels = [
        [f"{fg}: {count}" for fg, count in zip(*pair)]
        for pair in func_group_count.values
    ]
    func_group_labels = ["\n".join(fg_list) for fg_list in func_group_labels]

    wedges, *_ = ax.pie(
        experimental_count["poly_id"],
        # labels=experimental_count.index.to_list(),
        # wedgeprops=dict(width=0.5, edgecolor="w"),
        autopct="%1.0f%%",
        # startangle=90,
    )

    kw = {
        "arrowprops": {"arrowstyle": "-"},
        "bbox": {"boxstyle": "square,pad=0.3", "fc": "w", "ec": "k", "lw": 0.72},
        "zorder": 0,
        "va": "center",
    }
    for i, wedge in enumerate(wedges):
        ang = (wedge.theta2 - wedge.theta1) / 2 + wedge.theta1
        x, y = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(
            func_group_labels[i],
            xy=(x, y),
            xytext=(1.35 * np.sign(x), 1.4 * y),
            horizontalalignment=horizontalalignment,
            **kw,
        )

    ax.set_title("Distribution of Monomer Functional Groups", y=-0.01)
    plt.legend(wedges, experimental_count.index.to_list(), loc="upper right")

    # plt.show()
    plt.savefig("assets/functional_groups.svg")


async def draw():
    client = AsyncIOMotorClient(CLIENT_URL)
    database = client["tropic"]

    await init_beanie(
        database=database,
        document_models=[
            PolymerisationDocument,
            MonomerSummaryDocument,
        ],
    )

    await draw_dataset()


if __name__ == "__main__":
    asyncio.run(draw())
