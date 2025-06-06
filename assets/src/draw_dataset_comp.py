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
from matplotlib.patches import ConnectionPatch

CLIENT_URL = "mongodb://localhost:27017"


async def draw_dataset():

    polys = await Polymerisation.find_all().to_list()
    data = {
        "poly_id": [poly.polymerisation_id for poly in polys],
        "is_experimental": [poly.parameters.is_experimental for poly in polys],
        "type": [poly.type for poly in polys],
    }
    df = pd.DataFrame(data)

    plt.style.use("seaborn-v0_8-paper")
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.subplots_adjust(wspace=0)

    experimental_count = df.groupby("is_experimental").count()
    experimental_count["poly_id"] = (
        experimental_count["poly_id"] / experimental_count["poly_id"].sum()
    )
    patches, *_ = ax1.pie(
        experimental_count["poly_id"],
        autopct="%1.1f%%",
        startangle=(180 * experimental_count["poly_id"][1]),
        labels=["comp", "exp"],
        explode=[0.1, 0.0],
    )
    ax1.set_title("Experimental vs Computational")
    type_count = (
        df.groupby(["is_experimental", "type"])
        .count()
        .drop(labels=[True])
        .reset_index(level=0)
    )
    type_count["poly_id"] = type_count["poly_id"] / type_count["poly_id"].sum()
    bar_bottom = 1
    bar_width = 0.2
    for i, (bar_height, label) in enumerate(
        reversed([*zip(type_count["poly_id"], type_count.index)])
    ):
        bar_bottom -= bar_height
        bc = ax2.bar(
            0,
            bar_height,
            bar_width,
            bottom=bar_bottom,
            label=label,
            alpha=(0.1 + 0.25 * i),
        )
        ax2.bar_label(bc, labels=[f"{bar_height:.0%}"], label_type="center")

    ax2.set_title("Type")
    ax2.legend()
    ax2.axis("off")
    ax2.set_xlim(-2.5 * bar_width, 2.5 * bar_width)

    theta1, theta2 = patches[0].theta2, patches[0].theta1
    center, r = patches[0].center, patches[0].r
    bar_height = sum(type_count["poly_id"])

    # draw top connecting line
    x = r * np.cos(np.pi / 180 * theta2) + center[0]
    y = r * np.sin(np.pi / 180 * theta2) + center[1]
    con = ConnectionPatch(
        xyA=(-bar_width / 2, bar_height),
        coordsA=ax2.transData,
        xyB=(x, y),
        coordsB=ax1.transData,
    )
    con.set_color("0")
    con.set_linewidth(2)
    ax2.add_artist(con)

    # draw bottom connecting line
    x = r * np.cos(np.pi / 180 * theta1) + center[0]
    y = r * np.sin(np.pi / 180 * theta1) + center[1]
    con = ConnectionPatch(
        xyA=(-bar_width / 2, 0),
        coordsA=ax2.transData,
        xyB=(x, y),
        coordsB=ax1.transData,
    )
    con.set_color("0")
    ax2.add_artist(con)
    con.set_linewidth(2)

    # plt.show()
    plt.savefig("assets/dataset_comp.svg")


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
