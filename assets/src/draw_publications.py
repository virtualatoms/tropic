import asyncio
import pandas as pd
import matplotlib.pyplot as plt
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.api import DATABASE_NAME, DATABASE_URL
from roppy.api.documents import PolymerisationDocument, MonomerSummaryDocument

# TODO: This script needs to be updated to use the new database schema


async def draw_publications():

    polys = await PolymerisationDocument.find_all().to_list()
    data = {
        "poly_id": [poly.polymerisation_id for poly in polys],
        "date": [poly.metadata.date for poly in polys],
    }

    plt.style.use("seaborn-v0_8-paper")

    df = pd.DataFrame(data)
    df["date"] = pd.to_datetime(df["date"])

    fig, ax1 = plt.subplots()
    ax1.hist(df["date"], bins=120, alpha=0.5)
    ax1.set_xlabel("date")
    ax1.set_ylabel("data count")

    cumulative_count = df.groupby("date").count().cumsum()
    ax2 = ax1.twinx()
    ax2.plot(cumulative_count.index, cumulative_count["poly_id"])
    ax2.set_ylabel("cumulative data count")

    plt.savefig("assets/publications.svg")


async def draw():

    client = AsyncIOMotorClient(DATABASE_URL)
    await init_beanie(
        database=client[DATABASE_NAME],
        document_models=[
            PolymerisationDocument,
            MonomerSummaryDocument,
        ],
    )

    await draw_publications()


if __name__ == "__main__":

    asyncio.run(draw())
