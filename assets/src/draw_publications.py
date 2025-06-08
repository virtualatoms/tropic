import asyncio
import pandas as pd
import matplotlib.pyplot as plt
from beanie import init_beanie
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.api.documents import PolymerisationDocument, MonomerSummaryDocument

CLIENT_URL = "mongodb://localhost:27017"


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
    ax1.hist(df["date"], bins=60, alpha=0.5)
    ax1.set_xlabel("date")
    ax1.set_ylabel("data count")

    cumulative_count = df.groupby("date").count().cumsum()
    ax2 = ax1.twinx()
    ax2.plot(cumulative_count.index, cumulative_count["poly_id"])
    ax2.set_ylabel("cumulative data count")

    plt.savefig("assets/publications.svg")


async def draw():
    client = AsyncIOMotorClient(CLIENT_URL)
    database = client["roppy"]

    await init_beanie(
        database=database,
        document_models=[
            PolymerisationDocument,
            MonomerSummaryDocument,
        ],
    )

    await draw_publications()


if __name__ == "__main__":
    asyncio.run(draw())
