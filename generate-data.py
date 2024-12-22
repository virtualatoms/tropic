from roppy.api.documents import MonomerSummaryDocument as MonomerSummary
from beanie import init_beanie
import motor.motor_asyncio


async def init():
    client = motor.motor_asyncio.AsyncIOMotorClient("mongodb://localhost:27017")
    await init_beanie(database=client.roppy, document_models=[MonomerSummary])

    await MonomerSummary(
        monomer_id="1",
        smiles="CCS=O",
        ring_size=2,
        has_exp=False,
        has_calc=True,
    ).save()
    await MonomerSummary(
        monomer_id="2",
        smiles="CC=O",
        ring_size=5,
        has_exp=True,
        has_calc=True,
    ).save()
    await MonomerSummary(
        monomer_id="3",
        smiles="CCC1CCC(=O)O1",
        ring_size=2,
        has_exp=False,
        has_calc=False,
    ).save()
    await MonomerSummary(
        monomer_id="4",
        smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        ring_size=6,
        has_exp=True,
        has_calc=False,
    ).save()


if __name__ == "__main__":
    import asyncio

    asyncio.run(init())
