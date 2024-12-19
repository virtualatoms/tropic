from documents import MonomerSummary
from beanie import init_beanie
import motor.motor_asyncio


async def init():
    client = motor.motor_asyncio.AsyncIOMotorClient("mongodb://localhost:27017")
    await init_beanie(database=client.roppy, document_models=[MonomerSummary])

    await MonomerSummary(
        monomer_id="1",
        smiles="CCS=O",
        ring_size=2,
        average_delta_h=-1.0,
        average_delta_s=-2.0,
    ).save()
    await MonomerSummary(
        monomer_id="2",
        smiles="CC=O",
        ring_size=5,
        average_delta_h=1.0,
        average_delta_s=2.0,
    ).save()
    await MonomerSummary(
        monomer_id="3",
        smiles="CCC1CCC(=O)O1",
        ring_size=2,
        average_delta_h=-20.0,
        average_delta_s=10.0,
    ).save()
    await MonomerSummary(
        monomer_id="4",
        smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        ring_size="many",
        average_delta_h=-22.0,
        average_delta_s=0.01,
    ).save()


if __name__ == "__main__":
    import asyncio

    loop = asyncio.new_event_loop()
    tasks = [init()]
    loop.run_until_complete(asyncio.wait(tasks))
    loop.close()
