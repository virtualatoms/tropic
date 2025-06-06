import asyncio
from motor.motor_asyncio import AsyncIOMotorClient
from roppy.core.models import Polymerisation
from roppy.core.build import poly_id_to_poly

DATABASE_URL = "mongodb://localhost:27017"
DATABASE_NAME = "roppy"


def startup_event():
    client = AsyncIOMotorClient(DATABASE_URL)
    database = client[DATABASE_NAME]


async def insert_polymers(poly_id_to_poly: dict[int, Polymerisation]):
    result = await database.polymerisations.insert_many(poly_id_to_poly)
    print(result.inserted_ids)


if __name__ == "__main__":
    client = AsyncIOMotorClient(DATABASE_URL)
    database = client[DATABASE_NAME]

    loop = client.get_io_loop()
    loop.run_until_complete(insert_polymers)

    i
