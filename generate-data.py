from roppy.api.documents import MonomerSummaryDocument as MonomerSummary
from beanie import init_beanie
import motor.motor_asyncio


async def init():
    client = motor.motor_asyncio.AsyncIOMotorClient("mongodb://localhost:27017")
    await init_beanie(database=client.roppy, document_models=[MonomerSummary])

    await MonomerSummary(
        monomer_id="monomer-1",
        smiles="CN(C)CCc1c[nH]c2ccccc12",
        ring_size=2,
        has_exp=False,
        has_calc=True,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-2",
        smiles="c1ccc(C2(N3CCCCC3)CCCCC2)cc1",
        ring_size=5,
        has_exp=True,
        has_calc=True,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-3",
        smiles="CCC1CCC(=O)O1",
        ring_size=2,
        has_exp=False,
        has_calc=False,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-4",
        smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        ring_size=6,
        has_exp=True,
        has_calc=False,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-5",
        smiles="CNC(C)Cc1ccc2c(c1)OCO2",
        ring_size=3,
        has_exp=False,
        has_calc=True,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-6",
        smiles="COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",
        ring_size=6,
        has_exp=True,
        has_calc=True,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-7",
        smiles="CN[C@@H](C)Cc1ccccc1",
        ring_size=3,
        has_exp=False,
        has_calc=False,
    ).save()
    await MonomerSummary(
        monomer_id="monomer-8",
        smiles="CCN(CC)C(=O)[C@@H]1C=C2c3cccc4[nH]cc(c34)C[C@H]2N(C)C1",
        ring_size=4,
        has_exp=True,
        has_calc=False,
    ).save()


if __name__ == "__main__":
    import asyncio

    asyncio.run(init())
