from roppy.api.documents import MonomerSummaryDocument as MonomerSummary
from roppy.api.documents import PolymerisationSummaryDocument as PolymerisationSummary
from beanie import init_beanie
import motor.motor_asyncio

poly_data = {
    "reaction_smiles": "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",
    "polymer": {
        "polyinfo_id": "polyinfo_123",
        "polygenome_id": "polygenome_456",
    },
    "experimental_data": [
        {
            "state": "l-g",
            "solvent_conc": 0.5,
            "delta_h": -10.0,
            "delta_s": -20.0,
            "critical_temperature": 300.0,
            "equilibrium_temperature": 298.0,
            "doi": "10.1234/5678",
            "polymerisation_id": "poly_123",
        },
        {
            "state": "s-s",
            "solvent_conc": 2,
            "delta_h": -50.0,
            "delta_s": 20.0,
            "critical_temperature": 20.0,
            "equilibrium_temperature": 10.0,
            "doi": "10.1234/5678",
            "polymerisation_id": "poly_125",
        },
    ],
    "polymerisation_summary_id": "poly_summary_123",
}

rop_data = {
    "reaction_smiles": "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",
    "polymer": {
        "polyinfo_id": "polyinfo_123",
        "polygenome_id": "polygenome_456",
    },
    "experimental_data": [
        {
            "state": "lg",
            "solvent_conc": 0.5,
            "delta_h": -10.0,
            "delta_s": -20.0,
            "critical_temperature": 300.0,
            "equilibrium_temperature": 298.0,
            "doi": "10.1234/5678",
            "polymerisation_id": "poly_123",
            "initiator_smiles": "CC",
            "number_of_units": 10,
        },
        {
            "state": "ss",
            "solvent_conc": 2,
            "delta_h": -50.0,
            "delta_s": 20.0,
            "critical_temperature": 20.0,
            "equilibrium_temperature": 10.0,
            "doi": "10.1234/5678",
            "polymerisation_id": "poly_125",
            "initiator_smiles": "C=O",
            "number_of_units": 2,
        },
    ],
    "computational_data": [
        {
            "state": "solid",
            "solvent_model": "implicit",
            "solvent": "water",
            "delta_h": -10.0,
            "delta_s": -20.0,
            "critical_temperature": 300.0,
            "equilibrium_temperature": 298.0,
            "doi": "10.1234/5678",
            "polymerisation_id": "poly_123",
            "initiator_smiles": "CC",
            "number_of_units": 10,
        },
        {
            "state": "gas",
            "solvent_model": "implicit",
            "solvent": "dcm",
            "delta_h": -20.0,
            "delta_s": 10.0,
            "critical_temperature": 450.0,
            "equilibrium_temperature": 200.0,
            "doi": "10.1234/5672",
            "polymerisation_id": "poly_125",
            "initiator_smiles": "CC",
            "number_of_units": 1,
        },
    ],
    "polymerisation_summary_id": "poly_summary_123",
}
monomer_info = {
    "inchi": "InChI=1S/C20H26N2O5/c1-24-19(23)20(25-2)18(22)17-15-10-6-7-11-16(15)26-14-9-5-4-8-13(14)12-21(17)3/h6-7,10-11,13-14H,4-5,8-9,12H2,1-3H3/t13-,14+",
    "iupac_name": "methyl (2S,3S,4R,5R)-3,4-dihydroxy-5-((2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)tetrahydro-2H-pyran-2-yl)oxy-2-(hydroxymethyl)cyclohexyl (2S)-2-piperidinecarboxylate",
    "common_name": "Cocaine",
    "xyz": "12\nMethacrylic acid\nC\t-1.34\t-0.84\t0.14\nC\t-0.19\t0.09\t-0.01\nC\t-0.45\t1.38\t0.03\nC\t1.18\t-0.41\t-0.20\nO\t1.44\t-1.64\t-0.24\nO\t2.25\t0.47\t-0.34\nH\t-2.10\t-0.68\t-0.66\nH\t-1.88\t-0.58\t1.09\nH\t-1.04\t-1.90\t0.11\nH\t-1.44\t1.77\t0.16\nH\t0.37\t2.07\t-0.08\nH\t3.19\t0.25\t0.02\n",
}


async def init():
    client = motor.motor_asyncio.AsyncIOMotorClient("mongodb://localhost:27017")
    await init_beanie(
        database=client.roppy, document_models=[MonomerSummary, PolymerisationSummary]
    )

    ro = PolymerisationSummary(**rop_data)
    await ro.save()

    po = PolymerisationSummary(**poly_data)
    await po.save()

    await MonomerSummary(
        monomer_id="monomer-1",
        smiles="CN(C)CCc1c[nH]c2ccccc12",
        ring_size=2,
        has_exp=False,
        has_calc=True,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-2",
        smiles="c1ccc(C2(N3CCCCC3)CCCCC2)cc1",
        ring_size=5,
        has_exp=True,
        has_calc=True,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-3",
        smiles="CCC1CCC(=O)O1",
        ring_size=2,
        has_exp=False,
        has_calc=False,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-4",
        smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        ring_size=6,
        has_exp=True,
        has_calc=False,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-5",
        smiles="CNC(C)Cc1ccc2c(c1)OCO2",
        ring_size=3,
        has_exp=False,
        has_calc=True,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-6",
        smiles="COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",
        ring_size=6,
        has_exp=True,
        has_calc=True,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-7",
        smiles="CN[C@@H](C)Cc1ccccc1",
        ring_size=3,
        has_exp=False,
        has_calc=False,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()
    await MonomerSummary(
        monomer_id="monomer-8",
        smiles="CCN(CC)C(=O)[C@@H]1C=C2c3cccc4[nH]cc(c34)C[C@H]2N(C)C1",
        ring_size=4,
        has_exp=True,
        has_calc=False,
        monomer_info=monomer_info,
        polymerisation=po,
        ring_opening=[ro],
    ).save()


if __name__ == "__main__":
    import asyncio

    asyncio.run(init())
