# Polymerisation summary objects are input through the template given to experimentalists

# from a list of Polymerisation objects generate a list of MonomerSummary objects
# which each contain polymermisation summary objects

from roppy.core.models import *
from typing import Optional, Any

# data format extracted from database
data: list[dict[str, Any]] = [
    {
        "monomer_smiles": "CCN(CC)C(=O)[C@@H]1C=C2c3cccc4[nH]cc(c34)C[C@H]2N(C)C1",
        "initiator_smiles": "CNC(C)Cc1ccc2c(c1)OCO2",
        "product": {
            "number_of_units": 1.0,
        },
        "is_experimental": True,
        "parameters": {
            "temperature": None,
            "pressure": None,
            "solvent": None,
            "solvent_conc": None,
            "monomer_state": None,
            "polymer_state": None,
            "method": None,
            "solvent_model": None,
            "functional": None,
            "basis_set": None,
            "dispersion": None,
            "forcefield": None,
        },
        "thermo": {
            "delta_h": None,
            "delta_s": None,
            "ceiling_temperature": None,
        },
        "comment": None,
        "doi": None,
        "url": None,
    }
]


# building hashmaps from database
def build_from_db(
    polymerisation_data: list[dict[str, Any]],
) -> tuple[dict[int, Polymerisation], dict[int, Monomer], dict[int, Initiator]]:

    monomer_smiles_to_monomer_id: dict[str, int] = dict()

    initiator_smiles_to_initiator_id: dict[str, int] = dict()

    monomer_id_to_poly_id: dict[int, list[int]] = dict()

    initiator_id_to_poly_id: dict[int, list[int]] = dict()

    # iterate to populate hashmaps
    for poly_id, data in enumerate(polymerisation_data):

        # map monomer smiles to monomer id
        if data["monomer_smiles"] not in monomer_smiles_to_monomer_id:
            monomer_smiles_to_monomer_id[data["monomer_smiles"]] = poly_id

        # map monomer id to polymerisations
        monomer_id = monomer_smiles_to_monomer_id[data["monomer_smiles"]]
        if monomer_id in monomer_id_to_poly_id:
            monomer_id_to_poly_id[monomer_id].append(poly_id)
        else:
            monomer_id_to_poly_id[monomer_id] = [poly_id]

        # if polymerisation has an initiator
        if data["initiator_smiles"]:

            # map initiator smiles to initiator id
            if data["initiator_smiles"] not in initiator_smiles_to_initiator_id:
                initiator_smiles_to_initiator_id[data["initiator_smiles"]] = poly_id

            # map initiator id to polymerisations
            initiator_id = initiator_smiles_to_initiator_id[data["initiator_smiles"]]
            if initiator_id in initiator_id_to_poly_id:
                initiator_id_to_poly_id[initiator_id].append(poly_id)
            else:
                initiator_id_to_poly_id[initiator_id] = [poly_id]

    print(monomer_smiles_to_monomer_id)
    print(initiator_smiles_to_initiator_id)
    print(monomer_id_to_poly_id)
    print(initiator_id_to_poly_id)

    poly_id_to_poly: dict[int, Polymerisation] = dict()

    monomer_id_to_monomer: dict[int, Monomer] = dict()

    initiator_id_to_initiator: dict[int, Initiator] = dict()

    for poly_id, data in enumerate(polymerisation_data):

        monomer_id = monomer_smiles_to_monomer_id[data["monomer_smiles"]]
        if data["initiator_smiles"]:
            initiator_id = initiator_smiles_to_initiator_id[data["initiator_smiles"]]
        else:
            initiator_id = None

        poly = Polymerisation(
            polymerisation_id=poly_id,
            product=Product(**data["product"]),
            is_experimental=data["is_experimental"],
            parameters=Parameters(**data["parameters"]),
            thermo=Thermo(**data["thermo"]),
            comment=data["comment"],
            doi=data["doi"],
            url=data["url"],
            monomer_id=monomer_id,
            initiator_id=initiator_id,
        )

        poly_id_to_poly[poly_id] = poly

    for monomer_smiles, monomer_id in monomer_smiles_to_monomer_id.items():

        monomer = Monomer(
            monomer_id=monomer_id,
            smiles=monomer_smiles,
            polymerisations=monomer_id_to_poly_id[monomer_id],
        )

        monomer_id_to_monomer[monomer_id] = monomer

    for initiator_smiles, initiator_id in initiator_smiles_to_initiator_id.items():

        initiator = Initiator(
            initiator_id=initiator_id,
            smiles=initiator_smiles,
            polymerisations=initiator_id_to_poly_id[initiator_id],
        )

        initiator_id_to_initiator[initiator_id] = initiator

    print(poly_id_to_poly)
    print(monomer_id_to_monomer)
    print(initiator_id_to_initiator)

    return poly_id_to_poly, monomer_id_to_monomer, initiator_id_to_initiator


if __name__ == "__main__":

    a, b, c = build_from_db(data)
