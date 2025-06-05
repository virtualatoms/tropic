import erdantic as erd
from roppy.core.models import MonomerSummary, Polymerisation

# TODO: Split into separate files for each figure

if __name__ == "__main__":

    # Create Erdantic Models
    poly_diagram_base = erd.create(Polymerisation, terminal_models=[Polymerisation])
    poly_diagram = erd.create(Polymerisation)
    monomer_diagram = erd.create(MonomerSummary, terminal_models=[Polymerisation])

    # remove id and revision_id
    diagrams = [poly_diagram_base, monomer_diagram, poly_diagram]
    fields = ["Polymerisation", "Monomer", "Initiator", "MonomerSummary"]
    for diagram in diagrams:
        for field in fields:
            try:
                del diagram.models[f"roppy.core.models.{field}"].fields["id"]
                del diagram.models[f"roppy.core.models.{field}"].fields["revision_id"]
            except:
                pass

    # separate polymerisation from its attributes
    del poly_diagram.models["roppy.core.models.Polymerisation"]
    poly_diagram.edges.clear()

    # tidy up monomerSummary graph
    del monomer_diagram.models["roppy.core.models.Monomer"]
    del monomer_diagram.edges[
        "roppy.core.models.MonomerSummary-monomer-roppy.core.models.Monomer"
    ]
    del monomer_diagram.models["roppy.core.models.DataRow"]
    del monomer_diagram.edges[
        "roppy.core.models.MonomerSummary-data-roppy.core.models.DataRow"
    ]

    # TODO:Sorting fields

    # Output diagrams
    poly_diagram_base.draw(
        "assets/Polymerisation_base.svg",
        graph_attr={
            "layout": "dot",
            "center": True,
            "normalise": 0,
            # "dpi": 300,
        },
        node_attr={"fontsize": 8},
    )
    poly_diagram.draw(
        "assets/Polymerisation_attr.svg",
        graph_attr={
            "layout": "dot",
            "rankdir": "TB",
            "ratio": 0.5,
            # "pack": 1,
            "packmode": "node",
            # "dpi": 300,
        },
        node_attr={"fontsize": 8},
    )
    monomer_diagram.draw(
        "assets/MonomerSummary.svg",
        graph_attr={
            "layout": "circo",
            "rankdir": "LR",
            "packmode": "node",
            # "dpi": 300,
        },
        node_attr={"fontsize": 8},
    )
