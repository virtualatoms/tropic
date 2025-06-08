import erdantic as erd
from tropic.core.models import MonomerSummary, Polymerisation

# TODO: Split into separate files for each figure

if __name__ == "__main__":
    # Create Erdantic Models
    poly_base_diagram = erd.create(Polymerisation, terminal_models=[Polymerisation])
    poly_attr_diagram = erd.create(Polymerisation)
    monomer_diagram = erd.create(MonomerSummary, terminal_models=[Polymerisation])

    # remove id and revision_id
    diagrams = [poly_base_diagram, monomer_diagram, poly_attr_diagram]
    fields = ["Polymerisation", "Monomer", "Initiator", "MonomerSummary"]
    for diagram in diagrams:
        for field in fields:
            try:
                del diagram.models[f"tropic.core.models.{field}"].fields["id"]
                del diagram.models[f"tropic.core.models.{field}"].fields["revision_id"]
            except IndexError:
                pass

    # separate polymerisation from its attributes
    del poly_attr_diagram.models["tropic.core.models.Polymerisation"]
    poly_attr_diagram.edges.clear()

    # tidy up monomerSummary graph
    del monomer_diagram.models["tropic.core.models.Monomer"]
    del monomer_diagram.edges[
        "tropic.core.models.MonomerSummary-monomer-tropic.core.models.Monomer"
    ]
    del monomer_diagram.models["tropic.core.models.DataRow"]
    del monomer_diagram.edges[
        "tropic.core.models.MonomerSummary-data-tropic.core.models.DataRow"
    ]

    # TODO:Sorting fields

    # Output diagrams
    poly_base_diagram.draw(
        "assets/Polymerisation_base.svg",
        graph_attr={
            "layout": "dot",
            "center": True,
            "normalise": 0,
            # "dpi": 300,
        },
        node_attr={"fontsize": 8},
    )
    poly_attr_diagram.draw(
        "assets/Polymerisation_attr.svg",
        graph_attr={
            "layout": "dot",
            "rankdir": "TB",
            "ratio": 0.5,
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
