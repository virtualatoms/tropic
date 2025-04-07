import dash_mantine_components as dmc
from dash import dcc
from roppy.web.utils import smiles_to_image
from roppy.web.components.table import get_table
from dash import html


def get_polymer_info_card(polyinfo):
    table_data = [
        ("PolyInfo", dcc.Link(polyinfo["polyinfo_id"], href=f"https://polyinfo.com/{polyinfo['polyinfo_id']}")),
        ("PolyGenome", dcc.Link(polyinfo["polygenome_id"], href=f"https://polygenome.com/{polyinfo['polygenome_id']}")),
    ]
    return dmc.Card([get_table(table_data)], withBorder=True, shadow="sm", radius="md")


def get_polymer_reaction_card(data):
    reaction_data = smiles_to_image(data["reaction_smiles"], size=(500, 200))
    reaction_img = dmc.Image(
        radius="md", h=200, src=f"data:image/svg+xml;base64,{reaction_data}"
    )
    return dmc.Card(reaction_img, withBorder=True, shadow="sm", radius="md")


def get_polymerisation_section(data):

    poly_info_card = get_polymer_info_card(data["polymer"])
    poly_reaction_card = get_polymer_reaction_card(data)

    summary = dmc.Grid(
        [
            dmc.GridCol(poly_reaction_card, span=8),
            dmc.GridCol(poly_info_card, span=4),
        ],
        gutter="xl",
    )

    poly = dmc.Stack(
        [
            html.H1(
                "Polymerisation", id="poly", style={"margin-top": 0, "padding-top": 00}
            ),
            summary,
            html.H2(
                "Experimental",
                id="poly-exp",
                style={"margin-bottom": 10, "padding-bottom": 0},
            ),
            html.H2("Computational", id="poly-comp"),
        ],
        gap=20,
    )
    return poly
