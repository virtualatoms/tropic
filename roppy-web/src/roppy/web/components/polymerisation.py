import dash_mantine_components as dmc
from roppy.web.utils import smiles_to_image
from roppy.web.components.table import get_table
from dash import html


def get_polymerisation_section(data):
    reaction_data = smiles_to_image(data["smiles"], size=(500, 200))
    reaction_img = dmc.Image(
        radius="md", h=200, src=f"data:image/svg+xml;base64,{reaction_data}"
    )

    poly_info = data["polymerisation"]
    poly = dmc.Stack(
        [
            html.H1(
                "Polymerisation", id="poly", style={"margin-top": 0, "padding-top": 00}
            ),
            dmc.Grid(
                [
                    dmc.GridCol(
                        dmc.Card(
                            reaction_img,
                            withBorder=True,
                            shadow="sm",
                            radius="md",
                        ),
                        span=8,
                    ),
                    dmc.GridCol(
                        dmc.Card(
                            [get_table([("a", "b")])],
                            withBorder=True,
                            shadow="sm",
                            radius="md",
                        ),
                        span=4,
                    ),
                ],
                gutter="xl",
            ),
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
