import dash
from dash import html
import dash_mantine_components as dmc
from roppy.web.utils import smiles_to_image
from roppy.web.header import get_breadcrumbs

dash.register_page(__name__, path_template="/monomers/<monomer_id>")


def layout(monomer_id=None, **kwargs):
    breadcrumbs = get_breadcrumbs(["Home", "Monomer Search", f"{monomer_id}"])

    img_data = smiles_to_image("C=O", size=(200, 200))
    img = dmc.Image(
        radius="md",
        h=200,
        src=f"data:image/svg+xml;base64,{img_data}",
    )
    monomer_card1 = dmc.Card(
        [img],
        withBorder=True,
        shadow="sm",
        radius="md",
    )
    monomer_card2 = dmc.Card(
        [img],
        withBorder=True,
        shadow="sm",
        radius="md",
    )
    monomer_card3 = dmc.Card(
        [img],
        withBorder=True,
        shadow="sm",
        radius="md",
    )

    # dmc.Card(
    #     children=[
    #         dmc.CardTitle("Monomer"),
    #         dmc.CardContent(
    #             children=[
    #                 dmc.Text("Monomer ID:"),
    #                 dmc.Text("SMILES:"),
    #                 dmc.Text("InChI:"),
    #                 dmc.Text("InChI Key:"),
    #                 dmc.Text("IUPAC Name:"),
    #                 dmc.Text("Common Name:"),
    #             ]
    #         ),
    #     ]
    # )

    grid = dmc.Grid(
        [
            dmc.GridCol(monomer_card1, span=3),
            dmc.GridCol(monomer_card2, span=3),
            dmc.GridCol(monomer_card3, span=3),
        ],
        pt=50,
        gutter="xl",
    )
    return [breadcrumbs, grid]
    # return html.Div(
    #     [
    #         html.H1(f"Monomers page {monomer_id}"),
    #     ]
    # )
