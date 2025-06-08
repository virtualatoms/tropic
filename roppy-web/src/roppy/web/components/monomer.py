import dash_mantine_components as dmc
from dash import html

from roppy.web.components.molview import get_mol_viewer_card
from roppy.web.components.table import get_table
from roppy.web.components.toc import create_toc_sidebar
from roppy.web.utils import smiles_to_image

HIDDEN = {"height": 0, "width": 0, "overflow": "hidden"}


def get_monomer_info_card(data):
    pubchem_link = (
        html.A(
            data["monomer"]["pubchem_cid"],
            href=f"https://pubchem.ncbi.nlm.nih.gov/compound/{data['monomer']['pubchem_cid']}",
            className="mantine-focus-auto m_849cf0da mantine-Text-root mantine-Anchor-root",
            style={"text-decoration": "underline"},
            target="_blank",
        )
        if data["monomer"]["pubchem_cid"]
        else "N/A"
    )
    table_data = [
        ("SMILES", data["monomer"]["smiles"]),
        ("InChI", data["monomer"]["inchi"]),
        ("IUPAC Name", data["monomer"]["iupac_name"]),
        # ("Common Name", data["monomer"]["common_name"]),
        ("Molecular Weight", f"{data['monomer']['molecular_weight']:.2f} g/mol"),
        ("Ring Size", data["monomer"]["ring_size"]),
        ("CID", pubchem_link),
    ]

    header = dmc.CardSection(
        dmc.Center([dmc.Text("Monomer Information", fw=700)]),
        withBorder=True,
        inheritPadding=True,
        py="xs",
        style={"backgroundColor": dmc.DEFAULT_THEME["colors"]["gray"][0]},
    )

    return dmc.Card(
        [header, get_table(table_data)], withBorder=True, shadow="sm", radius="md"
    )


def get_monomer_page_summary(data):
    viewer_card = get_mol_viewer_card(data["monomer"]["xyz"])
    info_card = get_monomer_info_card(data)

    return dmc.Grid(
        [
            html.H1("Summary", id="summary", style=HIDDEN),
            dmc.GridCol(viewer_card, span=6),
            dmc.GridCol(info_card, span=6),
        ],
        gutter="xl",
        pb=10,
    )


def get_monomer_logo(smiles: str, monomer_id: str):
    img_data = smiles_to_image(smiles, size=(-1, 120))
    img = dmc.Image(radius="md", src=f"data:image/svg+xml;base64,{img_data}")
    card = dmc.Card(
        [img],
        withBorder=True,
        radius="lg",
        style={
            "border-color": dmc.DEFAULT_THEME["colors"]["blue"][5],
            "border-width": 3,
        },
    )
    badge = dmc.Center(
        dmc.Badge(monomer_id, size="lg", radius="lg", mt=-35, style={"zIndex": 1})
    )
    return dmc.Box(dmc.Stack([card, badge]), mb=20, px=0, w="90%")


def get_monomer_toc(data):
    header = get_monomer_logo(data["monomer"]["smiles"], data["monomer_id"])
    return create_toc_sidebar(header)
