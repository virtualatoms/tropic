import requests
from dash import html, register_page
import dash_mantine_components as dmc
from roppy.web.components.breadcrumbs import get_breadcrumbs
from roppy.web.components.toc import register_toc_callbacks
from roppy.web.components.molview import register_molview_callbacks
from roppy.web.components.monomer import get_monomer_toc, get_monomer_page_summary
from roppy.web.components.polymerisation import get_polymerisation_section


register_page(__name__, path_template="/monomers/<monomer_id>")

register_toc_callbacks()
register_molview_callbacks()


def layout(monomer_id="monomer-1", **kwargs):
    breadcrumbs = get_breadcrumbs(["Home", "Monomer Search", f"{monomer_id}"])

    response = requests.get(f"http://localhost:8000/monomers/{monomer_id}")
    data = response.json()

    if data.get("detail", "") == "Monomer not found":
        return html.H1("Monomer not found")

    toc = get_monomer_toc(data)
    monomer_summary = get_monomer_page_summary(data)

    page = [monomer_summary]

    # if data["polymerisation"]:
    #     poly_section = get_polymerisation_section(data["polymerisation"])
    #     page.extend([dmc.Divider(mt=60, mb=40), poly_section])

    # if data["ring_opening"]:
    #     rop_section = get_rop_section(data)
    #     page.extend([dmc.Divider(mt=60, mb=40), rop_section])

    content = dmc.Grid(
        [dmc.GridCol(toc, span=3), dmc.GridCol(html.Div(page, id="content"), span=9)],
        pt=40,
        gutter="xl",
    )

    return [breadcrumbs, content]
