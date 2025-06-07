import dash_mantine_components as dmc
import requests
from dash import html, register_page

from roppy.web import SETTINGS
from roppy.web.components.breadcrumbs import get_breadcrumbs
from roppy.web.components.molview import register_molview_callbacks
from roppy.web.components.monomer import get_monomer_page_summary, get_monomer_toc
from roppy.web.components.polymerisation import get_polymerisation_section
from roppy.web.components.toc import register_toc_callbacks

register_page(__name__, path_template="/monomers/<monomer_id>")

register_toc_callbacks()
register_molview_callbacks()


def layout(monomer_id="monomer-1", **_):
    breadcrumbs = get_breadcrumbs(["Home", "Monomer Search", f"{monomer_id}"])

    response = requests.get(
        f"{SETTINGS.API_ENDPOINT}/monomers/{monomer_id}",
        timeout=SETTINGS.REQUEST_TIMEOUT,
    )
    data = response.json()

    if data.get("detail", "") == "Monomer not found":
        return html.H1("Monomer not found")

    toc = get_monomer_toc(data)
    monomer_summary = get_monomer_page_summary(data)

    polymerisation_section = get_polymerisation_section(data)
    page = [monomer_summary, polymerisation_section]

    # if data["has_exp"]:
    #     poly_section = get_polymerisation_section(data["data"], exp=True)
    #     page.extend([dmc.Divider(mt=60, mb=40), poly_section])

    # if data["has_comp"]:
    #     poly_section = get_polymerisation_section(data["data"], exp=False)
    #     page.extend([dmc.Divider(mt=60, mb=40), poly_section])

    content = dmc.Grid(
        [dmc.GridCol(toc, span=3), dmc.GridCol(html.Div(page, id="content"), span=9)],
        pt=60,
        gutter="xl",
        mb=50,
    )

    return [breadcrumbs, content]
