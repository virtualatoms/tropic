import requests
import numpy as np
import dash_mantine_components as dmc
from urllib.parse import quote
from dash import Input, Output, State, callback, no_update, register_page, html
from dash import clientside_callback
from dash_iconify import DashIconify
from roppy.web.utils import smiles_to_image
from roppy.web.components.breadcrumbs import get_breadcrumbs
from roppy.web.components.sidebar import get_search_sidebar
from roppy.web.components.chart import get_search_chart
from roppy.web.components.searchtable import get_search_table, PAGE_SIZE
from roppy.web.components.draw import get_draw_molecule   


register_page(__name__)


def layout(**kwargs):
    breadcrumbs = get_breadcrumbs(["Home", "Monomer Search"])
    table = get_search_table()
    chart = get_search_chart()
    side_bar = get_search_sidebar()
    draw_modal = get_draw_molecule()
    
    tabs = dmc.Tabs(
        [
            dmc.TabsList(
                [
                    dmc.TabsTab(
                        "Table", value="table", leftSection=DashIconify(icon="tabler:table")
                    ),
                    dmc.TabsTab(
                        "Analysis",
                        value="analysis",
                        leftSection=DashIconify(icon="tabler:chart-dots"),
                    ),
                ]
            ),
            dmc.TabsPanel(table, value="table", pt=20),
            dmc.TabsPanel(chart, value="analysis", pt=20),
        ],
        value="table",
    )

    grid = dmc.Grid(
        [dmc.GridCol(side_bar, span=3, pt=72), dmc.GridCol(tabs, span=9)],
        pt=50,
        gutter="xl",
    )

    input_bar = dmc.Group(
        [
            dmc.TextInput(id="smiles-input", placeholder="Enter SMILES", w=500),
            dmc.Button("Draw Molecule", id="draw-button"),
        ],
        pt=50,
    )

    layout = [
        breadcrumbs,
        dmc.Center(input_bar),
        grid,
        draw_modal,
        html.Div(id='table-trigger', **{'data-label': {}}),
    ]

    return layout
    
    
@callback(
    Output("molecule-drawer-modal", "opened"),
    [
        Input("draw-button", "n_clicks"),
        Input("close-drawer", "n_clicks"),
        Input("apply-drawn-molecule", "n_clicks"),
    ],
    [State("molecule-drawer-modal", "opened")],
    prevent_initial_call=True,
)
def toggle_modal(draw_clicks, close_clicks, apply_clicks, is_open):
    if draw_clicks or close_clicks or apply_clicks:
        return not is_open
    return is_open

    
@callback(
    Output("smiles-input", "value"),
    [Input("apply-drawn-molecule", "n_clicks")],
    [State("jsme", "eventSmiles")],
    prevent_initial_call=True,
)
def update_smiles_from_drawing(n_clicks, drawn_smiles):
    if n_clicks and drawn_smiles:
        return drawn_smiles
    return no_update
    
@callback(
    Output("results-table", "getRowsRequest"),
    Output("results-table", "rowData"),
    Input("smiles-input", "value"),
    Input("ring-size-slider", "value"),
    Input("has-comp", "value"),
    Input("has-exp", "value"),
)
def reset_table(*_):
    return {"startRow": 1, "endRow": 5, "sortModel": [], "filterModel": {}}, []


@callback(
    Output("results-table", "getRowsResponse"),
    Input("smiles-input", "value"),
    Input("ring-size-slider", "value"),
    Input("has-comp", "value"),
    Input("has-exp", "value"),
    Input("results-table", "getRowsRequest"),
)
def update_table(smiles, ring_size_range, has_comp, has_exp, rows_request):
    # if rows_request is None:
    #     return no_update

    query = []
    if smiles:
        query.append(f"search={quote(smiles)}")
    query.append(f"monomer__ring_size__gte={ring_size_range[0]}")
    query.append(f"monomer__ring_size__lte={ring_size_range[1]}")

    if has_comp != "both":
        query.append(f"has_calc={has_comp == 'yes'}")

    if has_exp != "both":
        query.append(f"has_exp={has_exp == 'yes'}")

    # handle pagination
    page = rows_request["startRow"] // PAGE_SIZE + 1 if rows_request else 1
    query.append(f"size=5&page={page}")

    # handle sort
    if rows_request and rows_request["sortModel"]:
        sort_col = rows_request["sortModel"][0]["colId"]
        sort_dir = "-" if rows_request["sortModel"][0]["sort"] == "asc" else ""
        query.append(f"order_by={sort_dir}{sort_col}")

    query = "&".join(query)
    response = requests.get(f"http://localhost:8000/monomers?{query}")
    results = response.json()

    yes_no_mapping = {True: "Yes", False: "No"}
    table_data = []
    for result in results["items"]:
        img_data = smiles_to_image(result["monomer"]["smiles"])
        structure_cell = (
            f'<img src="data:image/svg+xml;base64,{img_data}" style="max-width:150px; display: block; margin: auto;">'
            if img_data
            else ""
        )

        table_data.append(
            {
                "structure": structure_cell,
                "monomer_id": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="monomers/{result["monomer_id"]}">{result["monomer_id"]}</a>',
                "smiles": result["monomer"]["smiles"],
                "ring_size": result["monomer"]["ring_size"],
                "has_exp": yes_no_mapping[result["has_exp"]],
                "has_calc": yes_no_mapping[result["has_calc"]],
            }
        )

    return {"rowData": np.array(table_data), "rowCount": results["total"]}



# this hack is needed to make AG Grid work properly 

clientside_callback(
    """async () => {
        const gridApi = await dash_ag_grid.getApiAsync("results-table");
        gridApi.refreshInfiniteCache();
        return dash_clientside.no_update
    }""",
    Output('table-trigger', 'data-label'),
    Input("smiles-input", "value"),
    Input("ring-size-slider", "value"),
    Input("has-comp", "value"),
    Input("has-exp", "value"),
)