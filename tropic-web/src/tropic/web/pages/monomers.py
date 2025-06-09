from urllib.parse import quote

import dash_mantine_components as dmc
import numpy as np
import requests
from dash import (
    Input,
    Output,
    State,
    callback,
    clientside_callback,
    html,
    no_update,
    register_page,
)
from dash_iconify import DashIconify

from tropic.web import SETTINGS
from tropic.web.components.breadcrumbs import get_breadcrumbs
from tropic.web.components.chart import get_search_chart
from tropic.web.components.draw import get_draw_molecule
from tropic.web.components.references import get_reference_table_data, get_references
from tropic.web.components.searchtable import SEARCH_NUM_ROWS, get_search_table
from tropic.web.components.sidebar import get_search_sidebar
from tropic.web.utils import smiles_to_image

register_page(__name__)

FILTERS = (
    ("smiles-input", "value"),
    ("ring-size-slider", "value"),
    ("mol-weight-slider", "value"),
    ("has-comp", "value"),
    ("has-exp", "value"),
)
FILTER_INPUTS = [Input(k, v) for k, v in FILTERS]
FILTER_STATES = [State(k, v) for k, v in FILTERS]


def layout(**_):
    breadcrumbs = get_breadcrumbs(["Home", "Monomer Search"])
    table = get_search_table()
    chart = get_search_chart()
    side_bar = get_search_sidebar()
    draw_modal = get_draw_molecule()
    references = get_references()

    tabs = dmc.Tabs(
        [
            dmc.TabsList(
                [
                    dmc.TabsTab(
                        "Table",
                        value="table",
                        leftSection=DashIconify(icon="tabler:table"),
                    ),
                    dmc.TabsTab(
                        "Analysis",
                        value="analysis",
                        leftSection=DashIconify(icon="tabler:chart-dots"),
                    ),
                    dmc.TabsTab(
                        "References",
                        value="references",
                        leftSection=DashIconify(icon="tabler:book"),
                    ),
                ]
            ),
            dmc.TabsPanel(table, value="table", pt=10),
            dmc.TabsPanel(chart, value="analysis", pt=10),
            dmc.TabsPanel(references, value="references", pt=10),
        ],
        id="tabs",
        value="table",
    )

    grid = dmc.Grid(
        [dmc.GridCol(side_bar, span=3, pt=108), dmc.GridCol(tabs, span=9)],
        pt=30,
        gutter="xl",
    )

    input_bar = dmc.Group(
        [
            dmc.TextInput(id="smiles-input", placeholder="Enter SMILES", w=500),
            dmc.Button("Draw Molecule", id="draw-button"),
        ],
        pt=30,
    )

    return [
        breadcrumbs,
        dmc.Center(input_bar),
        draw_modal,
        html.Div(id="table-trigger", **{"data-label": {}}),
        grid,
    ]


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
    *FILTER_INPUTS,
)
def reset_table(*_):
    return {"startRow": 1, "endRow": 5, "sortModel": [], "filterModel": {}}, []


@callback(
    Output("results-table", "getRowsResponse"),
    Input("tabs", "value"),
    Input("results-table", "getRowsRequest"),
    *FILTER_INPUTS,
)
def update_table(tabs, rows_request, *filter_args):
    if tabs != "table":
        return no_update

    query = [_build_query(*filter_args)]

    # handle pagination
    page = rows_request["startRow"] // SEARCH_NUM_ROWS + 1 if rows_request else 1
    query.append(f"size={SEARCH_NUM_ROWS}&page={page}")

    # handle sort
    if rows_request and rows_request["sortModel"]:
        sort_col = rows_request["sortModel"][0]["colId"]
        sort_dir = "-" if rows_request["sortModel"][0]["sort"] == "asc" else ""
        query.append(f"order_by={sort_dir}{sort_col}")

    query = "&".join(query)
    response = requests.get(
        f"{SETTINGS.API_ENDPOINT}/monomers?{query}", timeout=SETTINGS.REQUEST_TIMEOUT
    )
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


@callback(Output("analysis-chart", "data"), Input("tabs", "value"), *FILTER_INPUTS)
def update_chart(tabs, *filter_args):
    # TODO: Currently this repeats the API call from the table update.
    # Ideally, we should refactor to avoid duplicate requests.
    # currently one difference is that the chart wants all the data, not just the first page

    if tabs != "analysis":
        return no_update

    query = _build_query(*filter_args)
    query += "&size=1000"
    response = requests.get(f"{SETTINGS.API_ENDPOINT}/monomers?{query}", timeout=2)
    results = response.json()

    exp = []
    comp = []
    for result in results["items"]:
        for row in result["data"]:
            if row["delta_s"] and row["delta_h"]:
                delta_g = row["delta_h"] - row["delta_s"] * 298.15
            else:
                delta_g = None

            if row["is_experimental"]:
                exp.append(
                    {
                        "ring_size": result["monomer"]["ring_size"],
                        "delta_h": row["delta_h"],
                        "delta_s": row["delta_s"],
                        "delta_g": delta_g,
                        "ceiling_temperature": row["ceiling_temperature"],
                    }
                )
            else:
                comp.append(
                    {
                        "ring_size": result["monomer"]["ring_size"],
                        "delta_h": row["delta_h"],
                        "delta_s": row["delta_s"],
                        "delta_g": delta_g,
                        "ceiling_temperature": row["ceiling_temperature"],
                    }
                )

    return [
        {"color": "red.5", "name": "Computational", "data": comp},
        {"color": "blue.5", "name": "Experimental", "data": exp},
    ]


@callback(
    Output("references-table", "children"), Input("tabs", "value"), *FILTER_INPUTS
)
def update_references(tabs, *filter_args):
    if tabs != "references":
        return no_update

    query = _build_query(*filter_args)
    query += "&size=1000"
    response = requests.get(f"{SETTINGS.API_ENDPOINT}/monomers?{query}", timeout=2)
    results = response.json()

    refs = {}
    for result in results["items"]:
        for row in result["data"]:
            if row["doi"] not in refs and row["doi"] and row["formatted_reference"]:
                refs[row["doi"]] = row["formatted_reference"]

    return get_reference_table_data(*zip(*refs.items()))


@callback(
    Output("download-data", "data"),
    Input("export-data-button", "n_clicks"),
    State("export-data-select", "value"),
    *FILTER_STATES,
    prevent_initial_call=True,
)
def export(n_clicks, file_type, *filter_args):
    if n_clicks is None:
        return no_update

    query = _build_query(*filter_args)
    response = requests.get(
        f"{SETTINGS.API_ENDPOINT}/monomers?{query}", timeout=SETTINGS.REQUEST_TIMEOUT
    )
    if response.status_code != 200:
        return no_update

    data = get_export_data(response.json(), file_type)

    return {
        "content": data,
        "filename": f"monomers.{file_type}",
        "type": f"text/{file_type}",
    }


def get_export_data(data, file_type):
    if file_type == "csv":
        import csv
        from io import StringIO

        output = StringIO()
        writer = csv.writer(output)
        writer.writerow(
            [
                "Monomer ID",
                "SMILES",
                "Ring Size",
                "Type",
                "Is Experimental",
                "State",
                "Is Experimental",
                "Initial Monomer Conc",
                "Bulk Monomer Conc",
                "Solvent",
                "ΔH (kJ/mol)",
                "ΔS (kJ/mol·K)",
                "Ceiling Temperature (K)",
                "Repeating Units",
                "DFT Method",
                "Year",
                "Reference",
            ]
        )
        for item in data["items"]:
            for row in item["data"]:
                writer.writerow(
                    [
                        item["monomer_id"],
                        item["monomer"]["smiles"],
                        item["monomer"]["ring_size"],
                        row["type"],
                        row["is_experimental"],
                        row["state_summary"],
                        row["is_experimental"],
                        row["initial_monomer_conc"],
                        row["bulk_monomer_conc"],
                        row["solvent"],
                        row["delta_h"],
                        row["delta_s"],
                        row["ceiling_temperature"],
                        row["repeating_units"],
                        row["method"],
                        row["year"],
                        row["formatted_reference"],
                    ]
                )
        return output.getvalue()

    if file_type == "json":
        import json

        return json.dumps(data, indent=4)

    return ""


@callback(
    Output("analysis-chart", "dataKey"),
    Output("analysis-chart", "yAxisLabel"),
    Input("chart-select-y", "value"),
)
def update_chart_y_axis(selected_y):
    data_key = {"x": "ring_size", "y": selected_y} if selected_y else no_update
    y_axis_labels = {
        "delta_h": r"ΔH / kj/mol",
        "delta_s": r"ΔS / kj/mol",
        "delta_g": r"ΔG / kj/mol",
        "ceiling_temperature": "Ceiling Temperature / K",
    }
    y_axis_label = y_axis_labels.get(selected_y, "Value")
    return data_key, y_axis_label


def _build_query(smiles, ring_size_range, molecular_weight_range, has_comp, has_exp):
    query = []
    if smiles:
        query.append(f"search={quote(smiles)}")
    query.append(f"monomer__ring_size__gte={ring_size_range[0]}")

    query.append(f"monomer__molecular_weight__gte={molecular_weight_range[0]}")
    query.append(f"monomer__molecular_weight__lte={molecular_weight_range[1]}")

    if ring_size_range[1] is not None and ring_size_range[1] < 15:
        query.append(f"monomer__ring_size__lte={ring_size_range[1]}")

    if has_comp != "both":
        query.append(f"has_calc={has_comp == 'yes'}")

    if has_exp != "both":
        query.append(f"has_exp={has_exp == 'yes'}")

    return "&".join(query)


# this hack is needed to make AG Grid work properly

clientside_callback(
    """async () => {
        const gridApi = await dash_ag_grid.getApiAsync("results-table");
        gridApi.refreshInfiniteCache();
        return dash_clientside.no_update
    }""",
    Output("table-trigger", "data-label"),
    *FILTER_INPUTS,
)
