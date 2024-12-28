import requests
import numpy as np
import dash_bio as dashbio
import dash_mantine_components as dmc
import dash_ag_grid as dag
from urllib.parse import quote
from dash import Input, Output, State, callback, no_update, register_page
from dash_iconify import DashIconify
from dash import dcc
from roppy.web.utils import smiles_to_image
from roppy.web.header import get_breadcrumbs

register_page(__name__, path="/monomers")

PAGE_SIZE = 5

data = [
    {
        "color": "blue.5",
        "name": "Experimental",
        "data": [
            {"ring_size": x, "delta_h": y, "delta_s": z}
            for x, y, z in zip(range(1, 11), range(10, 20), range(20, 30))
        ],
    },
    {
        "color": "red.5",
        "name": "Computational",
        "data": [
            {"ring_size": x, "delta_h": y, "delta_s": z}
            for x, y, z in zip(range(1, 11), range(20, 30), range(30, 40))
        ],
    },
]

chart = dmc.ScatterChart(
    h=500,
    data=data,
    dataKey={"x": "ring_size", "y": "delta_h"},
    xAxisLabel="Ring Size",
    yAxisLabel=r"ΔH",
    pointLabels="xy",
    tickLine="xy",
    mt=20,
)

table = dag.AgGrid(
    id="results-table",
    defaultColDef={
        "resizable": True,
        "cellStyle": {"align-items": "center", "height": "100%", "display": "flex"},
    },
    columnDefs=[
        {"headerName": "Structure", "field": "structure", "cellRenderer": "markdown"},
        {"headerName": "Monomer ID", "field": "monomer_id", "cellRenderer": "markdown"},
        {"headerName": "SMILES", "field": "smiles"},
        {"headerName": "Ring Size", "field": "ring_size"},
        {"headerName": "Has Exp", "field": "has_exp"},
        {"headerName": "Has Comp", "field": "has_calc"},
    ],
    columnSize="responsiveSizeToFit",
    columnSizeOptions={
        "columnLimits": [
            {"key": "structure", "minWidth": 200},
            {"key": "smiles", "minWidth": 250},
        ]
    },
    rowModelType="infinite",
    dashGridOptions={
        "animateRows": False,
        "domLayout": "autoHeight",
        "suppressCellFocus": True,
        "pagination": True,
        "paginationPageSizeSelector": False,
        "cacheBlockSize": PAGE_SIZE,
        "paginationPageSize": PAGE_SIZE,
        "rowHeight": 100,
    },
    style={"height": None},
    dangerously_allow_code=True,
)


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

radio_data = [["yes", "Yes"], ["no", "No"], ["both", "Both"]]
yes_no_mapping = {True: "Yes", False: "No"}

side_bar = dmc.Card(
    children=[
        dmc.CardSection(
            dmc.Text("Filters", fw=700),
            withBorder=True,
            inheritPadding=True,
            py="xs",
        ),
        dmc.Space(h=10),
        dmc.Text("Ring Size", size="sm", style={"font-weight": "500"}),
        dmc.RangeSlider(
            id="ring-size-slider",
            value=[1, 15],
            marks=[
                {"value": 1, "label": "1"},
                {"value": 5, "label": "5"},
                {"value": 10, "label": "10"},
                {"value": 15, "label": "15+"},
            ],
            min=1,
            max=15,
            minRange=1,
            mb=25,
            mt=10,
        ),
        dmc.Space(h=10),
        dmc.Text("Molecular Weight (g/mol)", size="sm", style={"font-weight": "500"}),
        dmc.RangeSlider(
            id="mol-weight-slider",
            value=[10, 1000],
            marks=[
                {"value": i, "label": str(i)} for i in [10, 200, 400, 600, 800, 1000]
            ],
            min=10,
            max=1000,
            minRange=10,
            mb=25,
            mt=10,
        ),
        dmc.Space(h=10),
        dmc.MultiSelect(
            label="Functional Groups",
            placeholder="Select functional groups",
            id="functional-groups",
            value=[],
            data=[
                {"value": "ester", "label": "Ester"},
                {"value": "aldehyde", "label": "Aldehyde"},
                {"value": "ketone", "label": "Ketone"},
                {"value": "amide", "label": "Amide"},
                {"value": "ether", "label": "Ether"},
                {"value": "amine", "label": "Amine"},
            ],
            # mb=10,
            clearable=True,
        ),
        dmc.Space(h=10),
        dmc.RadioGroup(
            children=dmc.Group([dmc.Radio(l, value=k) for k, l in radio_data], my=10),
            id="has-exp",
            value="both",
            label="Has Experimental Data",
            size="sm",
        ),
        dmc.Space(h=10),
        dmc.RadioGroup(
            children=dmc.Group([dmc.Radio(l, value=k) for k, l in radio_data], my=10),
            id="has-comp",
            value="both",
            label="Has Computational Data",
            size="sm",
        ),
    ],
    withBorder=True,
    shadow="sm",
    radius="md",
)

grid = dmc.Grid(
    [dmc.GridCol(side_bar, span=3, pt=72), dmc.GridCol(tabs, span=9)],
    pt=50,
    gutter="xl",
)

draw_modal = dmc.Modal(
    title="Draw Molecule",
    children=[
        dashbio.Jsme(id="jsme", smiles=""),
        dmc.Space(h=20),
        dmc.Group(
            [
                dmc.Button("Apply", id="apply-drawn-molecule"),
                dmc.Button(
                    "Cancel",
                    id="close-drawer",
                    color="red",
                    variant="outline",
                ),
            ],
            justify="flex-end",
        ),
    ],
    id="molecule-drawer-modal",
    size=650,
)

input_bar = dmc.Group(
    [
        dmc.TextInput(id="smiles-input", placeholder="Enter SMILES", w=500),
        dmc.Button("Draw Molecule", id="draw-button"),
    ],
    pt=50,
)

breadcrumbs = get_breadcrumbs(["Home", "Monomer Search"])

layout = [
    breadcrumbs,
    dmc.Center(input_bar),
    grid,
    draw_modal,
]


@callback(
    Output("results-table", "getRowsRequest"),
    Input("smiles-input", "value"),
    Input("ring-size-slider", "value"),
    Input("has-comp", "value"),
    Input("has-exp", "value"),
)
def reset_table(*_):
    return {"startRow": 1, "endRow": 5, "sortModel": [], "filterModel": {}}


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
    query.append(f"ring_size__gte={ring_size_range[0]}")
    query.append(f"ring_size__lte={ring_size_range[1]}")

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

    table_data = []
    for result in results["items"]:
        img_data = smiles_to_image(result["smiles"])
        structure_cell = (
            f'<img src="data:image/svg+xml;base64,{img_data}" style="max-width:150px; display: block; margin: auto;">'
            if img_data
            else ""
        )

        table_data.append(
            {
                "structure": structure_cell,
                "monomer_id": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="monomers/{result["monomer_id"]}">{result["monomer_id"]}</a>',
                "smiles": result["smiles"],
                "ring_size": result["ring_size"],
                "has_exp": yes_no_mapping[result["has_exp"]],
                "has_calc": yes_no_mapping[result["has_calc"]],
            }
        )

    return {"rowData": np.array(table_data), "rowCount": results["total"]}


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
