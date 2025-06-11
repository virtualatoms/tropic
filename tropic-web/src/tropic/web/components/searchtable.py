import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash import dcc, html

SEARCH_NUM_ROWS = 5


def get_search_table():
    table = dag.AgGrid(
        id="results-table",
        defaultColDef={
            "resizable": True,
            "cellStyle": {"align-items": "center", "height": "100%", "display": "flex"},
        },
        columnDefs=[
            {
                "headerName": "Molecule",
                "field": "structure",
                "cellRenderer": "markdown",
            },
            {
                "headerName": "Monomer ID",
                "field": "monomer_id",
                "cellRenderer": "markdown",
            },
            {"headerName": "SMILES", "field": "smiles"},
            {"headerName": "Ring Size", "field": "ring_size"},
            {"headerName": "Has Exp", "field": "has_exp"},
            {"headerName": "Has Comp", "field": "has_comp", "resizable": False},
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
            "cacheBlockSize": SEARCH_NUM_ROWS,
            "domLayout": "autoHeight",
            "suppressCellFocus": True,
            "pagination": True,
            "paginationPageSizeSelector": False,
            "paginationPageSize": SEARCH_NUM_ROWS,
            "rowHeight": 100,
        },
        style={"height": None},
        dangerously_allow_code=True,
    )
    export = dmc.Group(
        [
            dmc.Text("Export as:", size="sm", mr=-5),
            dmc.Select(
                # label="Export format",
                placeholder="Select one",
                id="export-data-select",
                value="csv",
                data=[
                    {"value": "csv", "label": "CSV"},
                    {"value": "json", "label": "JSON"},
                ],
                w=120,
            ),
            dmc.Button(
                "Export",
                id="export-data-button",
                variant="outline",
                color="blue",
                size="sm",
                styles={"root": {"fontWeight": 500}},
            ),
        ],
        justify="flex-end",
        mb=10,
    )

    return html.Div([export, table, dcc.Download(id="download-data")])
