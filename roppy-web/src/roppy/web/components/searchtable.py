import dash_ag_grid as dag

PAGE_SIZE = 5

def get_search_table():
    return dag.AgGrid(
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
            "paginationPageSize": PAGE_SIZE,
            "rowHeight": 100,
        },
        style={"height": None},
        dangerously_allow_code=True,
    )