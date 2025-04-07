from dash import Dash, html, dcc, dash_table, Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc
from beanie import init_beanie
import motor.motor_asyncio
from documents import MonomerSummary
import json
import base64
from rdkit import Chem
from rdkit.Chem import Draw
import io
import numpy as np
import dash_bio as dashbio


app = Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.LUX,
        "https://cdn.jsdelivr.net/gh/AnnMarieW/dash-bootstrap-templates/dbc.min.css"
    ],
)


def smiles_to_image(smiles, size=(150, 150)):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=size)
            buffer = io.BytesIO()
            img.save(buffer, format="PNG")
            encoded_image = base64.b64encode(buffer.getvalue()).decode()
            return encoded_image
    except Exception as e:
        print(f"Error converting SMILES to image: {e}")
        return None
    return None


app.layout = html.Div(
    [
        dbc.Modal(
            [
                dbc.ModalHeader("Draw Molecule"),
                html.Div(
                    [
                        dashbio.Jsme(
                            id="jsme",
                            smiles="",
                        ),
                    ],
                    className="mx-auto py-4",
                ),
                dbc.ModalFooter(
                    [
                        dbc.Button("Cancel", id="close-drawer", className="ml-auto"),
                        dbc.Button("Apply", id="apply-drawn-molecule", color="primary"),
                    ]
                ),
            ],
            id="molecule-drawer-modal",
            size="lg",
        ),
        html.Center([html.H1("RopPy", className="mt-4")]),
        html.Div(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Input(
                                    id="smiles-input",
                                    type="text",
                                    placeholder="Enter SMILES notation",
                                    className="mb-3",
                                ),
                            ],
                            width=5,
                        ),
                        dbc.Col(
                            [
                                dbc.Button(
                                    "Draw",
                                    id="draw-button",
                                    color="primary",
                                    className="mb-3",
                                ),
                                html.Div(
                                    id="molecule-drawer-output",
                                    style={"display": "none"},
                                ),
                            ],
                            width="auto",
                        ),
                    ],
                    justify="center",
                    className="mt-5",
                ),
            ]
        ),
        dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader("Filters"),
                                        dbc.CardBody(
                                            [
                                                html.Label("Ring Size"),
                                                dcc.RangeSlider(
                                                    id="ring-size-slider",
                                                    min=0,
                                                    max=5,
                                                    step=1,
                                                    marks={
                                                        0: "0",
                                                        1: "1",
                                                        2: "2",
                                                        3: "3",
                                                        4: "4",
                                                        5: "many",
                                                    },
                                                    value=[0, 5],
                                                ),
                                            ]
                                        ),
                                    ]
                                )
                            ],
                            width=3,
                        ),
                        dbc.Col(
                            [
                                dash_table.DataTable(
                                    id="results-table",
                                    columns=[
                                        {
                                            "name": "Structure",
                                            "id": "structure",
                                            "presentation": "markdown",
                                        },
                                        {"name": "Monomer ID", "id": "monomer_id"},
                                        {"name": "SMILES", "id": "smiles"},
                                        {
                                            "name": "Ring Size",
                                            "id": "ring_size",
                                            "type": "numeric",
                                        },
                                        {
                                            "name": "Avg ΔH",
                                            "id": "average_delta_h",
                                            "type": "numeric",
                                        },
                                        {
                                            "name": "Avg ΔS",
                                            "id": "average_delta_s",
                                            "type": "numeric",
                                        },
                                    ],
                                    cell_selectable=False,
                                    sort_action="native",  # Enable native sorting
                                    sort_mode="single",  # Allow sorting by one column at a time
                                    markdown_options={"html": True},
                                    style_table={"overflowX": "auto"},
                                    style_cell={"textAlign": "left", "padding": "10px"},
                                    style_cell_conditional=[
                                        {
                                            "if": {"column_id": "structure"},
                                            "width": "150px",
                                            "padding": "5px",
                                        },
                                        {
                                            "if": {"column_id": "average_delta_h"},
                                            "textAlign": "right",
                                        },
                                        {
                                            "if": {"column_id": "average_delta_s"},
                                            "textAlign": "right",
                                        },
                                        {
                                            "if": {"column_id": "ring_size"},
                                            "textAlign": "right",
                                        },
                                    ],
                                    style_header={
                                        "backgroundColor": "rgb(230, 230, 230)",
                                        "fontWeight": "bold",
                                        # 'cursor': 'pointer'  # Add pointer cursor for sortable headers
                                    },
                                )
                            ]
                        ),
                    ],
                    className="mt-4 dbc",
                )
            ]
        ),
    ]
)


async def update_table_async(smiles, ring_size_range):
    # Initialize database connection
    client = motor.motor_asyncio.AsyncIOMotorClient("mongodb://localhost:27017")
    await init_beanie(database=client.roppy, document_models=[MonomerSummary])

    query = {}
    if smiles:
        query["smiles"] = {"$regex": smiles, "$options": "i"}

    if ring_size_range[1] == 5:
        # include many
        ring_query = {
            "$or": [
                {"ring_size": {"$gte": ring_size_range[0], "$lte": ring_size_range[1]}},
                {"ring_size": "many"},
            ]
        }
    else:
        ring_query = {
            "$or": [
                {"ring_size": {"$gte": ring_size_range[0], "$lte": ring_size_range[1]}}
            ]
        }

    query.update(ring_query)

    results = await MonomerSummary.find(query).to_list()

    table_data = []
    for result in results:
        # Convert SMILES to image
        img_data = smiles_to_image(result.smiles)
        structure_cell = (
            f'<img src="data:image/png;base64,{img_data}" style="max-width:150px;">'
            if img_data
            else ""
        )

        table_data.append(
            {
                "structure": structure_cell,
                "monomer_id": result.monomer_id,
                "smiles": result.smiles,
                "ring_size": result.ring_size,
                "average_delta_h": round(result.average_delta_h, 2),
                "average_delta_s": round(result.average_delta_s, 2),
            }
        )

    return np.array(table_data)


@callback(
    Output("results-table", "data"),
    Input("smiles-input", "value"),
    Input("ring-size-slider", "value"),
)
def update_table(smiles, ring_size_range):
    import asyncio

    return asyncio.run(update_table_async(smiles, ring_size_range))


@callback(
    Output("molecule-drawer-modal", "is_open"),
    [
        Input("draw-button", "n_clicks"),
        Input("close-drawer", "n_clicks"),
        Input("apply-drawn-molecule", "n_clicks"),
    ],
    [State("molecule-drawer-modal", "is_open")],
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


if __name__ == "__main__":
    app.run_server(debug=True)
