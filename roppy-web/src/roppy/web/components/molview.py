import dash_mantine_components as dmc
from dash import (
    ClientsideFunction,
    Input,
    Output,
    State,
    callback,
    clientside_callback,
    html,
)


def get_mol_viewer_card(xyz):
    return dmc.Card(
        [
            html.Div(
                id="viewer-container",
                style={
                    "width": "100%",
                    "aspect-ratio": "1 / 1",
                    "position": "relative",
                },
            ),
            html.Div(id="molecule-data", style={"display": "none"}, children=xyz),
            html.Div(id="page-load-trigger"),
            html.Div(id="javascript-trigger"),
        ],
        withBorder=True,
        shadow="sm",
        radius="md",
        py=0,
        px=0,
    )


def register_molview_callbacks():
    @callback(
        Output("molecule-data", "children"),
        Input("page-load-trigger", "children"),
        State("molecule-data", "children"),
    )
    def update_molecule_data(_, data):
        # Sample XYZ data for water molecule
        return data

    clientside_callback(
        ClientsideFunction(namespace="molecule_viewer", function_name="setupViewer"),
        Output("javascript-trigger", "children"),
        Input("molecule-data", "children"),
    )
