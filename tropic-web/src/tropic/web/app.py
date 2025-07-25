import dash
import dash_mantine_components as dmc
from dash import Dash, _dash_renderer, html

from tropic.web.components.footer import FOOTER
from tropic.web.components.header import HEADER

_dash_renderer._set_react_version("18.2.0")

app = Dash(
    external_stylesheets=dmc.styles.ALL,
    use_pages=True,
    external_scripts=[
        "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.4.0/3Dmol-min.js",
        "https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js",
    ],
)
app.layout = dmc.MantineProvider(
    html.Div(
        [
            html.Div(
                dmc.Container([HEADER, dash.page_container], size="xl"),
                style={"flex": "1"},
            ),
            FOOTER,
        ],
        style={
            "display": "flex",
            "flexDirection": "column",
            "minHeight": "100vh",
        },
    )
)


def main():
    app.run(debug=True)
