import dash_mantine_components as dmc
import dash
from dash import Dash, _dash_renderer
from roppy.web.components.header import HEADER

_dash_renderer._set_react_version("18.2.0")
app = Dash(external_stylesheets=dmc.styles.ALL, use_pages=True)

app.layout = dmc.MantineProvider([dmc.Container([HEADER, dash.page_container], size="xl")])

def main():
    app.run_server(debug=True)
