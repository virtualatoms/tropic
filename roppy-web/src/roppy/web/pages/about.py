import dash
from dash import html
from roppy.web.header import get_breadcrumbs

dash.register_page(__name__)

breadcrumbs = get_breadcrumbs(["Home", "About"])

layout = html.Div(
    [
        breadcrumbs,
        html.H1(f"About"),
    ]
)
