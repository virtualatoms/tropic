import dash
from dash import html
from roppy.web.components.breadcrumbs import get_breadcrumbs

dash.register_page(__name__)

breadcrumbs = get_breadcrumbs(["Home", "API"])

layout = html.Div(
    [
        breadcrumbs,
        html.H1(f"API"),
    ]
)
