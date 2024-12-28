from dash import register_page, html

register_page(__name__, path="/")


layout = html.Div(
    [
        html.H1(f"Welcome"),
    ]
)
