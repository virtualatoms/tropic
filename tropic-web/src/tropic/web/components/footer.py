import dash_mantine_components as dmc
from dash import html
from dash_iconify import DashIconify

_LINK_STYLE = {
    "color": dmc.DEFAULT_THEME["colors"]["dark"][2],
    "font-size": dmc.DEFAULT_THEME["fontSizes"]["sm"],
    "textDecoration": "none",
}
FOOTER = html.Div(
    dmc.Container(
        dmc.Group(
            [
                dmc.Group(
                    gap="sm",
                    children=[
                        DashIconify(
                            icon="eos-icons:molecules-outlined",
                            color=dmc.DEFAULT_THEME["colors"]["dark"][1],
                            width=24,
                        ),
                        dmc.Text(
                            "TROPIC",
                            size="sm",
                            c="dimmed",
                            styles={"root": {"fontWeight": 600}},
                        ),
                    ],
                ),
                dmc.Group(
                    [
                        html.A("About", href="/about", style=_LINK_STYLE),
                        html.A("API", href="/api", style=_LINK_STYLE),
                        html.A(
                            "GitHub",
                            href="https://github.com/virtualatoms/tropic",
                            style=_LINK_STYLE,
                        ),
                    ],
                    justify="flex-end",
                ),
            ],
            justify="space-between",
            py=20,
            px=15,
        ),
        size="xl",
    ),
    style={"backgroundColor": dmc.DEFAULT_THEME["colors"]["gray"][0]},
    className="footer",
)
