import dash_mantine_components as dmc
from dash import dcc
from dash_iconify import DashIconify

HEADER = dmc.Box(
    children=[
        dmc.Container(
            style={"height": 60},
            pt="md",
            size="xl",
            children=[
                dmc.Group(
                    justify="space-between",
                    style={"height": "100%"},
                    px=0,
                    children=[
                        dmc.Group(
                            gap="md",
                            children=[
                                DashIconify(
                                    icon="eos-icons:molecules-outlined",
                                    color=dmc.DEFAULT_THEME["colors"]["blue"][5],
                                    width=30,
                                ),
                                dmc.Text(
                                    "ROPdb",
                                    size="xl",
                                    styles={"root": {"fontWeight": 700}},
                                ),
                            ],
                        ),
                        dmc.Group(
                            gap="sm",
                            children=[
                                dcc.Link(
                                    dmc.Button(
                                        "Home",
                                        variant="subtle",
                                        color="black",
                                        size="sm",
                                        styles={"root": {"fontWeight": 500}},
                                    ),
                                    href="/",
                                ),
                                dcc.Link(
                                    dmc.Button(
                                        "Search",
                                        variant="subtle",
                                        color="black",
                                        size="sm",
                                        styles={"root": {"fontWeight": 500}},
                                    ),
                                    href="/monomers",
                                ),
                                dcc.Link(
                                    dmc.Button(
                                        "About",
                                        variant="subtle",
                                        color="black",
                                        size="sm",
                                        styles={"root": {"fontWeight": 500}},
                                    ),
                                    href="/about",
                                ),
                                dcc.Link(
                                    dmc.Button(
                                        "API",
                                        variant="subtle",
                                        color="black",
                                        size="sm",
                                        styles={"root": {"fontWeight": 500}},
                                    ),
                                    href="/api",
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
        dmc.Divider(),
    ],
)
