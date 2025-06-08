import dash_mantine_components as dmc
from dash_iconify import DashIconify

FOOTER = dmc.Container(
    [
        dmc.Container(
            [
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
                                dmc.Anchor("Contact", href="#", c="dimmed", size="sm"),
                                dmc.Anchor("Privacy", href="#", c="dimmed", size="sm"),
                                dmc.Anchor("Blog", href="#", c="dimmed", size="sm"),
                                dmc.Anchor("Careers", href="#", c="dimmed", size="sm"),
                            ],
                            justify="flex-end",
                        ),
                    ],
                    justify="space-between",
                    pt=30,
                    pb=30,
                    mt=30,
                ),
            ],
            size="xl",
            px=10,
        )
    ],
    size="100%",
    style={"backgroundColor": dmc.DEFAULT_THEME["colors"]["gray"][0]},
    className="footer",
)
