import dash_bio as dashbio
import dash_mantine_components as dmc


def get_draw_molecule():
    return dmc.Modal(
        title="Draw Molecule",
        children=[
            dashbio.Jsme(id="jsme", smiles=""),
            dmc.Space(h=20),
            dmc.Group(
                [
                    dmc.Button("Apply", id="apply-drawn-molecule"),
                    dmc.Button(
                        "Cancel",
                        id="close-drawer",
                        color="red",
                        variant="outline",
                    ),
                ],
                justify="flex-end",
            ),
        ],
        id="molecule-drawer-modal",
        size=650,
    )
