import dash
from dash import html
import dash_mantine_components as dmc
from roppy.web.utils import smiles_to_image
from roppy.web.components.breadcrumbs import get_breadcrumbs
from roppy.web.components.toc import create_toc_sidebar, register_toc_callbacks


dash.register_page(__name__, path_template="/monomers/<monomer_id>")

register_toc_callbacks()

lorem = """Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean nec ullamcorper neque. Quisque non leo elementum, porttitor leo at, aliquam nisi. Proin elementum eleifend libero, vel aliquam diam sollicitudin eget. Nam elementum eget odio eget bibendum. Morbi ac tellus nibh. Nullam faucibus non ipsum ac blandit. In et sem justo.

Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Nunc sollicitudin, sapien non vehicula blandit, lectus tellus blandit libero, vel luctus augue libero et tortor. Fusce iaculis consectetur mi et aliquet. Vivamus lectus magna, finibus in placerat vel, blandit maximus risus. Nulla feugiat lobortis elementum. Mauris eu velit enim. Suspendisse eget efficitur ex. Etiam nec euismod sapien. Curabitur id neque nisl. Etiam lobortis in justo ut mollis. Integer tempus molestie imperdiet. Phasellus vel accumsan elit. Sed non tincidunt lectus. Suspendisse faucibus, lorem sit amet facilisis varius, felis lorem rutrum nulla, in tempor magna eros sit amet nunc.

Aliquam erat volutpat. Cras nec risus dolor. Sed dignissim libero vitae nisl condimentum posuere. Sed quis justo ut leo pretium maximus eu ac diam. Fusce id mauris et justo congue volutpat. Suspendisse tincidunt, enim in lobortis convallis, neque felis tincidunt velit, nec vulputate est felis sit amet magna. Etiam felis ipsum, laoreet in sollicitudin nec, vulputate a tortor. Phasellus non condimentum quam. Quisque sed porttitor ex.

Praesent consequat diam feugiat felis fringilla faucibus eu sed ligula. Aliquam venenatis velit quam, at consequat risus facilisis eu. Sed non vestibulum mi, a vestibulum leo. Suspendisse in ultricies ligula, sed eleifend massa. Nunc et nisl sem. In hac habitasse platea dictumst. Ut tincidunt gravida mauris, vel porta tellus porttitor vitae. Curabitur eros magna, volutpat nec fringilla aliquam, tristique at augue. Proin cursus arcu ex, eu accumsan tellus egestas ac. Fusce luctus velit ullamcorper iaculis ornare. Fusce tristique, est ut scelerisque pharetra, ante libero convallis massa, eget pretium risus mauris a velit. Sed vehicula mollis libero, vitae lacinia tellus suscipit in. In hac habitasse platea dictumst. Proin venenatis faucibus leo, non congue ex pretium ut.

Phasellus ac justo lorem. Curabitur suscipit ligula sem, eu consequat mauris convallis ac. Proin metus sem, pretium in blandit non, egestas vitae tortor. Cras eu mollis dolor. Vestibulum sit amet porta turpis, ut tempus ante. Suspendisse blandit justo ante. In molestie diam vel ex molestie congue. Aenean pretium, ex aliquam blandit aliquet, felis mi aliquet est, quis feugiat eros mi eu augue. Duis id arcu justo. Aliquam erat volutpat. Donec vitae porta est. Aliquam hendrerit sem vitae risus aliquam vehicula. Nam ac risus aliquet, tempus nibh sed, blandit felis."""

def layout(monomer_id=None, **kwargs):
    breadcrumbs = get_breadcrumbs(["Home", "Monomer Search", f"{monomer_id}"])

    img_data = smiles_to_image("COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C", size=(200, 200))
    img = dmc.Image(
        radius="md",
        h=200,
        src=f"data:image/svg+xml;base64,{img_data}",
    )

    monomer_card1 = dmc.Card(
        [img],
        withBorder=True,
        shadow="sm",
        radius="md",
    )
    monomer_card2 = dmc.Card(
        [img],
        withBorder=True,
        shadow="sm",
        radius="md",
    )
    monomer_card3 = dmc.Card(
        [img],
        withBorder=True,
        shadow="sm",
        radius="md",
    )

    # dmc.Card(
    #     children=[
    #         dmc.CardTitle("Monomer"),
    #         dmc.CardContent(
    #             children=[
    #                 dmc.Text("Monomer ID:"),
    #                 dmc.Text("SMILES:"),
    #                 dmc.Text("InChI:"),
    #                 dmc.Text("InChI Key:"),
    #                 dmc.Text("IUPAC Name:"),
    #                 dmc.Text("Common Name:"),
    #             ]
    #         ),
    #     ]
    # )
    # 
    content = html.Div([
        html.H1("Main Title", id="main-title"),
        dmc.Text(lorem),
        html.H1("Section 1", id="section-1"),
        dmc.Text(lorem),
        html.H1("Section 2", id="section-2"),
        dmc.Text(lorem),
    ], className="main-content", id="page-content")

    grid = dmc.Grid(
        [
            dmc.GridCol(create_toc_sidebar(), span=3),
            dmc.GridCol(content, span=3),
            dmc.GridCol(monomer_card3, span=3),
        ],
        pt=50,
        gutter="xl",
    )
    return [     
        # create_toc_sidebar(),
        breadcrumbs, 
        grid
    ]
