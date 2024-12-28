import dash_mantine_components as dmc

LINKS = {
    "Home": "/",
    "Monomer Search": "/monomers",
    "About": "/about",
    "API": "/api",
}
BC_STYLE = {"root": {"fontWeight": 700}}


def get_breadcrumbs(pages: list[str]):
    links = []
    for i, page in enumerate(pages):
        if page in LINKS and i < len(pages) - 1:
            links.append(dmc.Anchor(page, href=LINKS[page], styles=BC_STYLE))
        else:
            links.append(dmc.Text(page, styles=BC_STYLE))

    return dmc.Box(
        [dmc.Breadcrumbs(links)],
        py=10,
        px=20,
        mt=10,
        style={"backgroundColor": dmc.DEFAULT_THEME["colors"]["gray"][1]},
    )