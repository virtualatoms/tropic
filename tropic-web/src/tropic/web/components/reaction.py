import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash import html

from tropic.web.utils import reaction_to_image

LINK_SVG = '<svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24"><path fill="none" stroke="currentColor" stroke-linecap="round" stroke-linejoin="round" stroke-width="3" d="M13.5 10.5L21 3m-5 0h5v5m0 6v5a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5"/></svg>'

# custom header for ag-grid to allow html
HTML_HEADER = """
<div class="ag-cell-label-container" role="presentation">
    <span ref="eMenu" class="ag-header-icon ag-header-cell-menu-button"></span>
    <span ref="eFilterButton" class="ag-header-icon ag-header-cell-filter-button"></span>
    <div ref="eLabel" class="ag-header-cell-label" role="presentation">
        <span ref="eSortOrder" class="ag-header-icon ag-sort-order ag-hidden"></span>
        <span ref="eSortAsc" class="ag-header-icon ag-sort-ascending-icon ag-hidden"></span>
        <span ref="eSortDesc" class="ag-header-icon ag-sort-descending-icon ag-hidden"></span>
        <span ref="eSortMixed" class="ag-header-icon ag-sort-mixed-icon ag-hidden"></span>
        <span ref="eSortNone" class="ag-header-icon ag-sort-none-icon ag-hidden"></span>
        <span>{}</span>
        <span ref="eFilter" class="ag-header-icon ag-filter-icon"></span>
    </div>
</div>
"""

COMMON_COLUMNS = [
    {
        "headerName": "Reaction ID",
        "field": "reaction_id",
        "cellRenderer": "markdown",
        "width": 140,
        "maxWidth": 140,
    },
    {
        "headerName": "Type",
        "field": "type",
        "width": 100,
        "maxWidth": 100,
    },
    {
        "headerName": "",
        "field": "delta_h",
        "headerComponentParams": {
            "template": HTML_HEADER.format("ΔH<sub>p</sub> (kJ/mol)")
        },
    },
    {
        "headerName": "",
        "field": "delta_s",
        "headerComponentParams": {
            "template": HTML_HEADER.format("ΔS<sub>p</sub> (J/K/mol)")
        },
    },
    {
        "headerName": "",
        "field": "ceiling_temperature",
        "headerComponentParams": {"template": HTML_HEADER.format("T<sub>c</sub> (K)")},
    },
    {"headerName": "Year", "field": "year", "width": 120, "maxWidth": 120},
    {
        "headerName": "Ref.",
        "field": "ref",
        "cellRenderer": "markdown",
        "width": 80,
        "maxWidth": 80,
        "resizable": False,
    },
]
EXP_COLUMNS = (
    COMMON_COLUMNS[:2]
    + [
        {"headerName": "State", "field": "state", "width": 90, "maxWidth": 90},
        # {
        #     "headerName": "Solvent",
        #     "field": "solvent",
        # },
        # {
        #     "headerName": "Initial Mon. Conc.",
        #     "field": "initial_monomer_conc",
        # },
        # {
        #     "headerName": "Bulk Mon. Conc.",
        #     "field": "bulk_monomer_conc",
        # }
    ]
    + COMMON_COLUMNS[2:]
)
COMP_COLUMNS = (
    COMMON_COLUMNS[:2]
    + [
        {
            "headerName": "Method",
            "field": "method",
            "width": 120,
            "maxWidth": 120,
        },
        {
            "headerName": "# Units",
            "field": "repeating_units",
            "width": 110,
            "maxWidth": 110,
        },
    ]
    + COMMON_COLUMNS[2:]
)


def get_reaction_reaction_card(data):
    reaction_data = reaction_to_image(data["reaction_smiles"], size=(-1, -1))
    reaction_img = dmc.Image(h=200, src=f"data:image/svg+xml;base64,{reaction_data}")
    return dmc.Card(reaction_img, withBorder=True, shadow="sm", radius="md")


def get_reaction_section(data):
    exp_table_data = []
    comp_table_data = []
    for row in data["data"]:
        if row["ceiling_temperature"]:
            ceiling_temperature = f"{row['ceiling_temperature']:.2f} K"
        elif (
            row["delta_h"]
            and row["delta_s"]
            and row["delta_s"] < 0
            and row["delta_h"] < 0
        ):
            ceiling_temperature = f"{row['delta_h'] * 1000 / row['delta_s']:.2f} K"
        else:
            ceiling_temperature = ""

        common_data = {
            "reaction_id": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="reactions/{row["reaction_id"]}">{row["reaction_id"]}</a>',
            "type": row["type"],
            "delta_h": f"{row['delta_h']:.2f}" if row["delta_h"] else "",
            "delta_s": f"{row['delta_s']:.2f}" if row["delta_s"] else "",
            "ceiling_temperature": ceiling_temperature,
            "year": row["year"],
            "ref": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="https://doi.org/{row["doi"]}" target="_blank">{LINK_SVG}</a>',
        }
        if row["is_experimental"]:
            common_data |= {
                "state": row["state_summary"],
                # "solvent": "-",
                "initial_monomer_conc": f"{row['initial_monomer_conc']:.2f}"
                if row["initial_monomer_conc"]
                else "",
                "bulk_monomer_conc": f"{row['bulk_monomer_conc']:.2f}"
                if row["bulk_monomer_conc"]
                else "",
            }
            exp_table_data.append(common_data)
        else:
            common_data |= {
                "method": row["method"],
                "repeating_units": row["repeating_units"],
            }
            comp_table_data.append(common_data)

    reaction = []
    if exp_table_data:
        exp_table = get_table(
            exp_table_data, EXP_COLUMNS, section_id="reaction-exp-table"
        )
        exp_section = dmc.Stack(
            [
                html.H1(
                    "Experimental Data",
                    id="reaction",
                ),
                exp_table,
            ],
            gap=0,
        )
        reaction.append(exp_section)

    if comp_table_data:
        comp_table = get_table(
            comp_table_data, COMP_COLUMNS, section_id="reaction-comp-table"
        )
        comp_section = dmc.Stack(
            [
                html.H1("Computational Data", id="reaction-comp"),
                comp_table,
            ],
            gap=0,
        )
        reaction.append(comp_section)
    return dmc.Stack(reaction)


def get_table(data, columns, section_id=None):
    table = dag.AgGrid(
        id=section_id,
        columnDefs=columns,
        columnSize="responsiveSizeToFit",
        columnSizeOptions={
            "defaultMinWidth": 10,
        },
        rowData=data,
        dashGridOptions={
            "domLayout": "autoHeight",
            "suppressCellFocus": True,
        },
        style={"height": None},
        dangerously_allow_code=True,
        className="ag-theme-alpine compact",
    )

    return table


# def get_reaction_section(data):
#     poly_info_card = get_reaction_info_card(data["reaction"])
#     poly_reaction_card = get_reaction_reaction_card(data)

#     summary = dmc.Grid(
#         [
#             dmc.GridCol(reaction_reaction_card, span=8),
#             dmc.GridCol(reaction_info_card, span=4),
#         ],
#         gutter="xl",
#     )

#     reaction = dmc.Stack(
#         [
#             html.H1(
#                 "Reaction", id="reaction", style={"margin-top": 0, "padding-top": 0}
#             ),
#             summary,
#             html.H2(
#                 "Experimental",
#                 id="reaction-exp",
#                 style={"margin-bottom": 10, "padding-bottom": 0},
#             ),
#             html.H2("Computational", id="reaction-comp"),
#         ],
#         gap=20,
#     )
#     return reaction

# def get_poly_info_card(polyinfo):
#     table_data = [
#         (
#             "PolyInfo",
#             dcc.Link(
#                 polyinfo["polyinfo_id"],
#                 href=f"https://polyinfo.com/{polyinfo['polyinfo_id']}",
#             ),
#         ),
#         (
#             "PolyGenome",
#             dcc.Link(
#                 polyinfo["polygenome_id"],
#                 href=f"https://polygenome.com/{polyinfo['polygenome_id']}",
#             ),
#         ),
#     ]
#     return dmc.Card([get_table(table_data)], withBorder=True, shadow="sm", radius="md")
