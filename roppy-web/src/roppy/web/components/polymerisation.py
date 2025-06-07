import dash_mantine_components as dmc
import dash_ag_grid as dag
from dash import dcc, html
from datetime import datetime

from roppy.web.components.table import get_table
from roppy.web.utils import reaction_to_image

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
        "headerName": "Polymer ID",
        "field": "polymerisation_id",
        "cellRenderer": "markdown",
    },
    {"headerName": "Type", "field": "type"},
    {
        "headerName": "", 
        "field": "delta_h",  
        "headerComponentParams": {
            "template": HTML_HEADER.format("ΔH<sub>p</sub> (kJ/mol)")
        }
    },
    {
        "headerName": "", 
        "field": "delta_s",  
        "headerComponentParams": {
            "template": HTML_HEADER.format("ΔS<sub>p</sub> (kJ/mol)")
        }
    },
    {
        "headerName": "", 
        "field": "ceiling_temperature",  
        "headerComponentParams": {
            "template": HTML_HEADER.format("T<sub>c</sub> (K)")
        }
    },
    {"headerName": "Year", "field": "year"},
    {
        "headerName": "Ref.",
        "field": "ref",
        "cellRenderer": "markdown",
    },
]
EXP_COLUMNS = COMMON_COLUMNS[:2] + [
    {"headerName": "State", "field": "state"},
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
] + COMMON_COLUMNS[2:]

def get_polymer_info_card(polyinfo):
    table_data = [
        (
            "PolyInfo",
            dcc.Link(
                polyinfo["polyinfo_id"],
                href=f"https://polyinfo.com/{polyinfo['polyinfo_id']}",
            ),
        ),
        (
            "PolyGenome",
            dcc.Link(
                polyinfo["polygenome_id"],
                href=f"https://polygenome.com/{polyinfo['polygenome_id']}",
            ),
        ),
    ]
    return dmc.Card([get_table(table_data)], withBorder=True, shadow="sm", radius="md")


def get_polymer_reaction_card(data):
    reaction_data = reaction_to_image(data["reaction_smiles"], size=(-1, -1))
    reaction_img = dmc.Image(h=200, src=f"data:image/svg+xml;base64,{reaction_data}")
    return dmc.Card(reaction_img, withBorder=True, shadow="sm", radius="md")

def get_polymerisation_section(data):
    exp_table_data = []
    comp_table_data = []
    for row in data["data"]:
        common_data = {
            "polymerisation_id": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="polymers/{row["polymerisation_id"]}">{row["polymerisation_id"]}</a>',
            "type": row["type"],
            "delta_h": f'{row["delta_h"]:.2f}' if row["delta_h"] else "",
            "delta_s": f'{row["delta_s"]:.2f}' if row["delta_s"] else "",
            "ceiling_temperature": f'{row["ceiling_temperature"]:.2f} K' if row["ceiling_temperature"] else "",
            "year": row["year"],
            "ref": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="https://doi.org/{row["doi"]}">{LINK_SVG}</a>',
            # "doi": f'<a class="mantine-focus-auto  m_849cf0da mantine-Text-root mantine-Anchor-root" data-underline="always" href="https://doi.org/{row["doi"]}">{row["doi"]}</a>',
        }
        if row["is_experimental"]:
            common_data |= {
                "state": row["state_summary"],
                # "solvent": "-",
                "initial_monomer_conc": f'{row["initial_monomer_conc"]:.2f}' if row["initial_monomer_conc"] else "",
                "bulk_monomer_conc": f'{row["bulk_monomer_conc"]:.2f}' if row["bulk_monomer_conc"] else "",
            }
            exp_table_data.append(common_data)
        else:
            comp_table_data.append(common_data)
    
    poly = []
    if exp_table_data:
        exp_table = get_table(exp_table_data, EXP_COLUMNS, id="poly-exp-table")
        exp_section = dmc.Stack(
            [
                html.H1(
                    "Experimental Data", id="poly", 
                ),
                exp_table,
            ],
            gap=0,
        )
        poly.append(exp_section)
    
    if comp_table_data:
        comp_table = get_table(comp_table_data, COMMON_COLUMNS, id="poly-comp-table")
        comp_section = dmc.Stack(
            [
                html.H1("Computational Data", id="poly-comp"),
                comp_table,
            ],
            gap=0,
        )
        poly.append(comp_section)
    return dmc.Stack(poly)

def get_table(data, columns, id=None):
    table = dag.AgGrid(
        id=id,
        columnDefs=columns,
        columnSize="responsiveSizeToFit",
        rowData=data,
        dashGridOptions={
            "domLayout": "autoHeight",
            "suppressCellFocus": True,
        },
        style = {"height": None},
        dangerously_allow_code=True,
    )

    return table



# def get_polymerisation_section(data):
#     poly_info_card = get_polymer_info_card(data["polymer"])
#     poly_reaction_card = get_polymer_reaction_card(data)

#     summary = dmc.Grid(
#         [
#             dmc.GridCol(poly_reaction_card, span=8),
#             dmc.GridCol(poly_info_card, span=4),
#         ],
#         gutter="xl",
#     )

#     poly = dmc.Stack(
#         [
#             html.H1(
#                 "Polymerisation", id="poly", style={"margin-top": 0, "padding-top": 0}
#             ),
#             summary,
#             html.H2(
#                 "Experimental",
#                 id="poly-exp",
#                 style={"margin-bottom": 10, "padding-bottom": 0},
#             ),
#             html.H2("Computational", id="poly-comp"),
#         ],
#         gap=20,
#     )
#     return poly
