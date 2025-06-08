import dash_mantine_components as dmc
from dash import html


def get_references():
    rows = [dmc.TableTr(["", ""])]
    head = dmc.TableThead(dmc.TableTr([dmc.TableTh("Reference"), dmc.TableTh("DOI")]))
    body = dmc.TableTbody(rows)
    return dmc.Table([head, body], id="references-table", mt=20)


def get_reference_table_data(dois, formatted_references):
    """
    Creates a table of references with formatted strings and DOIs.
    """
    rows = [
        dmc.TableTr(
            [
                dmc.TableTd(ref),
                dmc.TableTd(
                    html.A(
                        doi,
                        href=f"https://doi.org/{doi}",
                        className="mantine-focus-auto m_849cf0da mantine-Text-root mantine-Anchor-root",
                        style={"text-decoration": "underline"},
                    )
                ),
            ]
        )
        for ref, doi in zip(formatted_references, dois)
    ]

    return [
        dmc.TableThead(dmc.TableTr([dmc.TableTh("Reference"), dmc.TableTh("DOI")])),
        dmc.TableTbody(rows),
    ]
