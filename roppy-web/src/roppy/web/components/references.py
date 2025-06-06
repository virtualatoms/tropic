import dash_mantine_components as dmc


def get_references():
    rows = [dmc.TableTr(["", ""])]
    head = dmc.TableThead(dmc.TableTr([dmc.TableTh("Reference"), dmc.TableTh("DOI")]))
    body = dmc.TableTbody(rows)
    return dmc.Table([head, body], id="references-table")


def get_reference_table_data(dois, formatted_references):
    """
    Creates a table of references with formatted strings and DOIs.
    """
    rows = [
        dmc.TableTr(
            [
                dmc.TableTd(ref),
                dmc.TableTd(dmc.Anchor(doi, href=f"https://doi.org/{doi}")),
            ]
        )
        for ref, doi in zip(formatted_references, dois)
    ]

    return [
        dmc.TableThead(dmc.TableTr([dmc.TableTh("Reference"), dmc.TableTh("DOI")])),
        dmc.TableTbody(rows),
    ]
