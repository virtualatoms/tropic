import dash_mantine_components as dmc


def get_table(data, header=None):
    row_style = {} if header else {"font-weight": "bold"}
    rows = [
        dmc.TableTr([dmc.TableTd(name, style=row_style), dmc.TableTd(value)])
        for name, value in data
    ]
    body = dmc.TableTbody(rows)

    if header:
        head = dmc.TableThead(dmc.TableTr([dmc.TableTh(name) for name in header]))
        return dmc.Table([body, head])
    return dmc.Table([body])
