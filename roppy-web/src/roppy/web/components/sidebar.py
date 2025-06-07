import dash_mantine_components as dmc

radio_data = [["yes", "Yes"], ["no", "No"], ["both", "Both"]]


def get_search_sidebar():
    return dmc.Card(
        children=[
            dmc.CardSection(
                dmc.Text("Filters", fw=700),
                withBorder=True,
                inheritPadding=True,
                py="xs",
            ),
            dmc.Space(h=10),
            dmc.Text("Ring Size", size="sm", style={"font-weight": "500"}),
            dmc.RangeSlider(
                id="ring-size-slider",
                value=[1, 15],
                marks=[
                    {"value": 1, "label": "1"},
                    {"value": 5, "label": "5"},
                    {"value": 10, "label": "10"},
                    {"value": 15, "label": "15+"},
                ],
                min=1,
                max=15,
                minRange=1,
                mb=25,
                mt=10,
            ),
            dmc.Space(h=10),
            dmc.Text(
                "Molecular Weight (g/mol)", size="sm", style={"font-weight": "500"}
            ),
            dmc.RangeSlider(
                id="mol-weight-slider",
                value=[10, 500],
                marks=[
                    {"value": i, "label": str(i)}
                    for i in [10, 100, 200, 300, 400, 500]
                ],
                min=10,
                max=500,
                minRange=10,
                mb=25,
                mt=10,
            ),
            dmc.Space(h=10),
            dmc.MultiSelect(
                label="Functional Groups",
                placeholder="Select functional groups",
                id="functional-groups",
                value=[],
                data=[
                    {"value": "ester", "label": "Ester"},
                    {"value": "aldehyde", "label": "Aldehyde"},
                    {"value": "ketone", "label": "Ketone"},
                    {"value": "amide", "label": "Amide"},
                    {"value": "ether", "label": "Ether"},
                    {"value": "amine", "label": "Amine"},
                ],
                # mb=10,
                clearable=True,
            ),
            dmc.Space(h=10),
            dmc.RadioGroup(
                children=dmc.Group(
                    [dmc.Radio(l, value=k) for k, l in radio_data], my=10
                ),
                id="has-exp",
                value="both",
                label="Has Experimental Data",
                size="sm",
            ),
            dmc.Space(h=10),
            dmc.RadioGroup(
                children=dmc.Group(
                    [dmc.Radio(l, value=k) for k, l in radio_data], my=10
                ),
                id="has-comp",
                value="both",
                label="Has Computational Data",
                size="sm",
            ),
        ],
        withBorder=True,
        shadow="sm",
        radius="md",
    )
