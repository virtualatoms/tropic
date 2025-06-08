import dash_mantine_components as dmc
from dash import html


def get_search_chart():
    chart = dmc.ScatterChart(
        h=500,
        data=[
            {"color": "red.5", "name": "Computational", "data": []},
            {"color": "blue.5", "name": "Experimental", "data": []},
        ],
        dataKey={"x": "ring_size", "y": "delta_h"},
        xAxisLabel="Ring Size",
        yAxisLabel=r"ΔH / kj/mol",
        pointLabels="xy",
        tickLine="xy",
        mt=20,
        withLegend=True,
        id="analysis-chart",
    )
    y_select = dmc.Group(
        [
            dmc.Text("Y axis property:", size="sm", mr=-5),
            dmc.Select(
                # label="Y-axis property:",
                placeholder="Select one",
                id="chart-select-y",
                value="delta_h",
                data=[
                    {"value": "delta_h", "label": "ΔH"},
                    {"value": "delta_s", "label": "ΔS"},
                    {"value": "delta_g", "label": "ΔG"},
                    {"value": "ceiling_temperature", "label": "Ceiling Temperature"},
                ],
                w=200,
            ),
        ],
        justify="flex-end",
        mb=10,
    )

    return html.Div([y_select, chart])
