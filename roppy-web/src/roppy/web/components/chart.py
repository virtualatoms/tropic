import dash_mantine_components as dmc

def get_search_chart():
    data = [
        {
            "color": "blue.5",
            "name": "Experimental",
            "data": [
                {"ring_size": x, "delta_h": y, "delta_s": z}
                for x, y, z in zip(range(1, 11), range(10, 20), range(20, 30))
            ],
        },
        {
            "color": "red.5",
            "name": "Computational",
            "data": [
                {"ring_size": x, "delta_h": y, "delta_s": z}
                for x, y, z in zip(range(1, 11), range(20, 30), range(30, 40))
            ],
        },
    ]

    chart = dmc.ScatterChart(
        h=500,
        data=data,
        dataKey={"x": "ring_size", "y": "delta_h"},
        xAxisLabel="Ring Size",
        yAxisLabel=r"ΔH",
        pointLabels="xy",
        tickLine="xy",
        mt=20,
    )
    return chart