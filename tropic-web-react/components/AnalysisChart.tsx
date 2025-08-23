import { ScatterChart } from "@mantine/charts";
import { Select, Group, Text, LoadingOverlay, Box } from "@mantine/core";
import { useMemo, useState } from "react";
import { MonomerSummaryState } from "@/lib/types";

type ChartDatum = {
	ring_size: number;
	delta_h: number | null;
	delta_s: number | null;
	ceiling_temperature: number | null;
};

const yOptions = [
	{ value: "delta_h", label: "ΔH" },
	{ value: "delta_s", label: "ΔS" },
	{ value: "ceiling_temperature", label: "Ceiling Temperature" },
];

function yAxisLabel(key: keyof ChartDatum) {
	switch (key) {
		case "delta_h":
			return "ΔH / kJ/mol";
		case "delta_s":
			return "ΔS / J/mol·K";
		case "ceiling_temperature":
			return "Ceiling Temp / °C";
		default:
			return key;
	}
}


export default function AnalysisChart({ data, isLoading }: MonomerSummaryState) {
	const [yField, setYField] = useState<keyof ChartDatum>("delta_h");

	const chartData = useMemo(() => {
		const expData: ChartDatum[] = [];
		const compData: ChartDatum[] = [];

		for (const monomerSummary of data) {
			for (const d of monomerSummary.data) {
				if (d.delta_h === null || d.delta_s === null) continue;

				const datum: ChartDatum = {
					ring_size: monomerSummary.monomer.ring_size,
					delta_h: d.delta_h,
					delta_s: d.delta_s,
					ceiling_temperature: d.ceiling_temperature,
				};

				if (d.is_experimental) {
					expData.push(datum);
				} else {
					compData.push(datum);
				}
			}
		}
		return { experimental: expData, computational: compData };
	}, [data]);

	const formattedData = useMemo(() => {
		return [
			{
				name: "Comp",
				color: "red.5",
				data: chartData.computational.map((d) => ({
					x: d.ring_size,
					y: d[yField],
				})),
			},
			{
				name: "Exp",
				color: "blue.5",
				data: chartData.experimental.map((d) => ({
					x: d.ring_size,
					y: d[yField],
				})),
			},
		];
	}, [chartData, yField]);

	return (
		<Box pos="relative">
			<LoadingOverlay visible={isLoading} />
			<Group justify="flex-end" mb="md">
				<Text size="sm">Y axis property:</Text>
				<Select
					placeholder="Select one"
					data={yOptions}
					value={yField}
					onChange={(val) => setYField(val as keyof ChartDatum)}
					w={220}
				/>
			</Group>
			<ScatterChart
				h={500}
				data={formattedData}
				dataKey={{ x: "x", y: "y" }}
				xAxisLabel="Ring Size"
				yAxisLabel={yAxisLabel(yField)}
				withLegend
        		// pointLabels={undefined}
				// tickLine="xy"
        		// mt={20}
			/>
		</Box>
	);
}