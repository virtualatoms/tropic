import { ScatterChart } from "@mantine/charts";
import { Select, Group, Text } from "@mantine/core";
import { useMemo, useState } from "react";

type ChartDatum = {
  ring_size: number;
  delta_h: number;
  delta_s: number;
  delta_g: number;
  ceiling_temperature: number;
};

type Props = {
  experimental: ChartDatum[];
  computational: ChartDatum[];
};

const yOptions = [
  { value: "delta_h", label: "ΔH" },
  { value: "delta_s", label: "ΔS" },
  { value: "delta_g", label: "ΔG" },
  { value: "ceiling_temperature", label: "Ceiling Temperature" },
];

export default function AnalysisChart({
  experimental,
  computational,
}: Props) {
  const [yField, setYField] = useState<keyof ChartDatum>("delta_h");

  const formattedData = useMemo(() => {
    return [
      {
        name: "Computational",
        color: "red.5",
        data: computational.map((d) => ({
          x: d.ring_size,
          y: d[yField],
        })),
      },
      {
        name: "Experimental",
        color: "blue.5",
        data: experimental.map((d) => ({
          x: d.ring_size,
          y: d[yField],
        })),
      },
    ];
  }, [experimental, computational, yField]);

  return (
    <div>
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
        pointLabels="xy"
        tickLine="xy"
        withLegend
        id="analysis-chart"
      />
    </div>
  );
}

function yAxisLabel(key: keyof ChartDatum) {
  switch (key) {
    case "delta_h":
      return "ΔH / kJ/mol";
    case "delta_s":
      return "ΔS / J/mol·K";
    case "delta_g":
      return "ΔG / kJ/mol";
    case "ceiling_temperature":
      return "Ceiling Temp / °C";
    default:
      return key;
  }
}