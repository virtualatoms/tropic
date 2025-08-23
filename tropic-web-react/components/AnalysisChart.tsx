import { Select, Group, Text, LoadingOverlay, Box } from "@mantine/core";
import { useMemo, useState } from "react";
import { useRouter } from "next/router";
import { MonomerSummary, Reaction } from "@/lib/types";
import { Scatter } from "react-chartjs-2";
import {
  Chart as ChartJS,
  LinearScale,
  PointElement,
  LineElement,
  Tooltip,
  Legend,
} from "chart.js";

ChartJS.register(LinearScale, PointElement, LineElement, Tooltip, Legend);


export const processChartMonomerSummary = (data: MonomerSummary[]) => {
  const expData: ChartDatum[] = [];
  const compData: ChartDatum[] = [];

  for (const monomerSummary of data) {
    for (const d of monomerSummary.data) {
      if (d.delta_h === null || d.delta_s === null) continue;
      const datum = {
        ring_size: monomerSummary.monomer.ring_size,
        delta_h: d.delta_h,
        delta_s: d.delta_s,
        ceiling_temperature: d.ceiling_temperature,
        svg: monomerSummary.monomer.svg,
        id: monomerSummary.monomer.monomer_id,
        pageType: 'monomer',
      };
      if (d.is_experimental) {
        expData.push(datum);
      } else {
        compData.push(datum);
      }
    }
  }
  return { experimental: expData, computational: compData };
};

export const processChartReaction = (data: Reaction[]) => {
  const expData: ChartDatum[] = [];
  const compData: ChartDatum[] = [];

  for (const reaction of data) {
    if (reaction.thermo.delta_h === null || reaction.thermo.delta_s === null) continue;
    const datum = {
        ring_size: reaction.monomer.ring_size,
        delta_h: reaction.thermo.delta_h,
        delta_s: reaction.thermo.delta_s,
        ceiling_temperature: reaction.thermo.ceiling_temperature,
        svg: reaction.monomer.svg,
        id: reaction.reaction_id,
        pageType: 'reaction',
    };
    if (reaction.parameters.is_experimental) {
      expData.push(datum);
    } else {
      compData.push(datum);
    }
  }
  return { experimental: expData, computational: compData };
};

type ChartDatum = {
  ring_size: number;
  delta_h: number | null;
  delta_s: number | null;
  ceiling_temperature: number | null;
  svg: string;
  id: string; 
  pageType: 'monomer' | 'reaction';
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

type AnalysisChartProps<T> = {
  data: T[];
  isLoading: boolean;
  onProcessData: (data: T[]) => {
    experimental: ChartDatum[];
    computational: ChartDatum[];
  };
};

const externalTooltipHandler = (context) => {
  const { chart, tooltip } = context;
  let tooltipEl = chart.canvas.parentNode.querySelector("div#chartjs-tooltip");

  if (!tooltipEl) {
    tooltipEl = document.createElement("div");
    tooltipEl.id = "chartjs-tooltip";
    tooltipEl.style.background = "rgba(0, 0, 0, 0.7)";
    tooltipEl.style.borderRadius = "3px";
    tooltipEl.style.color = "white";
    tooltipEl.style.opacity = 1;
    tooltipEl.style.pointerEvents = "none";
    tooltipEl.style.position = "absolute";
    tooltipEl.style.transform = "translate(-50%, 0)";
    tooltipEl.style.transition = "all .1s ease";
    tooltipEl.style.padding = "6px";
    tooltipEl.style.width = "150px";
    chart.canvas.parentNode.appendChild(tooltipEl);
  }

  if (tooltip.opacity === 0) {
    tooltipEl.style.opacity = 0;
    return;
  }

  if (tooltip.body) {
    const dataPoint = tooltip.dataPoints[0].raw;
    const original = dataPoint.originalData;

    let innerHTML = `<img src="data:image/svg+xml;base64,${original.svg}" alt="Monomer" style="width: 100%; height: auto; background: white; border-radius: 2px; margin-bottom: 5px;" />`;
    innerHTML += `<div>Ring Size: ${original.ring_size}</div>`;
    innerHTML += `<div>ΔH: ${original.delta_h?.toFixed(2)} kJ/mol</div>`;

    if (original.delta_s !== null && original.delta_s !== undefined) {
      innerHTML += `<div>ΔS: ${original.delta_s.toFixed(2)} J/mol·K</div>`;
    }

    if (original.ceiling_temperature !== null && original.ceiling_temperature !== undefined) {
      innerHTML += `<div>Tc: ${original.ceiling_temperature.toFixed(2)} °C</div>`;
    }
    tooltipEl.innerHTML = innerHTML;
  }

  const { offsetLeft: positionX, offsetTop: positionY } = chart.canvas;
  tooltipEl.style.opacity = 1;
  tooltipEl.style.left = positionX + tooltip.caretX + "px";
  tooltipEl.style.top = positionY + tooltip.caretY + "px";
};

export function AnalysisChart<T>({
  data,
  isLoading,
  onProcessData,
}: AnalysisChartProps<T>) {
  const [yField, setYField] = useState<keyof ChartDatum>("delta_h");
  const router = useRouter(); 

  const chartData = useMemo(() => {
    return onProcessData(data);
  }, [data, onProcessData]);

  const chartJSData = useMemo(() => {
    return {
      datasets: [
        {
          label: "Experimental",
          data: chartData.experimental.map((d) => ({
            x: d.ring_size,
            y: d[yField],
            originalData: d, 
          })),
          pointRadius: 4,
          backgroundColor: "rgba(54, 162, 235, 0.6)",
        },
        {
          label: "Computational",
          data: chartData.computational.map((d) => ({
            x: d.ring_size,
            y: d[yField],
            originalData: d,
          })),
          pointRadius: 4,
          backgroundColor: "rgba(255, 99, 132, 0.6)",
        },
      ],
    };
  }, [chartData, yField]);

  const options = {
    responsive: true,
    maintainAspectRatio: false,
    scales: {
       x: {
        title: {
          display: true,
          text: "Ring Size",
          font: { size: 14 },
        },
        ticks: { font: { size: 14 } },
        min: 0,
      },
      y: {
        title: {
          display: true,
          text: yAxisLabel(yField),
          font: { size: 14 },
        },
        ticks: {
          font: { size: 14 },
          },
        },
      },
    onClick: (event: any, elements: any[]) => {
      if (elements.length === 0) return; 

      const element = elements[0];
      const dataPoint = chartJSData.datasets[element.datasetIndex].data[element.index];
      const { id, pageType } = dataPoint.originalData;

      if (id && pageType) {
        const url = pageType === 'monomer' ? `/monomers/${id}` : `/reactions/${id}`;
        router.push(url);
      }
    },
    plugins: {
      legend: { position: "top" as const, labels: {font: {size: 14}}},
      tooltip: {
        enabled: false, // Disable the default canvas tooltip
        external: externalTooltipHandler, // Enable our custom HTML tooltip
      },
    },
  };

  return (
    <Box pos="relative">
      <LoadingOverlay visible={isLoading} />
      <Group justify="flex-end" mb="md">
        <Text size="sm">Y axis property:</Text>
        <Select
          data={yOptions}
          value={yField}
          onChange={(val) => setYField(val as keyof ChartDatum)}
          w={220}
        />
      </Group>
      <Box h={500}>
        <Scatter options={options} data={chartJSData} />
      </Box>
    </Box>
  );
}