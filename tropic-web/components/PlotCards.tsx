import { Card, CardSection, Text, Center, Box } from "@mantine/core";
import { Vanthoff, ComputationalExtrapolation } from "@/lib/types";
import { useMemo } from "react";

import { Scatter } from "react-chartjs-2";
import {
  Chart as ChartJS,
  LinearScale,
  PointElement,
  LineElement,
  Tooltip,
  Legend,
} from "chart.js";
import annotationPlugin from "chartjs-plugin-annotation";

ChartJS.register(
  LinearScale,
  PointElement,
  LineElement,
  Tooltip,
  Legend,
  annotationPlugin,
);

export function ExtrapolationPlotCard({
  data,
}: {
  data: ComputationalExtrapolation | null;
}) {
  const chartJSData = useMemo(() => {
    if (!data?.inverse_repeating_units || !data?.delta_h) {
      return { datasets: [] };
    }
    const scatterDataset = {
      label: "Data Points",
      type: "scatter" as const,
      data: data.inverse_repeating_units.map((inv_n, i) => ({
        x: inv_n,
        y: data.delta_h![i],
      })),
      backgroundColor: "rgba(54, 162, 235, 0.6)",
    };
    const trendlineDataset = {
      label: "Trendline",
      type: "scatter" as const,
      data: [] as { x: number; y: number }[],
      borderColor: "rgba(255, 99, 132, 1)",
      borderWidth: 2,
      pointRadius: 0,
      fill: false,
      showLine: true,
    };
    if (typeof data.slope === "number" && typeof data.intercept === "number") {
      const xMax = Math.max(...data.inverse_repeating_units);
      trendlineDataset.data = [
        { x: 0, y: data.intercept },
        { x: xMax, y: data.slope * xMax + data.intercept },
      ];
    }
    return { datasets: [scatterDataset, trendlineDataset] };
  }, [data]);

  const options = {
    responsive: true,
    scales: {
      x: {
        type: "linear" as const,
        position: "bottom" as const,
        title: {
          display: true,
          text: "1 / Repeating Units",
          font: { size: 14 },
        },
        ticks: {
          font: { size: 14 },
        },
      },
      y: {
        title: { display: true, text: "ΔH (kJ/mol)", font: { size: 14 } },
        ticks: {
          font: { size: 14 },
        },
      },
    },
    plugins: {
      legend: {
        position: "top" as const,
        labels: { font: { size: 14 } },
      },
      tooltip: {
        callbacks: {
          label: function (context: any) {
            if (context.dataset.label !== "Data Points") return "";
            const xVal = context.raw.x;
            const yVal = context.raw.y;
            const repeatingUnits = 1 / xVal;
            return [
              `Repeating Units: ${repeatingUnits.toFixed(2)}`,
              `1 / Units: ${xVal.toFixed(4)}`,
              `Energy: ${yVal.toFixed(2)} kJ/mol`,
            ];
          },
        },
      },
      annotation: {
        annotations: {
          interceptLine: {
            type: "line" as const,
            yMin: typeof data?.intercept === "number" ? data.intercept : 0,
            yMax: typeof data?.intercept === "number" ? data.intercept : 0,
            borderColor: "rgb(75, 192, 192)",
            borderWidth: 2,
            borderDash: [6, 6],
          },
          summaryLabel: {
            type: "label" as const,
            content: [
              `Intercept (at x=0): ${data?.intercept?.toFixed(2) ?? "N/A"}`,
              `Std. Error: ${data?.std_err?.toFixed(2) ?? "N/A"}`,
            ],
            xValue: (context: any) => context.chart.scales.x.max,
            yValue: (context: any) => context.chart.scales.y.max,
            xAdjust: -100,
            yAdjust: 50,
            textAlign: "left" as const,
            backgroundColor: "rgba(245, 245, 245, 0.8)",
            color: "black",
            font: { size: 14 },
            padding: 6,
            borderRadius: 6,
          },
        },
      },
    },
  };
  return (
    <Card withBorder shadow="sm" radius="md">
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Computational Extrapolation</Text>
        </Center>
      </CardSection>
      <Box p="md" style={{ minHeight: 350 }}>
        <Scatter data={chartJSData} options={options} />
      </Box>
    </Card>
  );
}

export function VanthoffPlotCard({ data }: { data: Vanthoff | null }) {
  const chartJSData = useMemo(() => {
    if (
      !data?.inverse_temperature ||
      !data?.r_ln_equilibrium_concentration ||
      !data?.temperature
    ) {
      return { datasets: [] };
    }

    const scatterDataset = {
      label: "Data Points",
      data: data.inverse_temperature.map((inv_t, i) => ({
        x: inv_t,
        y: data.r_ln_equilibrium_concentration![i],
        temp: data.temperature![i], // Store original temp for tooltip
      })),
      backgroundColor: "rgba(54, 162, 235, 0.6)",
    };

    return {
      datasets: [scatterDataset],
    };
  }, [data]);

  const options = {
    responsive: true,
    scales: {
      x: {
        type: "linear" as const,
        position: "bottom" as const,
        title: {
          display: true,
          text: "1/T (K⁻¹)",
          font: { size: 14 },
        },
        ticks: { font: { size: 14 } },
      },
      y: {
        title: {
          display: true,
          text: "R·ln([M]e)",
          font: { size: 14 },
        },
        ticks: { font: { size: 14 } },
      },
    },
    plugins: {
      legend: {
        position: "top" as const,
        labels: { font: { size: 14 } },
      },
      // Custom tooltip configuration
      tooltip: {
        callbacks: {
          label: function (context: any) {
            const xVal = context.raw.x;
            const yVal = context.raw.y;
            const temp = context.raw.temp; // Retrieve original temperature

            return [
              `Temperature: ${temp.toFixed(2)} K`,
              `1/T: ${xVal.toFixed(4)} K⁻¹`,
              `R·ln([M]e): ${yVal.toFixed(2)}`,
            ];
          },
        },
      },
    },
  };

  return (
    <Card withBorder shadow="sm" radius="md">
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Van&apos;t Hoff Analysis</Text>
        </Center>
      </CardSection>
      <Box p="md" style={{ minHeight: 350 }}>
        <Scatter data={chartJSData} options={options} />
      </Box>
    </Card>
  );
}
