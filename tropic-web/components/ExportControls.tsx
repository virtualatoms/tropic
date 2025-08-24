import React, { useState } from "react";
import { Button, Select, Group, Text } from "@mantine/core";
import { MonomerSummary } from "@/lib/types";
import { monomerSummaryToCSV, triggerDownload } from "@/lib/export";

type ExportFormat = "csv" | "json";

// Define the props for the new generic component
type ExportControlsProps<T> = {
  data: T[];
  fileName: string; // e.g., "monomer_summary" or "reactions"
  onConvertToCSV: (data: T[]) => string; // The specific CSV conversion function
};

export function ExportControls<T>({
  data,
  fileName,
  onConvertToCSV,
}: ExportControlsProps<T>) {
  const [format, setFormat] = useState<ExportFormat>("csv");
  const [loading, setLoading] = useState(false);

  const handleExport = () => {
    setLoading(true);
    try {
      let fileContent: string;
      let mimeType: string;

      if (format === "csv") {
        fileContent = onConvertToCSV(data);
        mimeType = "text/csv";
      } else {
        // Create a deep copy and remove the svg property before stringifying
        const dataWithoutSVG = data.map((item: any) => {
          if (item.monomer && "svg" in item.monomer) {
            const { svg, ...restOfMonomer } = item.monomer;
            return { ...item, monomer: restOfMonomer };
          }
          return item;
        });
        fileContent = JSON.stringify(dataWithoutSVG, null, 2);
        mimeType = "application/json";
      }

      triggerDownload(fileContent, `${fileName}.${format}`, mimeType);
    } catch (error) {
      console.error("Export failed:", error);
    } finally {
      setLoading(false);
    }
  };

  return (
    <Group justify="flex-end" mb="md">
      <Text size="sm">Export as:</Text>
      <Select
        data={[
          { value: "csv", label: "CSV" },
          { value: "json", label: "JSON" },
        ]}
        value={format}
        onChange={(value) => setFormat(value as ExportFormat)}
        w={120}
      />
      <Button
        variant="outline"
        onClick={handleExport}
        loading={loading}
        disabled={data.length === 0}
      >
        Export
      </Button>
    </Group>
  );
}
