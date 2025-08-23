import React, { useState } from "react";
import { Button, Select, Group, Text } from "@mantine/core";
import { MonomerSummary } from "@/lib/types";
import { monomerSummaryToCSV, triggerDownload } from "@/lib/export";

type ExportFormat = "csv" | "json";

export function ExportControls({ data }: { data: MonomerSummary[] }) {
  const [format, setFormat] = useState<ExportFormat>("csv");
  const [loading, setLoading] = useState(false);

  const handleExport = () => {
    setLoading(true);
    try {
      let fileContent: string;
      let mimeType: string;

      if (format === "csv") {
        fileContent = monomerSummaryToCSV(data);
        mimeType = "text/csv";
      } else {
        fileContent = JSON.stringify(data, null, 2);
        mimeType = "application/json";
      }

      triggerDownload(fileContent, `monomer_summary.${format}`, mimeType);
    } catch (error) {
      console.error("Export failed:", error);
      // TODO: Show a notification
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
        disabled={data.length === 0} // Disable button if there's no data
      >
        Export
      </Button>
    </Group>
  );
}
