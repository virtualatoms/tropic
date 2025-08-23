import React from "react";
import { Tabs } from "@mantine/core";
import { IconTable, IconChartDots, IconBook } from "@tabler/icons-react";
import { processRefMonomerSummary, References } from "./References";
import { AnalysisChart, processChartMonomerSummary } from "./AnalysisChart";
import MonomerSearchTable from "./MonomerSearchTable";
import { MonomerSummaryState } from "@/lib/types";
import { monomerSummaryToCSV } from "@/lib/export";
import { ExportControls } from "./ExportControls";


export default function MonomerTabsView({ data, isLoading }: MonomerSummaryState) {
  return (
    <Tabs defaultValue="table">
      <Tabs.List>
        <Tabs.Tab value="table" leftSection={<IconTable size={16} />}>
          Table
        </Tabs.Tab>
        <Tabs.Tab value="analysis" leftSection={<IconChartDots size={16} />}>
          Analysis
        </Tabs.Tab>
        <Tabs.Tab value="references" leftSection={<IconBook size={16} />}>
          References
        </Tabs.Tab>
      </Tabs.List>

      <Tabs.Panel value="table" pt={10}>
		<ExportControls
			data={data}
			fileName="monomer_summary"
			onConvertToCSV={monomerSummaryToCSV}
		/>
        <MonomerSearchTable data={data} />
      </Tabs.Panel>

      <Tabs.Panel value="analysis" pt={10}>
		<AnalysisChart
		data={data}
		isLoading={isLoading}
		onProcessData={processChartMonomerSummary}
		/>
      </Tabs.Panel>

      <Tabs.Panel value="references" pt={10}>
		<References
		data={data}
		isLoading={isLoading}
		onProcessData={processRefMonomerSummary}
		/>
      </Tabs.Panel>
    </Tabs>
  );
}
