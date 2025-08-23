import React from "react";
import { Tabs } from "@mantine/core";
import { IconTable, IconChartDots, IconBook } from "@tabler/icons-react";
import { processRefReaction, References } from "./References";
import { AnalysisChart, processChartReaction } from "./AnalysisChart";
import { ReactionState } from "@/lib/types";
import { reactionToCSV } from "@/lib/export";
import { ExportControls } from "./ExportControls";
import ReactionSearchTable from "./ReactionSearchTable";


export default function ReactionTabsView({ data, isLoading }: ReactionState) {
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
			fileName="reactions"
			onConvertToCSV={reactionToCSV}
		/>
        <ReactionSearchTable data={data} />
      </Tabs.Panel>

      <Tabs.Panel value="analysis" pt={10}>
		<AnalysisChart
		data={data}
		isLoading={isLoading}
		onProcessData={processChartReaction}
		/>
      </Tabs.Panel>

      <Tabs.Panel value="references" pt={10}>
		<References
		data={data}
		isLoading={isLoading}
		onProcessData={processRefReaction}
		/>
      </Tabs.Panel>
    </Tabs>
  );
}
