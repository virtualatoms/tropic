import React from 'react';
import { Tabs } from '@mantine/core';
import { IconTable, IconChartDots, IconBook } from '@tabler/icons-react';
import { References } from './References';
import { MonomerFilters } from '@/lib/types';

export default function TabsView({filters}: {filters: MonomerFilters;}) {
	return (
		<Tabs defaultValue="table" color="blue">
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
				<h1>Table</h1>
			</Tabs.Panel>

			<Tabs.Panel value="analysis" pt={10}>
				<h1>Analysis</h1>
			</Tabs.Panel>

			<Tabs.Panel value="references" pt={10}>
				<References filters={filters} />
			</Tabs.Panel>
		</Tabs>
	);
}