import React, { useState, useEffect } from 'react';
import { Table, Anchor, Text, LoadingOverlay } from '@mantine/core';
import { MonomerFilters } from '@/lib/types';
import { fetchMonomerSummaries } from '@/lib/api';

interface Reference {
	doi: string;
	formatted: string;
}

export function References({ filters }: { filters: MonomerFilters }) {
	const [references, setReferences] = useState<Reference[]>([]);
	const [loading, setLoading] = useState(true);

	useEffect(() => {
		const fetchReferences = async () => {
			try {
                const monomerSummaries = await fetchMonomerSummaries(filters);

				const uniqueRefs = new Map<string, string>();
				for (const monomerSummary of monomerSummaries) {
                    for (const result of monomerSummary.data || []) {
                        const doi = result.doi;
						const formattedRef = result.formatted_reference;
						if (doi && formattedRef && !uniqueRefs.has(doi)) {
							uniqueRefs.set(doi, formattedRef);
						}
					}
				}

				const refsArray = Array.from(uniqueRefs, ([doi, formatted]) => ({
					doi,
					formatted,
				}));
				setReferences(refsArray);
			} catch (error) {
				console.error("Failed to fetch references:", error);
				setReferences([]);
			} finally {
				setLoading(false);
			}
		};

		fetchReferences();
	}, [filters]);

	const rows = references.map((ref) => (
		<Table.Tr key={ref.doi}>
			<Table.Td>{ref.formatted}</Table.Td>
			<Table.Td>
				<Anchor href={`https://doi.org/${ref.doi}`} target="_blank" size="sm">
					{ref.doi}
				</Anchor>
			</Table.Td>
		</Table.Tr>
	));

	return (
		<div style={{ position: 'relative' }}>
			<LoadingOverlay visible={loading} zIndex={1000} overlayProps={{ radius: "sm", blur: 2 }} />
			<Table mt={20}>
				<Table.Thead>
					<Table.Tr>
						<Table.Th>Reference</Table.Th>
						<Table.Th>DOI</Table.Th>
					</Table.Tr>
				</Table.Thead>
				<Table.Tbody>
                    {rows.length > 0 ? rows : (
                        <Table.Tr>
                            <Table.Td colSpan={2}>
                                <Text c="dimmed" align="center">
                                    No references found for the current filters.
                                </Text>
                            </Table.Td>
                        </Table.Tr>
                    )}
                </Table.Tbody>
			</Table>
		</div>
	);
}
