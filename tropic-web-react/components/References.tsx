import React, { useState, useEffect } from 'react';
import { Table, Anchor, Text, LoadingOverlay } from '@mantine/core';
import { API_ENDPOINT } from '@/lib/constants';
import { MonomerFilters } from '@/lib/types';

interface Reference {
	doi: string;
	formatted: string;
}


const buildQuery = (filters: MonomerFilters): string => {
	const params = new URLSearchParams();

	if (filters.smiles) {
		params.append('search', filters.smiles);
	}

	params.append('monomer__ring_size__gte', filters.ringSize[0].toString());
	if (filters.ringSize[1] < 15) {
		params.append('monomer__ring_size__lte', filters.ringSize[1].toString());
	}

	params.append('monomer__molecular_weight__gte', filters.molWeight[0].toString());
	params.append('monomer__molecular_weight__lte', filters.molWeight[1].toString());

	if (filters.hasComp !== 'both') {
		params.append('has_comp', (filters.hasComp === 'yes').toString());
	}
	if (filters.hasExp !== 'both') {
		params.append('has_exp', (filters.hasExp === 'yes').toString());
	}

    // This assumes your API can handle multiple functional_groups params
	filters.functionalGroups.forEach(group => {
		params.append('functional_groups', group);
	});

	params.append('size', '1000');
	return params.toString();
};

export function References({ filters }: { filters: MonomerFilters }) {
	const [references, setReferences] = useState<Reference[]>([]);
	const [loading, setLoading] = useState(true);

	useEffect(() => {
		const fetchReferences = async () => {
			setLoading(true);
			try {
				const query = buildQuery(filters);
				const response = await fetch(`${API_ENDPOINT}/reactions?${query}`);
				if (!response.ok) {
					throw new Error('Network response was not ok');
				}
				const results = await response.json();

				const uniqueRefs = new Map<string, string>();
				for (const result of results) {
					const doi = result.metadata?.doi;
					const formattedRef = result.metadata?.formatted_reference;
					if (doi && formattedRef && !uniqueRefs.has(doi)) {
						uniqueRefs.set(doi, formattedRef);
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
