import { Card, CardSection, Center, Text, Table, Anchor } from "@mantine/core";

interface MonomerData {
	monomer: {
		smiles: string;
		inchi: string;
		iupac_name: string;
		molecular_weight: number;
		ring_size: string;
		pubchem_cid?: string | null;
	};
}

export function MonomerInfoCard({ data }: { data: MonomerData }) {
	const { monomer } = data;

	const tableData = [
		["SMILES", monomer.smiles],
		["InChI", monomer.inchi],
		["IUPAC Name", monomer.iupac_name],
		["Molecular Weight", `${monomer.molecular_weight.toFixed(2)} g/mol`],
		["Ring Size", monomer.ring_size],
		[
			"CID",
			monomer.pubchem_cid ? (
				<Anchor
					href={`https://pubchem.ncbi.nlm.nih.gov/compound/${monomer.pubchem_cid}`}
					target="_blank"
					fz="sm"
				>
					{monomer.pubchem_cid}
				</Anchor>
			) : (
				"N/A"
			),
		],
	];
	const rows = tableData.map(([label, value]) => (
		<Table.Tr key={label}>
			<Table.Td>
				<Text fw="bold" fz="sm">{label}</Text>
			</Table.Td>
			<Table.Td><Text fz="sm">{value}</Text></Table.Td>
		</Table.Tr>
	));

	return (
		<Card withBorder shadow="sm" radius="md" pb={10}>
			<CardSection
				withBorder
				inheritPadding
				py="xs"
				bg="var(--mantine-color-gray-0)"
			>
				<Center>
					<Text fw={700}>Monomer Information</Text>
				</Center>
			</CardSection>

			<Table mt={10}>
				<Table.Tbody>{rows}</Table.Tbody>
			</Table>
		</Card>
	);
}
