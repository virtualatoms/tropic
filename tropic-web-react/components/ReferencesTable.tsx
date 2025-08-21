import { Anchor, Table, Title } from "@mantine/core";

type Props = {
  references: string[]; // formatted references
  dois: string[];
};

export default function ReferencesTable({ references, dois }: Props) {
  const rows = references.map((ref, index) => {
    const doi = dois[index];
    return (
      <Table.Tr key={index}>
        <Table.Td>{ref}</Table.Td>
        <Table.Td>
          <Anchor
            href={`https://doi.org/${doi}`}
            target="_blank"
            underline="always"
          >
            {doi}
          </Anchor>
        </Table.Td>
      </Table.Tr>
    );
  });

  return (
    <div style={{ marginTop: 20 }}>
      <Title order={4} mb="sm">
        References
      </Title>
      <Table striped highlightOnHover withTableBorder>
        <Table.Thead>
          <Table.Tr>
            <Table.Th>Reference</Table.Th>
            <Table.Th>DOI</Table.Th>
          </Table.Tr>
        </Table.Thead>
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>
    </div>
  );
}