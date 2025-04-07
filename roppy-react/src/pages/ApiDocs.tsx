import { Title, Text, Stack, Paper, Code, Table } from '@mantine/core';

export default function ApiDocs() {
  const endpoints = [
    {
      method: 'GET',
      path: '/api/monomers',
      description: 'Get all monomers',
      response: `[
  {
    "id": "string",
    "name": "string",
    "formula": "string",
    "molecular_weight": number
  }
]`,
    },
    {
      method: 'GET',
      path: '/api/monomers/{id}',
      description: 'Get a specific monomer by ID',
      response: `{
  "id": "string",
  "name": "string",
  "formula": "string",
  "molecular_weight": number,
  "smiles": "string",
  "description": "string"
}`,
    },
  ];

  return (
    <Stack>
      <Title order={1}>API Documentation</Title>

      <Paper p="md" withBorder>
        <Title order={2}>Base URL</Title>
        <Text>
          All API endpoints are relative to the base URL:
          <Code block>https://api.roppy.com/v1</Code>
        </Text>
      </Paper>

      <Paper p="md" withBorder>
        <Title order={2}>Authentication</Title>
        <Text>
          The API uses API keys for authentication. Include your API key in the request header:
          <Code block>Authorization: Bearer your-api-key</Code>
        </Text>
      </Paper>

      <Paper p="md" withBorder>
        <Title order={2}>Endpoints</Title>
        <Table>
          <Table.Thead>
            <Table.Tr>
              <Table.Th>Method</Table.Th>
              <Table.Th>Path</Table.Th>
              <Table.Th>Description</Table.Th>
              <Table.Th>Response</Table.Th>
            </Table.Tr>
          </Table.Thead>
          <Table.Tbody>
            {endpoints.map((endpoint) => (
              <Table.Tr key={endpoint.path}>
                <Table.Td>
                  <Code>{endpoint.method}</Code>
                </Table.Td>
                <Table.Td>
                  <Code>{endpoint.path}</Code>
                </Table.Td>
                <Table.Td>{endpoint.description}</Table.Td>
                <Table.Td>
                  <Code block>{endpoint.response}</Code>
                </Table.Td>
              </Table.Tr>
            ))}
          </Table.Tbody>
        </Table>
      </Paper>
    </Stack>
  );
} 