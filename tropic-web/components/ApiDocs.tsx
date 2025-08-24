import {
  Title,
  Text,
  Accordion,
  Table,
  Code,
  Stack,
  Divider,
  Badge,
  ThemeIcon,
  Group,
} from "@mantine/core";
import { CodeHighlight } from "@mantine/code-highlight";
import { IconApi, IconBook, IconBrandPython } from "@tabler/icons-react";
import { API_ENDPOINT } from "@/lib/constants";

const ModelFieldTable = ({ fields }) => {
  const rows = fields.map((field) => (
    <Table.Tr key={field.name}>
      <Table.Td>
        <Code fw={700}>{field.name}</Code>
      </Table.Td>
      <Table.Td>
        <Badge variant="outline" color="gray">
          {field.type}
        </Badge>
      </Table.Td>
      <Table.Td>{field.description}</Table.Td>
    </Table.Tr>
  ));

  return (
    <Table striped withTableBorder withColumnBorders>
      <Table.Thead>
        <Table.Tr>
          <Table.Th>Field Name</Table.Th>
          <Table.Th>Type</Table.Th>
          <Table.Th>Description</Table.Th>
        </Table.Tr>
      </Table.Thead>
      <Table.Tbody>{rows}</Table.Tbody>
    </Table>
  );
};

const models = {
  Reaction: [
    {
      name: "reaction_id",
      type: "string",
      description: "Unique display ID for the reaction.",
    },
    {
      name: "type",
      type: "string",
      description: 'Type of reaction (e.g., "ROP", "ROR").',
    },
    {
      name: "monomer",
      type: "Monomer",
      description: "The monomer object for the reaction.",
    },
    {
      name: "product",
      type: "Product",
      description: "The product object for the reaction.",
    },
    {
      name: "parameters",
      type: "Parameters",
      description: "Experimental/computational parameters.",
    },
    {
      name: "thermo",
      type: "Thermo",
      description: "Thermodynamic data and results.",
    },
    {
      name: "metadata",
      type: "Metadata",
      description: "Reaction references and metadata.",
    },
  ],
  Monomer: [
    {
      name: "monomer_id",
      type: "string",
      description: "Unique display ID for the monomer.",
    },
    {
      name: "smiles",
      type: "string",
      description: "SMILES notation of the molecular structure.",
    },
    {
      name: "inchi",
      type: "string | null",
      description: "InChI identifier for the molecule.",
    },
    {
      name: "molecular_weight",
      type: "float | null",
      description: "Molecular weight in g/mol.",
    },
    {
      name: "functional_group",
      type: "string | null",
      description: "Primary functional group.",
    },
    {
      name: "iupac_name",
      type: "string | null",
      description: "IUPAC name of the molecule.",
    },
    {
      name: "pubchem_cid",
      type: "int | null",
      description: "PubChem Compound ID (CID).",
    },
    {
      name: "ring_size",
      type: "int | null",
      description: "Size of the ring in the monomer structure.",
    },
    {
      name: "xyz",
      type: "string | null",
      description: "XYZ coordinates for the molecule.",
    },
  ],
  Product: [
    {
      name: "smiles",
      type: "string | null",
      description: "BigSMILES notation for the polymer.",
    },
    {
      name: "repeating_units",
      type: "float | null",
      description: "Number of repeating units.",
    },
    {
      name: "deg_of_poly",
      type: "float | null",
      description: "Average number of monomer units.",
    },
    {
      name: "dispersity",
      type: "float | null",
      description: "Molar mass distribution (Đ).",
    },
    {
      name: "n_avg_molar_mass",
      type: "float | null",
      description: "Number-average molar mass.",
    },
    {
      name: "m_avg_molar_mass",
      type: "float | null",
      description: "Mass-average molar mass.",
    },
  ],
  Parameters: [
    {
      name: "is_experimental",
      type: "boolean",
      description: "True for experimental, false for computational.",
    },
    {
      name: "temperature",
      type: "float | null",
      description: "Reaction temperature in Kelvin (K).",
    },
    {
      name: "pressure",
      type: "float | null",
      description: "Reaction pressure (if not standard).",
    },
    {
      name: "initiator_smiles",
      type: "string | null",
      description: "SMILES of the initiator.",
    },
    {
      name: "medium",
      type: "string | null",
      description: 'Medium of the reaction (e.g., "bulk", "solvent").',
    },
    {
      name: "method",
      type: "string | null",
      description: "Experimental or computational method used.",
    },
    {
      name: "state_summary",
      type: "string",
      description: "Formatted summary of monomer-polymer states.",
    },
  ],
  Thermo: [
    {
      name: "delta_h",
      type: "float | null",
      description: "Enthalpy of polymerization (kJ/mol).",
    },
    {
      name: "delta_s",
      type: "float | null",
      description: "Entropy of polymerization (J/mol·K).",
    },
    {
      name: "delta_g",
      type: "float | null",
      description: "Free energy of polymerization (kJ/mol).",
    },
    {
      name: "ceiling_temperature",
      type: "float | null",
      description: "Ceiling temperature (K).",
    },
    {
      name: "vanthoff",
      type: "Vanthoff | null",
      description: "Van't Hoff data for the reaction.",
    },
    {
      name: "extrapolation",
      type: "ComputationalExtrapolation | null",
      description: "Computational extrapolation data.",
    },
  ],
  Metadata: [
    { name: "year", type: "int | null", description: "Year of publication." },
    {
      name: "doi",
      type: "string | null",
      description: "Digital Object Identifier (DOI) for the source.",
    },
    {
      name: "formatted_reference",
      type: "string | null",
      description: "Full formatted reference string.",
    },
    {
      name: "comment",
      type: "string | null",
      description: "Additional comments or notes.",
    },
  ],
};

const filterParams = [
  { param: "reaction_id", desc: "Filter by a specific reaction ID." },
  { param: "type", desc: "Filter by reaction type (e.g., ROP, ROR)." },
  {
    param: "monomer__ring_size__gte",
    desc: "Filter by monomer ring size (greater than or equal to).",
  },
  {
    param: "monomer__molecular_weight__lte",
    desc: "Filter by monomer molecular weight (less than or equal to).",
  },
  {
    param: "parameters__is_experimental",
    desc: "Filter by experimental (true) or computational (false) data.",
  },
  {
    param: "parameters__method__in",
    desc: "Filter by one or more methods (e.g., dft, dsc).",
  },
  {
    param: "thermo__delta_h__gte",
    desc: "Filter by enthalpy of polymerization.",
  },
  { param: "metadata__year__gte", desc: "Filter by publication year." },
];

export default function ApiDocs() {
  const curlExample = `curl -X 'GET' '${API_ENDPOINT}/reactions?parameters__is_experimental=true&monomer__ring_size__gte=5' -H 'accept: application/json'`;
  const pythonInstallCode = `pip install tropic-client`;
  const pythonUsageCode = `from tropic.client import TropicClient

with TropicClient() as client:
    # Find all ROR reactions with a monomer ring size >= 10
    reactions = client.get_reactions(
        type="ROR",
        monomer__ring_size__gte=10
    )

    for reaction in reactions:
        print(reaction.reaction_id, reaction.thermo.delta_h)`;

  return (
    <Stack gap="xl">
      <Title order={1}>API Documentation</Title>
      <Text c="dimmed">
        This document provides the information you need to interact with the
        TROPIC Database programmatically. You can use the REST API directly or
        our convenient Python client.
      </Text>

      <Divider />

      <Stack gap="xs">
        <Group>
          <ThemeIcon size="lg" variant="light">
            <IconBook />
          </ThemeIcon>
          <Title order={2} id="data-models">
            Data Models
          </Title>
        </Group>
        <Text>
          The API returns data structured according to the following models. The
          primary model is <Code>Reaction</Code>, which contains several nested
          data objects.
        </Text>
        <Accordion variant="separated">
          {Object.entries(models).map(([name, fields]) => (
            <Accordion.Item key={name} value={name}>
              <Accordion.Control>
                <Text fw={500}>{name}</Text>
              </Accordion.Control>
              <Accordion.Panel>
                <ModelFieldTable fields={fields} />
              </Accordion.Panel>
            </Accordion.Item>
          ))}
        </Accordion>
      </Stack>

      <Divider />

      <Stack gap="xs">
        <Group>
          <ThemeIcon size="lg" variant="light">
            <IconApi />
          </ThemeIcon>
          <Title order={2} id="rest-api-usage">
            REST API Usage
          </Title>
        </Group>
        <Text>
          You can interact with the database by making HTTP requests to the REST
          API endpoints. The base URL for the API is <Code>{API_ENDPOINT}</Code>
          .
        </Text>
        <Title order={4} id="querying-reactions" mt={10}>
          Querying Reactions
        </Title>
        <Text>
          The main endpoint for retrieving data is <Code>/reactions</Code>. You
          can filter the results by appending query parameters to the URL.
          Nested fields are accessed using a double underscore (<Code>__</Code>
          ).
        </Text>
        <Title order={5} mt="sm">
          Example Request:
        </Title>
        <CodeHighlight code={curlExample} language="bash" />
        <Title order={5} mt="sm">
          Common Filter Parameters:
        </Title>
        <Table>
          <Table.Thead>
            <Table.Tr>
              <Table.Th>Parameter</Table.Th>
              <Table.Th>Description</Table.Th>
            </Table.Tr>
          </Table.Thead>
          <Table.Tbody>
            {filterParams.map((p) => (
              <Table.Tr key={p.param}>
                <Table.Td>
                  <Code>{p.param}</Code>
                </Table.Td>
                <Table.Td>{p.desc}</Table.Td>
              </Table.Tr>
            ))}
          </Table.Tbody>
        </Table>
      </Stack>

      <Divider />

      <Stack gap="xs">
        <Group>
          <ThemeIcon size="lg" variant="light">
            <IconBrandPython />
          </ThemeIcon>
          <Title order={2} id="python-client">
            Python Client
          </Title>
        </Group>
        <Text>
          For Python users, we provide a simple client library that handles API
          requests and data parsing for you.
        </Text>
        <Title order={4} id="installation" mt="md">
          Installation
        </Title>
        <CodeHighlight code={pythonInstallCode} language="bash" />
        <Title order={4} id="example-usage" mt="md">
          Usage Example
        </Title>
        <CodeHighlight code={pythonUsageCode} language="python" />
      </Stack>
    </Stack>
  );
}
