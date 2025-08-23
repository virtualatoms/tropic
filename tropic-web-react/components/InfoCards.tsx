import { Card, CardSection, Text, Table, Anchor, Center } from "@mantine/core";
import {
  Monomer,
  Reaction,
  Product,
  Parameters,
  Thermo,
  Metadata,
} from "@/lib/types";

const InfoRow = ({
  label,
  value,
}: {
  label: string;
  value: React.ReactNode;
}) => {
  if (value === null || value === undefined || value === "") {
    value = "";
  }
  return (
    <Table.Tr>
      <Table.Td>
        <Text fw="bold" fz="sm">
          {label}
        </Text>
      </Table.Td>
      <Table.Td>
        <Text fz="sm">{value}</Text>
      </Table.Td>
    </Table.Tr>
  );
};

export function MonomerInfoCard({ data }: { data: Monomer }) {
  const rows = [
    <InfoRow key="smiles" label="SMILES" value={data.smiles} />,
    <InfoRow key="inchi" label="InChI" value={data.inchi} />,
    <InfoRow key="iupac" label="IUPAC Name" value={data.iupac_name} />,
    <InfoRow
      key="mw"
      label="Molecular Weight"
      value={
        data.molecular_weight
          ? `${data.molecular_weight.toFixed(2)} g/mol`
          : null
      }
    />,
    <InfoRow key="ring_size" label="Ring Size" value={data.ring_size} />,
    <InfoRow
      key="cid"
      label="PubChem CID"
      value={
        data.pubchem_cid ? (
          <Anchor
            href={`https://pubchem.ncbi.nlm.nih.gov/compound/${data.pubchem_cid}`}
            target="_blank"
            fz="sm"
          >
            {data.pubchem_cid}
          </Anchor>
        ) : null
      }
    />,
  ];

  return (
    <Card withBorder shadow="sm" radius="md">
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Monomer Information</Text>
        </Center>
      </CardSection>
      <Table mt="xs">
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>
    </Card>
  );
}

export function ProductInfoCard({ data }: { data: Product }) {
  const rows = [
    <InfoRow key="smiles" label="SMILES" value={data.smiles} />,
    <InfoRow
      key="repeating_units"
      label="Repeating Units"
      value={data.repeating_units}
    />,
    <InfoRow
      key="deg_of_poly"
      label="Degree of Polymerization"
      value={data.deg_of_poly}
    />,
    <InfoRow key="dispersity" label="Dispersity (Đ)" value={data.dispersity} />,
    <InfoRow
      key="n_avg_molar_mass"
      label="Number-Avg Molar Mass"
      value={data.n_avg_molar_mass}
    />,
    <InfoRow
      key="m_avg_molar_mass"
      label="Mass-Avg Molar Mass"
      value={data.m_avg_molar_mass}
    />,
  ];

  return (
    <Card withBorder shadow="sm" radius="md">
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Product Information</Text>
        </Center>
      </CardSection>
      <Table mt="xs">
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>
    </Card>
  );
}

export function ParametersInfoCard({ data }: { data: Parameters }) {
  const rows = [
    <InfoRow
      key="is_experimental"
      label="Is Experimental"
      value={data.is_experimental ? "Yes" : "No"}
    />,
    <InfoRow
      key="temperature"
      label="Temperature (K)"
      value={data.temperature}
    />,
    <InfoRow key="pressure" label="Pressure (atm)" value={data.pressure} />,
    <InfoRow
      key="state_summary"
      label="State (Monomer-Polymer)"
      value={data.state_summary}
    />,
    <InfoRow
      key="initiator_smiles"
      label="Initiator SMILES"
      value={data.initiator_smiles}
    />,
    <InfoRow key="medium" label="Medium" value={data.medium} />,
    <InfoRow key="solvent" label="Solvent" value={data.solvent} />,
    <InfoRow key="method" label="Method" value={data.method} />,
    <InfoRow
      key="method_calc"
      label="Calculation Method"
      value={data.method_calc}
    />,
    <InfoRow key="functional" label="DFT Functional" value={data.functional} />,
    <InfoRow key="basis_set" label="Basis Set" value={data.basis_set} />,
  ];

  return (
    <Card withBorder shadow="sm" radius="md" mt={25}>
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Parameters</Text>
        </Center>
      </CardSection>
      <Table mt="xs">
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>
    </Card>
  );
}

export function ThermoInfoCard({ data }: { data: Thermo }) {
  const formatValueWithUncertainty = (
    value: number | null,
    std: number | null,
    unit: string,
  ) => {
    if (value === null) return null;
    let formattedValue = value.toFixed(2);
    if (std !== null) {
      formattedValue += ` ± ${std.toFixed(2)}`;
    }
    return `${formattedValue} ${unit}`;
  };

  const rows = [
    <InfoRow
      key="delta_h"
      label="ΔH"
      value={formatValueWithUncertainty(
        data.delta_h,
        data.delta_h_std,
        "kJ/mol",
      )}
    />,
    <InfoRow
      key="delta_s"
      label="ΔS"
      value={formatValueWithUncertainty(
        data.delta_s,
        data.delta_s_std,
        "J/mol·K",
      )}
    />,
    <InfoRow
      key="delta_g"
      label="ΔG"
      value={data.delta_g !== null ? `${data.delta_g.toFixed(2)} kJ/mol` : null}
    />,
    <InfoRow
      key="ceiling_temperature"
      label="Ceiling Temperature"
      value={
        data.ceiling_temperature !== null
          ? `${data.ceiling_temperature.toFixed(2)} °C`
          : null
      }
    />,
  ];

  return (
    <Card withBorder shadow="sm" radius="md" mt={25}>
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Thermodynamics</Text>
        </Center>
      </CardSection>
      <Table mt="xs">
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>
    </Card>
  );
}

export function MetadataInfoCard({ data }: { data: Metadata }) {
  const rows = [
    <InfoRow
      key="doi"
      label="Reference"
      value={
        data.doi ? (
          <Anchor href={`https://doi.org/${data.doi}`} target="_blank" fz="sm">
            {data.formatted_reference || data.doi}
          </Anchor>
        ) : (
          data.formatted_reference
        )
      }
    />,
    <InfoRow key="year" label="Year" value={data.year} />,
    <InfoRow key="comment" label="Comment" value={data.comment} />,
    <InfoRow key="flag" label="Flag" value={data.flag} />,
  ];

  return (
    <Card withBorder shadow="sm" radius="md" mt={25}>
      <CardSection withBorder inheritPadding py="xs">
        <Center>
          <Text fw={700}>Metadata</Text>
        </Center>
      </CardSection>
      <Table mt="xs">
        <Table.Tbody>{rows}</Table.Tbody>
      </Table>
    </Card>
  );
}
