import { ReactionFilters } from "@/lib/types";
import {
  Card,
  CardSection,
  Text,
  RangeSlider,
  MultiSelect,
  Radio,
  Group,
  Space,
  TextInput,
} from "@mantine/core";

interface ReactionSearchSidebarProps {
  filters: ReactionFilters;
  setFilters: (filters: ReactionFilters) => void;
}

const experimentalRadioData = [
  { value: "yes", label: "Yes" },
  { value: "no", label: "No" },
  { value: "both", label: "Both" },
];

const mediumOptions = [
  { value: "bulk", label: "Bulk" },
  { value: "solvent", label: "Solvent" },
  { value: "gas", label: "Gas" },
];

const reactionTypeOptions = [
  { value: "ROR", label: "ROR" },
  { value: "RER", label: "RER" },
  { value: "ROP", label: "ROP" },
  { value: "RCE", label: "RCE" },
  { value: "MIR", label: "MIR" },
  { value: "unspecified", label: "Unspecified" },
];

const methodOptions = [
    { value: "dft", label: "DFT" },
    { value: "adv_dft", label: "Advanced DFT" },
    { value: "ffmd", label: "FF MD" },
    { value: "aimd", label: "AI MD" },
    { value: "mlmd", label: "ML MD" },
    { value: "xtb", label: "xTB" },
    { value: "ml", label: "ML" },
    { value: "vant_hoff", label: "Van't Hoff" },
    { value: "dsc", label: "DSC" },
    { value: "nmr", label: "NMR" },
    { value: "calorimetry", label: "Calorimetry" },
    { value: "other", label: "Other" },
    { value: "unspecified", label: "Unspecified" },
];

export default function ReactionSearchSidebar({
  filters,
  setFilters,
}: ReactionSearchSidebarProps) {
  const currentYear = new Date().getFullYear();

  return (
    <Card withBorder shadow="sm" radius="md">
      <CardSection withBorder py="xs" inheritPadding>
        <Text fw={700}>Reaction Filters</Text>
      </CardSection>

      <Space h="md" />
      <TextInput
        label="Initiator SMILES"
        placeholder="e.g., CC(=O)OC"
        value={filters.initiator}
        onChange={(event) =>
          setFilters({ ...filters, initiator: event.currentTarget.value })
        }
      />

      <Space h="md" />
      <MultiSelect
        label="Reaction Type"
        placeholder="Select types"
        data={reactionTypeOptions}
        value={filters.reactionType}
        onChange={(value) => setFilters({ ...filters, reactionType: value })}
        clearable
      />

      <Space h="md" />
      <MultiSelect
        label="Medium"
        placeholder="Select mediums"
        data={mediumOptions}
        value={filters.medium}
        onChange={(value) => setFilters({ ...filters, medium: value })}
        clearable
      />

      <Space h="md" />
      <MultiSelect
        label="Method"
        placeholder="Select methods"
        data={methodOptions}
        value={filters.method}
        onChange={(value) => setFilters({ ...filters, method: value })}
        clearable
      />

      <Space h="md" />
      <Text size="sm" fw={500}>
        Year of Publication
      </Text>
      <RangeSlider
        min={1950}
        max={currentYear}
        value={filters.year}
        onChange={(value) => setFilters({ ...filters, year: value })}
        marks={[
          { value: 1950, label: "1950" },
          { value: 1980, label: "1980" },
          { value: 2010, label: "2010" },
          { value: currentYear, label: `${currentYear}` },
        ]}
        mb="xl"
        mt="xs"
      />

      <Radio.Group
        label="Is Experimental"
        size="sm"
        value={filters.isExperimental}
        onChange={(value) => setFilters({ ...filters, isExperimental: value as any })}
      >
        <Group my="xs">
          {experimentalRadioData.map((r) => (
            <Radio key={r.value} label={r.label} value={r.value} />
          ))}
        </Group>
      </Radio.Group>
    </Card>
  );
}