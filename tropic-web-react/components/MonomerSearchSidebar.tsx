import { MonomerFilters } from "@/lib/types";
import {
  Card,
  CardSection,
  Text,
  RangeSlider,
  MultiSelect,
  Radio,
  RadioGroup,
  Space,
  Group,
} from "@mantine/core";


interface MonomerSearchSidebarProps {
    filters: MonomerFilters;
    setFilters: (filters: MonomerFilters) => void;
}

const radioData = [
    { value: "yes", label: "Yes" },
    { value: "no", label: "No" },
    { value: "both", label: "Both" },
];

export default function MonomerSearchSidebar({ filters, setFilters }: MonomerSearchSidebarProps) {
    return (
        <Card withBorder shadow="sm" radius="md">
            <CardSection withBorder py="xs" inheritPadding>
                <Text fw={700}>Filters</Text>
            </CardSection>

            <Space h={10} />

            <Text size="sm" fw={500}>
                Ring Size
            </Text>
            <RangeSlider
                min={1}
                max={15}
                value={filters.ringSize}
                onChange={(value) => setFilters({ ...filters, ringSize: value })}
                marks={[
                    { value: 1, label: "1" },
                    { value: 5, label: "5" },
                    { value: 10, label: "10" },
                    { value: 15, label: "15+" },
                ]}
                mb={25}
                mt={10}
            />

            <Space h={10} />
            <Text size="sm" fw={500}>
                Molecular Weight (g/mol)
            </Text>
            <RangeSlider
                min={10}
                max={500}
                minRange={10}
                value={filters.molWeight}
                onChange={(value) => setFilters({ ...filters, molWeight: value })}
                marks={[10, 100, 200, 300, 400, 500].map((v) => ({
                    value: v,
                    label: v.toString(),
                }))}
                mb={25}
                mt={10}
            />

            <Space h={10} />
            <MultiSelect
                label="Functional Groups"
                placeholder="Select functional groups"
                data={[
                    { value: "ester", label: "Ester" },
                    { value: "aldehyde", label: "Aldehyde" },
                    { value: "ketone", label: "Ketone" },
                    { value: "amide", label: "Amide" },
                    { value: "ether", label: "Ether" },
                    { value: "amine", label: "Amine" },
                ]}
                value={filters.functionalGroups}
                onChange={(value) => setFilters({ ...filters, functionalGroups: value })}
                clearable
            />

            <Space h={10} />
            <RadioGroup
                label="Has Experimental Data"
                size="sm"
                value={filters.hasExp}
                onChange={(value) => setFilters({ ...filters, hasExp: value })}
            >
                <Group my={10}>
                    {radioData.map((r) => (
                        <Radio key={r.value} label={r.label} value={r.value} />
                    ))}
                </Group>
            </RadioGroup>

            <Space h={10} />
            <RadioGroup
                label="Has Computational Data"
                size="sm"
                value={filters.hasComp}
                onChange={(value) => setFilters({ ...filters, hasComp: value })}
            >
                <Group my={10}>
                    {radioData.map((r) => (
                        <Radio key={r.value} label={r.label} value={r.value} />
                    ))}
                </Group>
            </RadioGroup>
        </Card>
    );
}