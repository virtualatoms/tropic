import React, { useMemo } from "react";
import { Table, Anchor, Text, LoadingOverlay } from "@mantine/core";
import { MonomerSummaryState } from "@/lib/types";

interface Reference {
  doi: string;
  formatted: string;
}

export default function References({ data, isLoading }: MonomerSummaryState) {
  const references = useMemo<Reference[]>(() => {
    const uniqueRefs = new Map<string, string>();

    for (const summary of data) {
      for (const result of summary.data || []) {
        const doi = result.doi;
        const formattedRef = result.formatted_reference;
        if (doi && formattedRef && !uniqueRefs.has(doi)) {
          uniqueRefs.set(doi, formattedRef);
        }
      }
    }
    return Array.from(uniqueRefs, ([doi, formatted]) => ({ doi, formatted }));
  }, [data]);

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
    <div style={{ position: "relative" }}>
      <LoadingOverlay
        visible={isLoading}
        zIndex={1000}
        overlayProps={{ radius: "sm", blur: 2 }}
      />
      <Table mt={20}>
        <Table.Thead>
          <Table.Tr>
            <Table.Th>Reference</Table.Th>
            <Table.Th>DOI</Table.Th>
          </Table.Tr>
        </Table.Thead>
        <Table.Tbody>
          {rows.length > 0 ? (
            rows
          ) : (
            <Table.Tr>
              <Table.Td colSpan={2}>
                <Text c="dimmed" ta="center">
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
