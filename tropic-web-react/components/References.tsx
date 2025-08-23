import React, { useMemo } from "react";
import { Table, Anchor, Text, LoadingOverlay } from "@mantine/core";
import { MonomerSummary, Reaction } from "@/lib/types";

interface Reference {
  doi: string;
  formatted: string;
}

type ReferencesProps<T> = {
  data: T[];
  isLoading: boolean;
  onProcessData: (data: T[]) => Reference[];
};

export const processRefMonomerSummary = (
  data: MonomerSummary[],
): Reference[] => {
  const uniqueRefs = new Map<string, string>();
  for (const summary of data) {
    for (const result of summary.data || []) {
      const { doi, formatted_reference: formattedRef } = result;
      if (doi && formattedRef && !uniqueRefs.has(doi)) {
        uniqueRefs.set(doi, formattedRef);
      }
    }
  }
  return Array.from(uniqueRefs, ([doi, formatted]) => ({ doi, formatted }));
};

export const processRefReaction = (data: Reaction[]): Reference[] => {
  const uniqueRefs = new Map<string, string>();
  for (const reaction of data) {
    const { doi, formatted_reference: formattedRef } = reaction.metadata;
    if (doi && formattedRef && !uniqueRefs.has(doi)) {
      uniqueRefs.set(doi, formattedRef);
    }
  }
  return Array.from(uniqueRefs, ([doi, formatted]) => ({ doi, formatted }));
};

export function References<T>({
  data,
  isLoading,
  onProcessData,
}: ReferencesProps<T>) {
  const references = useMemo<Reference[]>(() => {
    return onProcessData(data);
  }, [data, onProcessData]);

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
