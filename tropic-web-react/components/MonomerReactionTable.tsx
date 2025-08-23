import React, { useMemo } from "react";
import { Stack, Title, Anchor } from "@mantine/core";
import { AgGridReact } from "ag-grid-react";
import {
  ModuleRegistry,
  AllCommunityModule,
  themeAlpine,
} from "ag-grid-community";
import Link from "next/link";

ModuleRegistry.registerModules([AllCommunityModule]);
// ModuleRegistry.registerModules([ClientSideRowModelModule]);

const theme = themeAlpine.withParams({
  headerBackgroundColor: "var(--mantine-color-gray-0)",
  borderColor: "var(--mantine-color-gray-3)",
  cellHorizontalPadding: "12px",
});

const CustomHeader = ({ displayName }: { displayName: string }) => (
  <span dangerouslySetInnerHTML={{ __html: displayName }} />
);

const ReactionIdCell = ({ value }: { value: string }) => (
  <Anchor href={`/reactions/${value}`} size="sm" component={Link}>
    {value}
  </Anchor>
);

const RefCell = ({ value }: { value: string }) => (
  <Anchor href={`https://doi.org/${value}`} target="_blank">
    <svg
      xmlns="http://www.w3.org/2000/svg"
      width="14"
      height="14"
      viewBox="0 0 24 24"
      style={{ display: "inline" }}
    >
      <path
        fill="none"
        stroke="currentColor"
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth="3"
        d="M13.5 10.5L21 3m-5 0h5v5m0 6v5a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5"
      />
    </svg>
  </Anchor>
);

const MonomerReactionTable = ({ data }: { data: any }) => {
  const COMMON_COLUMNS = useMemo(
    () => [
      {
        headerName: "Reaction ID",
        field: "reaction_id",
        cellRenderer: ReactionIdCell,
        maxWidth: 150,
      },
      {
        headerName: "Type",
        field: "type",
        maxWidth: 100,
      },
      {
        headerName: "ΔH<sub>p</sub> (kJ/mol)",
        field: "delta_h",
        headerComponent: CustomHeader,
      },
      {
        headerName: "ΔS<sub>p</sub> (J/K/mol)",
        field: "delta_s",
        headerComponent: CustomHeader,
      },
      {
        headerName: "T<sub>c</sub> (K)",
        field: "ceiling_temperature",
        headerComponent: CustomHeader,
      },
      {
        headerName: "Year",
        field: "year",
        maxWidth: 120,
      },
      {
        headerName: "Ref.",
        field: "ref",
        cellRenderer: RefCell,
        maxWidth: 80,
        resizable: false,
      },
    ],
    [],
  );

  const autoSizeStrategy = useMemo(() => {
    return {
      type: "fitGridWidth",
      defaultMinWidth: 10,
    };
  }, []);

  const EXP_COLUMNS = useMemo(
    () => [
      ...COMMON_COLUMNS.slice(0, 2),
      { headerName: "State", field: "state", width: 90, maxWidth: 90 },
      ...COMMON_COLUMNS.slice(2),
    ],
    [COMMON_COLUMNS],
  );

  const COMP_COLUMNS = useMemo(
    () => [
      ...COMMON_COLUMNS.slice(0, 2),
      { headerName: "Method", field: "method", width: 120, maxWidth: 120 },
      {
        headerName: "# Units",
        field: "repeating_units",
        width: 110,
        maxWidth: 110,
      },
      ...COMMON_COLUMNS.slice(2),
    ],
    [COMMON_COLUMNS],
  );

  const expTableData = [];
  const compTableData = [];

  data.data.forEach((row) => {
    let ceiling_temperature = "";
    if (row.ceiling_temperature) {
      ceiling_temperature = `${row.ceiling_temperature.toFixed(2)}`;
    } else if (
      row.delta_h &&
      row.delta_s &&
      row.delta_s < 0 &&
      row.delta_h < 0
    ) {
      ceiling_temperature = `${((row.delta_h * 1000) / row.delta_s).toFixed(2)}`;
    }

    if (row.is_experimental) {
      expTableData.push({
        reaction_id: row.reaction_id,
        type: row.type,
        delta_h: row.delta_h !== null ? row.delta_h.toFixed(2) : "",
        delta_s: row.delta_s !== null ? row.delta_s.toFixed(2) : "",
        ceiling_temperature,
        year: row.year,
        ref: row.doi,
        state: row.state_summary,
      });
    } else {
      compTableData.push({
        reaction_id: row.reaction_id,
        type: row.type,
        delta_h: row.delta_h !== null ? row.delta_h.toFixed(2) : "",
        delta_s: row.delta_s !== null ? row.delta_s.toFixed(2) : "",
        ceiling_temperature,
        year: row.year,
        ref: row.doi,
        method: row.method,
        repeating_units: row.repeating_units,
      });
    }
  });

  return (
    <Stack gap="md" mt={10}>
      {expTableData.length > 0 && (
        <Stack gap={10}>
          <Title order={1} id="exp">
            Experimental Data
          </Title>
          <AgGridReact
            columnDefs={EXP_COLUMNS}
            theme={theme}
            rowData={expTableData}
            domLayout="autoHeight"
            suppressCellFocus
            autoSizeStrategy={autoSizeStrategy}
            onGridSizeChanged={(params) => params.api.sizeColumnsToFit()}
          />
        </Stack>
      )}
      {compTableData.length > 0 && (
        <Stack gap={10}>
          <Title order={1} id="comp">
            Computational Data
          </Title>
          <AgGridReact
            theme={theme}
            columnDefs={COMP_COLUMNS}
            rowData={compTableData}
            domLayout="autoHeight"
            suppressCellFocus
            autoSizeStrategy={autoSizeStrategy}
            onGridSizeChanged={(params) => params.api.sizeColumnsToFit()}
          />
        </Stack>
      )}
    </Stack>
  );
};

export default MonomerReactionTable;
