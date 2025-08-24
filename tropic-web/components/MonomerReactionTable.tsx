import React, { useMemo } from "react";
import { Stack, Title } from "@mantine/core";
import { AgGridReact } from "ag-grid-react";
import {
  ModuleRegistry,
  ColDef,
  AllCommunityModule,
  themeAlpine,
} from "ag-grid-community";
import {
  ReactionIdCell,
  HtmlHeader,
  RefCell,
  decimalFormatter,
} from "./TableFormats";

ModuleRegistry.registerModules([AllCommunityModule]);

const theme = themeAlpine.withParams({
  headerBackgroundColor: "var(--mantine-color-gray-0)",
  borderColor: "var(--mantine-color-gray-3)",
  cellHorizontalPadding: "12px",
});

type GridRowData = {
  reaction_id: string;
  type: string;
  delta_h: number | null;
  delta_s: number | null;
  ceiling_temperature: number | null;
  state: string | null;
  method: string | null;
  year: number | null;
  ref: string | null;
  repeating_units: number | null;
};

const MonomerReactionTable = ({ data }: { data: any }) => {
  const COMMON_COLUMNS = useMemo<ColDef[]>(
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
        headerComponent: HtmlHeader,
        valueFormatter: (params) => decimalFormatter(params.data.delta_h),
      },
      {
        headerName: "ΔS<sub>p</sub> (J/K/mol)",
        field: "delta_s",
        headerComponent: HtmlHeader,
        valueFormatter: (params) => decimalFormatter(params.data.delta_s),
      },
      {
        headerName: "T<sub>c</sub> (K)",
        field: "ceiling_temperature",
        headerComponent: HtmlHeader,
        valueFormatter: (params) =>
          decimalFormatter(params.data.ceiling_temperature),
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

  const expTableData: GridRowData[] = [];
  const compTableData: GridRowData[] = [];

  data.data.forEach((row) => {
    if (row.is_experimental) {
      expTableData.push({
        reaction_id: row.reaction_id,
        type: row.type,
        delta_h: row.delta_h,
        delta_s: row.delta_s,
        ceiling_temperature: row.ceiling_temperature,
        year: row.year,
        ref: row.doi,
        state: row.state_summary,
        method: null,
        repeating_units: null,
      });
    } else {
      compTableData.push({
        reaction_id: row.reaction_id,
        type: row.type,
        delta_h: row.delta_h,
        delta_s: row.delta_s,
        ceiling_temperature: row.ceiling_temperature,
        year: row.year,
        ref: row.doi,
        method: row.method,
        repeating_units: row.repeating_units,
        state: null,
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
