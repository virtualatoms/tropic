import React, { useMemo, useRef } from "react";
import { Reaction } from "@/lib/types"; 

import { AgGridReact } from "ag-grid-react";
import {
  ModuleRegistry,
  ColDef,
  AllCommunityModule,
  themeAlpine,
} from "ag-grid-community";

import { ReactionIdCell, HtmlHeader, ImageRenderer, boolFormatter, decimalFormatter, idComparator, RefCell } from "./TableFormats";


ModuleRegistry.registerModules([AllCommunityModule]);

const theme = themeAlpine.withParams({
  headerBackgroundColor: "var(--mantine-color-gray-0)",
  borderColor: "var(--mantine-color-gray-3)",
  cellHorizontalPadding: "12px",
});


type GridRowData = {
  structure: string;
  reaction_id: string;
  is_experimental: boolean;
  type: string;
  delta_h: number | null;
  medium: string | null;
  method: string | null;
  year: number | null;
  ref: string | null; 
};


export default function ReactionSearchTable({ data }: { data: Reaction[] }) {
  const gridRef = useRef<AgGridReact>(null);

  const rowData = useMemo<GridRowData[]>(() => {
    return data.map((r) => ({
      structure: r.monomer.svg,
      reaction_id: r.reaction_id,
      type: r.type,
      delta_h: r.thermo.delta_h,
      medium: r.parameters.medium,
      method: r.parameters.method,
      is_experimental: r.parameters.is_experimental,
      year: r.metadata.year,
      ref: r.metadata.doi,
    }));
  }, [data]);

  const columnDefs = useMemo<ColDef[]>(
    () => [
      {
        headerName: "Molecule",
        field: "structure",
        cellRenderer: ImageRenderer,
        width: 180,
        autoHeight: true,
      },
      {
        headerName: "Reaction ID",
        field: "reaction_id",
        cellRenderer: ReactionIdCell,
        comparator: idComparator,
        width: 150,
      },
      { headerName: "Type", field: "type", width: 90 },
      {
        headerName: "ΔH<sub>p</sub> (kJ/mol)",
        field: "delta_h",
        headerComponent: HtmlHeader,
        valueFormatter: params => decimalFormatter(params.data.delta_h),
        width: 140,
      },
      {
        headerName: "Exp.",
        field: "is_experimental",
        valueFormatter: boolFormatter,
        width: 70,
      },
      {
        headerName: "Medium",
        field: "medium",
        width: 90,
      },
      {
        headerName: "Method",
        field: "method",
        width: 90,
      },
      { headerName: "Year", field: "year", width: 90 },
      {
        headerName: "Ref.",
        field: "ref",
        cellRenderer: RefCell,
        width: 90,
        resizable: false,
      },
    ],
    []
  );

  const defaultColDef = useMemo<ColDef>(
    () => ({
      resizable: true,
      sortable: true,
      filter: false,
      cellStyle: { display: "flex", alignItems: "center" },
    }),
    []
  );

  return (
    <AgGridReact
      ref={gridRef}
      theme={theme}
      rowData={rowData}
      columnDefs={columnDefs}
      defaultColDef={defaultColDef}
      domLayout="autoHeight"
      pagination={true}
      paginationPageSize={10}
      paginationPageSizeSelector={[10, 25, 50]}
      gridOptions={{ suppressCellFocus: true, suppressNoRowsOverlay: true }}
      animateRows={true}
      rowHeight={100}
      autoSizeStrategy={{ type: "fitGridWidth", defaultMinWidth: 10 }}
      onGridSizeChanged={(params) => params.api.sizeColumnsToFit()}
      initialState={{
        sort: { sortModel: [{ colId: "reaction_id", sort: "asc" }] },
      }}
    />
  );
}