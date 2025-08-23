import React, { useMemo, useRef } from "react";
import { MonomerSummary } from "@/lib/types";

import { AgGridReact } from "ag-grid-react";
import {
  ModuleRegistry,
  ColDef,
  AllCommunityModule,
  themeAlpine,
} from "ag-grid-community";

import { MonomerIdCell, ImageRenderer, boolFormatter,  idComparator } from "./TableFormats";

ModuleRegistry.registerModules([AllCommunityModule]);

const theme = themeAlpine.withParams({
  headerBackgroundColor: "var(--mantine-color-gray-0)",
  borderColor: "var(--mantine-color-gray-3)",
  cellHorizontalPadding: "12px",
});

type GridRowData = {
  structure: string;
  monomer_id: string;
  smiles: string;
  ring_size: number;
  has_exp: boolean;
  has_comp: boolean;
};

export default function MonomerSearchTable({
  data,
}: {
  data: MonomerSummary[];
}) {
  const gridRef = useRef<AgGridReact>(null);

  const rowData = useMemo<GridRowData[]>(() => {
    return data.map((summary) => ({
      structure: summary.monomer.svg,
      monomer_id: summary.monomer.monomer_id,
      smiles: summary.monomer.smiles,
      ring_size: summary.monomer.ring_size,
      has_exp: summary.has_exp,
      has_comp: summary.has_comp,
    }));
  }, [data]);

  const columnDefs = useMemo<ColDef[]>(
    () => [
      {
        headerName: "Molecule",
        field: "structure",
        cellRenderer: ImageRenderer,
        width: 200,
        autoHeight: true,
      },
      {
        headerName: "Monomer ID",
        field: "monomer_id",
        cellRenderer: MonomerIdCell,
        comparator: idComparator,
        width: 120,
      },
      { headerName: "SMILES", field: "smiles", minWidth: 250 },
      { headerName: "Ring Size", field: "ring_size", width: 90 },
      {
        headerName: "Has Exp",
        field: "has_exp",
        valueFormatter: boolFormatter,
        width: 90,
      },
      {
        headerName: "Has Comp",
        field: "has_comp",
        valueFormatter: boolFormatter,
        width: 90,
        resizable: false,
      },
    ],
    [],
  );

  const defaultColDef = useMemo<ColDef>(
    () => ({
      resizable: true,
      sortable: true,
      filter: false,
      cellStyle: {
        display: "flex",
        alignItems: "center",
      },
    }),
    [],
  );

  return (
    <>
      <AgGridReact
        ref={gridRef}
        theme={theme}
        rowData={rowData}
        columnDefs={columnDefs}
        defaultColDef={defaultColDef}
        domLayout="autoHeight"
        pagination={true}
        paginationPageSize={8}
        paginationPageSizeSelector={[8, 20, 50, 100]}
        animateRows={true}
        gridOptions={{ suppressCellFocus: true, suppressNoRowsOverlay: true }}
        rowHeight={100}
        autoSizeStrategy={{ type: "fitGridWidth", defaultMinWidth: 10 }}
        onGridSizeChanged={(params) => params.api.sizeColumnsToFit()}
        initialState={{
          sort: { sortModel: [{ colId: "monomer_id", sort: "asc" }] },
        }}
      />
    </>
  );
}
