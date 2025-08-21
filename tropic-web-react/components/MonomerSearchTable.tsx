import React, { useCallback, useMemo, useRef } from "react";
import { AgGridReactProps } from "ag-grid-react";
import { Select, Button, Anchor } from "@mantine/core";
import { AgGridReact } from "ag-grid-react";
import {
	ModuleRegistry,
	AllCommunityModule,
	themeAlpine,
} from "ag-grid-community";

ModuleRegistry.registerModules([AllCommunityModule]);

const theme = themeAlpine.withParams({
	headerBackgroundColor: "var(--mantine-color-gray-0)",
	borderColor: "var(--mantine-color-gray-3)",
  	cellHorizontalPadding: "12px",
});

const CustomHeader = ({ displayName }: { displayName: string }) => (
	<span dangerouslySetInnerHTML={{ __html: displayName }} />
);

const MonomerIdCell = ({ value }: { value: string }) => (
	<Anchor href={`/monomers/${value}`} size="sm">
		{value}
	</Anchor>
);

type Molecule = {
  structure: string;
  monomer_id: string;
  smiles: string;
  ring_size: number;
  has_exp: boolean;
  has_comp: boolean;
};

type Props = {
  getRows: (params: {
    startRow: number;
    endRow: number;
  }) => Promise<{ rows: Molecule[]; totalCount: number }>;
};

const SEARCH_NUM_ROWS = 5;

export default function MonomerSearchTable({ getRows }: Props) {
  const gridRef = useRef<AgGridReact>(null);

  const columnDefs = useMemo<AgGridReactProps["columnDefs"]>(
    () => [
      {
        headerName: "Molecule",
        field: "structure",
        cellRenderer: MarkdownRenderer,
        minWidth: 200,
      },
      {
        headerName: "Monomer ID",
        field: "monomer_id",
        cellRenderer: CustomHeader,
      },
      { headerName: "SMILES", field: "smiles", minWidth: 250 },
      { headerName: "Ring Size", field: "ring_size" },
      { headerName: "Has Exp", field: "has_exp", valueFormatter: boolFormatter },
      { headerName: "Has Comp", field: "has_comp", valueFormatter: boolFormatter },
    ],
    []
  );

  const datasource = useMemo(() => {
    return {
      getRows: async (params: any) => {
        const { startRow, endRow } = params.request;
        const res = await getRows({ startRow, endRow });

        params.successCallback(res.rows, res.totalCount);
      },
    };
  }, [getRows]);

  const onGridReady = useCallback((params: any) => {
    params.api.setDatasource(datasource);
  }, [datasource]);

  return (
    <div className="space-y-4">
      <div className="flex justify-end gap-2 items-center">
        <span className="text-sm text-muted-foreground">Export as:</span>
        {/* <Select>
          <Select.Trigger className="w-[120px]" />
          <Select.Content>
            <Select.Item value="csv">CSV</Select.Item>
            <Select.Item value="json">JSON</Select.Item>
          </Select.Content>
        </Select>
        <Button variant="outline" size="sm">Export</Button> */}
      </div>

      <div className="ag-theme-alpine" style={{ height: "500px", width: "100%" }}>
        <AgGridReact
          ref={gridRef}
          columnDefs={columnDefs}
          defaultColDef={{
            resizable: true,
            cellStyle: {
              display: "flex",
              alignItems: "center",
              height: "100%",
            },
          }}
          domLayout="autoHeight"
          rowModelType="infinite"
          cacheBlockSize={SEARCH_NUM_ROWS}
          pagination={true}
          paginationPageSize={SEARCH_NUM_ROWS}
          animateRows={false}
          onGridReady={onGridReady}
          theme={theme}
        />
      </div>
    </div>
  );
}

function MarkdownRenderer(params: any) {
  // Basic Markdown support (replace this with a markdown lib if needed)
  return (
    <div
      className="prose text-sm"
      dangerouslySetInnerHTML={{ __html: markdownToHtml(params.value) }}
    />
  );
}

function markdownToHtml(md: string): string {
  return md.replace(/\[(.*?)\]\((.*?)\)/g, `<a href="$2" target="_blank">$1</a>`);
}

function boolFormatter(params: any) {
  return params.value ? "Yes" : "No";
}