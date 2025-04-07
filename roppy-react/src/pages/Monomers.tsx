import { useState, useEffect } from 'react';
import { Title, TextInput, Stack, Group } from '@mantine/core';
import { AgGridReact } from 'ag-grid-react';
import 'ag-grid-community/styles/ag-grid.css';
import 'ag-grid-community/styles/ag-theme-alpine.css';
import { IconSearch } from '@tabler/icons-react';

interface Monomer {
  id: string;
  name: string;
  formula: string;
  molecular_weight: number;
  ring_size: number;
  has_exp: boolean;
  has_comp: boolean;
}

export default function Monomers() {
  const [monomers, setMonomers] = useState<Monomer[]>([]);
  const [search, setSearch] = useState('');

  useEffect(() => {
    // TODO: Replace with actual API call
    const fetchMonomers = async () => {
      try {
        const response = await fetch('/api/monomers');
        const data = await response.json();
        setMonomers(data);
      } catch (error) {
        console.error('Error fetching monomers:', error);
      }
    };

    fetchMonomers();
  }, []);

  const filteredMonomers = monomers.filter((monomer) =>
    monomer.name.toLowerCase().includes(search.toLowerCase())
  );

  const columnDefs = [
    { headerName: 'Structure', field: 'structure' },
    { headerName: 'Monomer ID', field: 'id' },
    { headerName: 'SMILES', field: 'formula' },
    { headerName: 'Ring Size', field: 'ring_size' },
    { headerName: 'Has Exp', field: 'has_exp' },
    { headerName: 'Has Comp', field: 'has_comp' },
  ];

  return (
    <Stack>
      <Group justify="space-between">
        <Title order={1}>Monomers</Title>
        <TextInput
          placeholder="Search monomers..."
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          leftSection={<IconSearch size="1rem" />}
        />
      </Group>

      <div className="ag-theme-alpine" style={{ height: 400, width: '100%' }}>
        <AgGridReact
          rowData={filteredMonomers}
          columnDefs={columnDefs}
          pagination={true}
          paginationPageSize={5}
        />
      </div>
    </Stack>
  );
} 