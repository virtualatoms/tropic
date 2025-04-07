import { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import { Title, Stack, Group, Text, Paper } from '@mantine/core';
import { create3DmolViewer } from '../components/MolViewer';

interface Monomer {
  id: string;
  name: string;
  formula: string;
  molecular_weight: number;
  smiles: string;
  description: string;
}

export default function MonomerDetail() {
  const { id } = useParams();
  const [monomer, setMonomer] = useState<Monomer | null>(null);
  const [viewer, setViewer] = useState<any>(null);

  useEffect(() => {
    const fetchMonomer = async () => {
      try {
        const response = await fetch(`/api/monomers/${id}`);
        const data = await response.json();
        setMonomer(data);
      } catch (error) {
        console.error('Error fetching monomer:', error);
      }
    };

    fetchMonomer();
  }, [id]);

  useEffect(() => {
    if (monomer && !viewer) {
      const container = document.getElementById('molviewer');
      if (container) {
        const newViewer = create3DmolViewer(container);
        newViewer.addModel(monomer.smiles, 'sdf');
        newViewer.setStyle({ stick: {} });
        newViewer.zoomTo();
        newViewer.render();
        setViewer(newViewer);
      }
    }
  }, [monomer, viewer]);

  if (!monomer) {
    return <Text>Loading...</Text>;
  }

  return (
    <Stack>
      <Title order={1}>{monomer.name}</Title>
      
      <Group align="start" grow>
        <Stack>
          <Paper p="md" withBorder>
            <Title order={3}>Properties</Title>
            <Text>Formula: {monomer.formula}</Text>
            <Text>Molecular Weight: {monomer.molecular_weight.toFixed(2)}</Text>
          </Paper>

          <Paper p="md" withBorder>
            <Title order={3}>Description</Title>
            <Text>{monomer.description}</Text>
          </Paper>
        </Stack>

        <Paper p="md" withBorder>
          <Title order={3}>3D Structure</Title>
          <div id="molviewer" style={{ width: '100%', height: '400px' }} />
        </Paper>
      </Group>
    </Stack>
  );
} 