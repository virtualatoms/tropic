import { Title, Text, Stack, Paper } from '@mantine/core';

export default function About() {
  return (
    <Stack>
      <Title order={1}>About Roppy</Title>
      
      <Paper p="md" withBorder>
        <Title order={2}>What is Roppy?</Title>
        <Text>
          Roppy is a web application designed to help researchers and scientists explore
          and analyze polymer structures. It provides a user-friendly interface for
          browsing monomer databases, visualizing molecular structures, and accessing
          polymer-related data.
        </Text>
      </Paper>

      <Paper p="md" withBorder>
        <Title order={2}>Features</Title>
        <Text>
          • Comprehensive monomer database with detailed information
        </Text>
        <Text>
          • Interactive 3D molecular visualization
        </Text>
        <Text>
          • Search and filter capabilities
        </Text>
        <Text>
          • RESTful API for programmatic access
        </Text>
      </Paper>

      <Paper p="md" withBorder>
        <Title order={2}>Technology Stack</Title>
        <Text>
          • Frontend: React with TypeScript and Mantine UI
        </Text>
        <Text>
          • 3D Visualization: 3Dmol.js
        </Text>
        <Text>
          • Backend: Python FastAPI
        </Text>
      </Paper>
    </Stack>
  );
} 