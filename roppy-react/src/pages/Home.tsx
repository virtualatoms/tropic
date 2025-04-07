import { Title, Text, Container, Stack } from '@mantine/core';

export default function Home() {
  return (
    <Stack>
      <Title order={1}>Welcome to Roppy</Title>
      <Text size="lg">
        Roppy is a web application for exploring and analyzing polymer structures.
        Use the sidebar to navigate to different sections of the application.
      </Text>
      <Container>
        <Title order={2}>Features</Title>
        <Stack>
          <Text>
            • Browse and search through a database of monomers
          </Text>
          <Text>
            • View detailed information about each monomer
          </Text>
          <Text>
            • Access the API documentation for integration
          </Text>
        </Stack>
      </Container>
    </Stack>
  );
} 