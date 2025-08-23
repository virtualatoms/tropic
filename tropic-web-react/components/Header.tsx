import Link from "next/link";
import { Box, Container, Group, Text, Button, Divider } from "@mantine/core";
import { Icon } from "@iconify/react";

export function Header() {
  return (
    <Box>
      <Container size="xl" pt="md" h="60">
        <Group justify="space-between">
          <Group gap="md">
            <Icon
              icon="eos-icons:molecules-outlined"
              color="var(--mantine-color-blue-5)"
              width="30"
            />
            <Text size="xl" fw={700}>
              TROPIC
            </Text>
          </Group>
          <Group gap="sm">
            <Button
              variant="subtle"
              color="black"
              size="sm"
              fw={500}
              component={Link}
              href="/"
            >
              Home
            </Button>
            <Button
              variant="subtle"
              color="black"
              size="sm"
              fw={500}
              component={Link}
              href="/monomers"
            >
              Monomers
            </Button>
            <Button
              variant="subtle"
              color="black"
              size="sm"
              fw={500}
              component={Link}
              href="/reactions"
            >
              Reactions
            </Button>
            <Button
              variant="subtle"
              color="black"
              size="sm"
              fw={500}
              component={Link}
              href="/about"
            >
              About
            </Button>
            <Button
              variant="subtle"
              color="black"
              size="sm"
              fw={500}
              component={Link}
              href="/docs"
            >
              API
            </Button>
          </Group>
        </Group>
      </Container>
      <Divider />
    </Box>
  );
}
