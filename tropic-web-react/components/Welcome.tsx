import {
  Container,
  Title,
  Text,
  Button,
  Stack,
  ThemeIcon,
} from "@mantine/core";
import { IconFlask } from "@tabler/icons-react";
import Link from "next/link";

export function Welcome() {
  return (
    <Container size="md" style={{ paddingTop: "80px", paddingBottom: "80px" }}>
      <Stack align="center" gap="xl">
        <ThemeIcon
          size={80}
          radius={80}
          variant="gradient"
          gradient={{ from: "blue", to: "cyan" }}
        >
          <IconFlask style={{ width: "50px", height: "50px" }} />
        </ThemeIcon>

        <Title
          order={1}
          ta="center"
          style={{
            fontSize: "clamp(2rem, 5vw, 3.5rem)",
            fontWeight: 900,
          }}
        >
          Welcome to the TROPIC Polymerization Thermodynamics Database
        </Title>

        <Text size="lg" c="dimmed" ta="center" maw={600}>
          An open-source repository of experimental and computational
          thermodynamic data for polymerization reactions. Explore, search, and
          analyze a curated collection of monomer and reaction properties to
          accelerate your research.
        </Text>

        <Button
          component={Link}
          href="/monomers"
          size="lg"
          variant="gradient"
          gradient={{ from: "blue", to: "cyan" }}
          radius="xl"
          mt="md"
        >
          Start Exploring Now
        </Button>
      </Stack>
    </Container>
  );
}
