import Link from "next/link";
import { Icon } from "@iconify/react";
import {
  Container,
  Group,
  Stack,
  Burger,
  Box,
  Text,
  Button,
  Divider,
  Paper,
  Collapse,
} from "@mantine/core";
import { useDisclosure } from "@mantine/hooks";

const links = [
  { link: "/", label: "Home" },
  { link: "/monomers", label: "Monomers" },
  { link: "/reactions", label: "Reactions" },
  { link: "/about", label: "About" },
  { link: "/docs", label: "API" },
];

export function Header() {
  const [opened, { toggle }] = useDisclosure(false);

  const items = links.map((link) => (
    <Button
      key={link.label}
      component={Link}
      href={link.link}
      variant="subtle"
      color="black"
      size={opened ? "lg" : "sm"}
      fw={500}
      fullWidth={opened}
    >
      {link.label}
    </Button>
  ));

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
          <Group gap="sm" visibleFrom="sm">
            {items}
          </Group>
          <Burger opened={opened} onClick={toggle} hiddenFrom="sm" size="sm" />
        </Group>
      </Container>
      <Collapse in={opened}>
        <Paper shadow="md">
          <Stack align="center" gap={0} py="sm">
            {items}
          </Stack>
        </Paper>
      </Collapse>
      <Divider />
    </Box>
  );
}

// export function Header() {

//   return (
//     <Box>
//       <header style={{ height: 60 }}>
//         <Container size="xl" style={{ height: '100%' }}>
//           <Group justify="space-between" style={{ height: '100%' }}>
//             <Group gap="md">
//               <Icon
//                 icon="eos-icons:molecules-outlined"
//                 color="var(--mantine-color-blue-5)"
//                 width="30"
//               />
//               <Text size="xl" fw={700}>
//                 TROPIC
//               </Text>
//             </Group>

//             {/* Desktop Navigation */}
//             <Group gap="sm" visibleFrom="sm">
//               {items}
//             </Group>

//             {/* Mobile Burger Menu */}
//             <Burger opened={opened} onClick={toggle} hiddenFrom="sm" size="sm" />
//           </Group>
//         </Container>
//       </header>

//       {/* Mobile Navigation Dropdown */}
//       <Collapse in={opened}>
//         <Paper shadow="md">
//           <Stack align="center" gap={0} py="sm">
//             {items}
//           </Stack>
//         </Paper>
//       </Collapse>
//       <Divider />
//     </Box>
//   );
// }
// export function Header() {
//   const [opened, { toggle }] = useDisclosure(false);

//   const items = links.map((link) => (
//     <Button
//       key={link.label}
//       component={Link}
//       href={link.link}
//       variant="subtle"
//       color="black"
//       size="sm"
//       fw={500}
//     >
//       {link.label}
//     </Button>
//   ));

//   return (
//     <Box>
//       <Container size="xl" pt="md" h="60">
//         <Group justify="space-between">
//           <Group gap="md">
//             <Icon
//               icon="eos-icons:molecules-outlined"
//               color="var(--mantine-color-blue-5)"
//               width="30"
//             />
//             <Text size="xl" fw={700}>
//               TROPIC
//             </Text>
//           </Group>
//           <Group gap="sm">
// <Button
//   variant="subtle"
//   color="black"
//   size="sm"
//   fw={500}
//   component={Link}
//   href="/monomers"
// >
//               Home
//             </Button>
//             <Button
//               variant="subtle"
//               color="black"
//               size="sm"
//               fw={500}
//               component={Link}
//               href="/monomers"
//             >
//               Monomers
//             </Button>
//             <Button
//               variant="subtle"
//               color="black"
//               size="sm"
//               fw={500}
//               component={Link}
//               href="/reactions"
//             >
//               Reactions
//             </Button>
//             <Button
//               variant="subtle"
//               color="black"
//               size="sm"
//               fw={500}
//               component={Link}
//               href="/about"
//             >
//               About
//             </Button>
//             <Button
//               variant="subtle"
//               color="black"
//               size="sm"
//               fw={500}
//               component={Link}
//               href="/docs"
//             >
//               API
//             </Button>
//           </Group>
//         </Group>
//       </Container>
//       <Divider />
//     </Box>
//   );
// }
