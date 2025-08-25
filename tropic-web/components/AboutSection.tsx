import {
  Container,
  Title,
  Text,
  Stack,
  Divider,
  Grid,
  Card,
  Image,
  ThemeIcon,
  Alert,
  Blockquote,
  Table,
  Group,
  Code,
} from "@mantine/core";
import {
  IconChartInfographic,
  IconFlask,
  IconListDetails,
  IconInfoCircle,
} from "@tabler/icons-react";
import { InlineMath, BlockMath } from "react-katex";
import "katex/dist/katex.min.css";

export default function AboutSection() {
  return (
      <Stack gap="xl">
        <Stack gap="xs">
          <Title order={1}>About the TROPIC Database</Title>
          <Text>
            TROPIC (Thermodynamics of Ring-Opening Polymerisation Informatics
            Collection) is an open-source database of experimental and
            computational thermodynamic parameters for ring-opening
            polymerization (ROP) reactions. Our goal is to facilitate the
            data-driven discovery of novel, functional polymers that are primed
            for chemical recycling.
          </Text>
        </Stack>

        <Divider />

        <Stack gap="xs">
          <Group>
            <ThemeIcon size="lg" variant="light">
              <IconListDetails />
            </ThemeIcon>
            <Title order={2} id="data-collection">
              Data Collection & Curation
            </Title>
          </Group>
          <Text>
            The data in this database has been manually collected from
            peer-reviewed academic literature. Each entry is linked to its
            source via a Digital Object Identifier (DOI), allowing for easy
            access to the original publication and its metadata.
          </Text>
          <Text>
            We include key experimental conditions, such as the reaction medium
            and monomer concentration, as these factors significantly influence
            thermodynamic values. To ensure data quality, we use a flagging
            system to highlight entries that are incomplete or appear
            inconsistent.
          </Text>
          <Alert
            icon={<IconInfoCircle size="1rem" />}
            title="Community Contributions"
            color="blue"
            mt="md"
          >
            TROPIC is an evolving resource. We encourage the community to
            contribute by submitting new or existing data. Template spreadsheets
            are available to streamline the submission process, helping us
            quickly verify and integrate new information.
          </Alert>
        </Stack>

        <Divider />

        <Stack gap="xs">
          <Title order={2} id="molecule-classes">
            Molecule Classes
          </Title>
          <Text>
            The current dataset focuses on several classes of polar cyclic
            monomers (PCMs) used in ring-opening polymerizations.
          </Text>
          <Grid gutter="md" mt="md">
            <Grid.Col span={8}>
              <Image
                radius="md"
                src="/about-classes.png"
                alt="Classes of polar cyclic monomers"
                my="md"
              />
            </Grid.Col>
            <Grid.Col span={4}>
              <Card withBorder p={0}>
                <Table striped verticalSpacing="xs">
                  <Table.Thead>
                    <Table.Tr>
                      <Table.Th>Abbreviation</Table.Th>
                      <Table.Th>Monomer Class</Table.Th>
                    </Table.Tr>
                  </Table.Thead>
                  <Table.Tbody>
                    <Table.Tr>
                      <Table.Td>L</Table.Td>
                      <Table.Td>Lactone</Table.Td>
                    </Table.Tr>
                    <Table.Tr>
                      <Table.Td>CC</Table.Td>
                      <Table.Td>Cyclic Carbonate</Table.Td>
                    </Table.Tr>
                    <Table.Tr>
                      <Table.Td>Lm</Table.Td>
                      <Table.Td>Lactam</Table.Td>
                    </Table.Tr>
                    <Table.Tr>
                      <Table.Td>tL, tnL, dtL</Table.Td>
                      <Table.Td>Sulfur-containing Lactone Derivatives</Table.Td>
                    </Table.Tr>
                    <Table.Tr>
                      <Table.Td>CtC, CdtC, etc.</Table.Td>
                      <Table.Td>
                        Sulfur-containing Cyclic Carbonate Derivatives
                      </Table.Td>
                    </Table.Tr>
                  </Table.Tbody>
                </Table>
              </Card>
            </Grid.Col>
          </Grid>
        </Stack>

        <Divider />

        <Stack gap="xs">
          <Group>
            <ThemeIcon size="lg" variant="light">
              <IconFlask />
            </ThemeIcon>
            <Title order={2} id="reaction-types">
              Reaction Types
            </Title>
          </Group>
          <Text>
            The database classifies reactions into four distinct categories to
            provide clear context for the thermodynamic data.
          </Text>
          <Stack gap="lg" mt="md">
            <Card withBorder radius="md" padding="md">
              <Grid>
                <Grid.Col span={{ base: 12, sm: 8 }}>
                  <Text fw={700}>ROP (Ring-Opening Polymerization)</Text>
                  <Text size="sm" mt={4}>
                    An experimental reaction where cyclic monomers form
                    long-chain linear polymers.
                  </Text>
                </Grid.Col>
                <Grid.Col span={{ base: 12, sm: 4 }}>
                  <Image
                    src="/rop.png"
                    radius="sm"
                    alt="Ring-Opening Polymerization"
                  />
                </Grid.Col>
              </Grid>
            </Card>

            <Card withBorder radius="md" padding="md">
              <Grid>
                <Grid.Col span={{ base: 12, sm: 8 }}>
                  <Text fw={700}>RCE (Ring-Chain Equilibrium)</Text>
                  <Text size="sm" mt={4}>
                    An experimental equilibrium between all ring sizes
                    (including monomer) and the final polymer.
                  </Text>
                </Grid.Col>
                <Grid.Col span={{ base: 12, sm: 4 }}>
                  <Image
                    src="/rce.png"
                    radius="sm"
                    alt="Ring-Chain Equilibrium"
                  />
                </Grid.Col>
              </Grid>
            </Card>

            <Card withBorder radius="md" padding="md">
              <Grid>
                <Grid.Col span={{ base: 12, sm: 8 }}>
                  <Text fw={700}>ROR (Ring-Opening Reaction)</Text>
                  <Text size="sm" mt={4}>
                    A model reaction, computational or experimental, where one
                    or a few monomers are opened to form a short linear chain.
                  </Text>
                </Grid.Col>
                <Grid.Col span={{ base: 12, sm: 4 }}>
                  <Image
                    src="/ror.png"
                    radius="sm"
                    alt="Ring-Opening Reaction"
                  />
                </Grid.Col>
              </Grid>
            </Card>

            <Card withBorder radius="md" padding="md">
              <Grid>
                <Grid.Col span={{ base: 12, sm: 8 }}>
                  <Text fw={700}>RER (Ring-Expansion Reaction)</Text>
                  <Text size="sm" mt={4}>
                    A computational model where cyclic monomers form short
                    cyclic polymer chains or &quot;loops.&quot;
                  </Text>
                </Grid.Col>
                <Grid.Col span={{ base: 12, sm: 4 }}>
                  <Image
                    src="/rer.png"
                    radius="sm"
                    alt="Ring-Expansion Reaction"
                  />
                </Grid.Col>
              </Grid>
            </Card>
          </Stack>
        </Stack>

        <Divider />

        <Stack gap="xs">
          <Group>
            <ThemeIcon size="lg" variant="light">
              <IconChartInfographic />
            </ThemeIcon>
            <Title order={2} id="thermodynamic-principles">
              Thermodynamic Principles
            </Title>
          </Group>
          <Text>
            The ceiling temperature (<InlineMath math="T_c" />) is a critical
            parameter for chemical recyclability. It represents the temperature
            at which polymerization and depolymerization are in equilibrium. It
            is derived from the enthalpy (<InlineMath math="\Delta H_p" />) and
            entropy (<InlineMath math="\Delta S_p" />) of polymerization.
          </Text>
          <Blockquote mt="xs">
            The relationship is defined by the Dainton equation, where{" "}
            <InlineMath math="[M]_{eq}" /> is the monomer concentration at
            equilibrium:
            <Text component="div" ta="center" ff="monospace" mt="sm">
              <BlockMath math="T_c = \frac{\Delta H_p}{\Delta S_p^\circ + R \cdot \ln([M]_{eq})}" />
            </Text>
          </Blockquote>
        </Stack>

        <Divider />

        <Stack gap="xs">
          <Title order={2} id="methodologies">
            Methodologies
          </Title>
          <Grid gutter="xl" mt="md">
            <Grid.Col span={{ base: 12, md: 6 }}>
              <Stack>
                <Title order={4} ta="center">
                  Experimental Methods
                </Title>
                <Card withBorder radius="md">
                  <Text fw={700}>Van&apos;t Hoff Analysis</Text>
                  <Text size="sm" mt={5}>
                    Determining <InlineMath math="\Delta H_p" /> and{" "}
                    <InlineMath math="\Delta S_p" /> from the monomer
                    equilibrium concentration at various temperatures, often
                    measured via NMR spectroscopy.
                  </Text>
                </Card>
                <Card withBorder radius="md">
                  <Text fw={700}>Differential Scanning Calorimetry (DSC)</Text>
                  <Text size="sm" mt={5}>
                    A thermal analysis technique used to measure heat flow
                    associated with polymer transitions.
                  </Text>
                </Card>
                <Card withBorder radius="md">
                  <Text fw={700}>Calorimetry</Text>
                  <Text size="sm" mt={5}>
                    Direct measurement of heat changes during polymerization.
                  </Text>
                </Card>
                <Card withBorder radius="md">
                  <Text fw={700}>NMR Spectroscopy</Text>
                  <Text size="sm" mt={5}>
                    Used to monitor monomer conversion and determine equilibrium
                    concentrations.
                  </Text>
                </Card>
              </Stack>
            </Grid.Col>

            <Grid.Col span={{ base: 12, md: 6 }}>
              <Stack>
                <Title order={4} ta="center">
                  Computational Methods
                </Title>
                <Card withBorder radius="md">
                  <Text fw={700}>Density Functional Theory (DFT)</Text>
                  <Text size="sm" mt={5}>
                    Used for Ring-Opening Reaction (ROR) models to calculate
                    thermodynamic parameters. The specific <Code>functional</Code> and
                    <Code>basis_set</Code> are recorded.
                  </Text>
                </Card>
                <Card withBorder radius="md">
                  <Text fw={700}>Molecular Dynamics (MD)</Text>
                  <Text size="sm" mt={5}>
                    Used to identify stable conformers before DFT calculations.
                    Methods include classical (<Code>ffmd</Code>), ab initio (<Code>aimd</Code>), and
                    machine-learned (<Code>mlmd</Code>) dynamics.
                  </Text>
                </Card>
                <Card withBorder radius="md">
                  <Text fw={700}>Semi-Empirical Methods</Text>
                  <Text size="sm" mt={5}>
                    Quantum mechanical methods like <Code>xTB</Code> that offer a balance
                    between speed and accuracy.
                  </Text>
                </Card>
              </Stack>
            </Grid.Col>
          </Grid>
        </Stack>
      </Stack>
  );
}
