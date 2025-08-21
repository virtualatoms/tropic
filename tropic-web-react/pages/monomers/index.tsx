import { Breadcrumbs } from "../../components/Breadcrumbs";
import {
  Container,
  TextInput,
  Button,
  Center,
  Group,
  Grid,
  GridCol,
} from "@mantine/core";
import { useState } from "react";
import SidebarFilters from "../../components/MonomerSearchSidebar";
// import TabsView from "../../components/TabsView";
import DrawModal from "../../components/DrawModal";

export default function MonomerSearchPage() {
  const [smiles, setSmiles] = useState("");
  const [modalOpen, setModalOpen] = useState(false);
  const [filters, setFilters] = useState<{
    smiles: string;
    ringSize: [number, number];
    molWeight: [number, number];
    hasExp: string;
    hasComp: string;
    functionalGroups: string[];
  }>({
    smiles: "",
    ringSize: [1, 15],
    molWeight: [10, 500],
    hasExp: "both",
    hasComp: "both",
    functionalGroups: [],
  });

  return (
    <>
      <Breadcrumbs pages={["Home", "API"]} />
      <Container fluid>
        <Center my="md">
          <Group>
            <TextInput
              placeholder="Enter SMILES"
              value={smiles}
              onChange={(e) => setSmiles(e.currentTarget.value)}
              w={500}
            />
            <Button onClick={() => setModalOpen(true)}>Draw Molecule</Button>
          </Group>
        </Center>

        <DrawModal opened={modalOpen} onApply={() => {}} onClose={() => setModalOpen(false)} />

        <Grid gutter="xl" pt="xl">
          <GridCol span={3}>
            <SidebarFilters filters={filters} setFilters={setFilters} />
          </GridCol>
          <GridCol span={9}>
            {/* <TabsView filters={filters} /> */}
          </GridCol>
        </Grid>
      </Container>
    </>
  );
}