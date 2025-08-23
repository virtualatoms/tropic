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
import { useState, useEffect } from "react";
import SidebarFilters from "@/components/MonomerSearchSidebar";
import TabsView from "@/components/TabsView";
import DrawModal from "@/components/DrawModal";
import { MonomerFilters, MonomerSummary } from "@/lib/types";
import { fetchMonomerSummaries } from "@/lib/api";

export default function MonomerSearchPage() {
  const [modalOpen, setModalOpen] = useState(false);
  const [filters, setFilters] = useState<MonomerFilters>({
    smiles: "",
    ringSize: [1, 15],
    molWeight: [10, 500],
    hasExp: "both",
    hasComp: "both",
    functionalGroups: [],
  });

  const [monomerData, setMonomerData] = useState<MonomerSummary[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    const fetchData = async () => {
      setIsLoading(true);
      try {
        const data = await fetchMonomerSummaries(filters);
        setMonomerData(data || []); // Ensure data is always an array
      } catch (error) {
        console.error("Failed to fetch monomer summaries:", error);
        setMonomerData([]);
      } finally {
        setIsLoading(false);
      }
    };

    fetchData();
  }, [filters]);

  return (
    <>
      <Breadcrumbs pages={["Home", "Monomer Search"]} />
      <Container fluid>
        <Center my="md">
          <Group>
            <TextInput
              placeholder="Enter SMILES"
              value={filters.smiles}
              onChange={(e) => setFilters({ ...filters, smiles: e.currentTarget.value })}
              w={500}
            />
            <Button onClick={() => setModalOpen(true)}>Draw Molecule</Button>
          </Group>
        </Center>

        <DrawModal
          opened={modalOpen}
          onApply={() => {}}
          onClose={() => setModalOpen(false)}
        />

        <Grid gutter="xl" pt="xl">
          <GridCol span={3}>
            <SidebarFilters filters={filters} setFilters={setFilters} />
          </GridCol>
          <GridCol span={9}>
            <TabsView data={monomerData} isLoading={isLoading} />
          </GridCol>
        </Grid>
      </Container>
    </>
  );
}
