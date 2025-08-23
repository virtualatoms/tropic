import { Breadcrumbs } from "@/components/Breadcrumbs";
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
import ReactionSearchSidebar from "@/components/ReactionSearchSidebar";
import ReactionTabsView from "@/components/ReactionTabsView";
import DrawModal from "@/components/DrawModal";
import { ReactionFilters, Reaction } from "@/lib/types";
import { fetchReactions } from "@/lib/api";

export default function ReactionSearchPage() {
  const [modalOpen, setModalOpen] = useState(false);
  const currentYear = new Date().getFullYear();

  const [filters, setFilters] = useState<ReactionFilters>({
    smiles: "",
    ringSize: [1, 15],
    molecularWeight: [10, 500],
    isExperimental: "both",
    initiator: "",
    medium: [],
    reactionType: [],
    method: [],
    year: [1950, currentYear],
  });

  const [monomerData, setReactionData] = useState<Reaction[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    const fetchData = async () => {
      setIsLoading(true);
      try {
        const data = await fetchReactions(filters);
        setReactionData(data || []); // Ensure data is always an array
      } catch (error) {
        console.error("Failed to fetch reaction summaries:", error);
        setReactionData([]);
      } finally {
        setIsLoading(false);
      }
    };

    fetchData();
  }, [filters]);

  return (
    <>
      <Breadcrumbs pages={["Home", "Reaction Search"]} />
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
            <ReactionSearchSidebar filters={filters} setFilters={setFilters} />
          </GridCol>
          <GridCol span={9}>
            <ReactionTabsView data={monomerData} isLoading={isLoading} />
          </GridCol>
        </Grid>
      </Container>
    </>
  );
}
