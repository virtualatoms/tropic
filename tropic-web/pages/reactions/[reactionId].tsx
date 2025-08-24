import { useRouter } from "next/router";
import { useEffect, useState } from "react";
import { Grid, GridCol, Loader, Center, Alert } from "@mantine/core";
import { IconAlertCircle } from "@tabler/icons-react";
import { Breadcrumbs } from "@/components/Breadcrumbs";
import { MonomerLogo } from "@/components/MonomerLogo";
import { TableOfContents } from "@/components/TableOfContents/TableOfContents";
import { API_ENDPOINT } from "@/lib/constants";
import { Reaction } from "@/lib/types";
import ReactionSummary from "@/components/ReactionSummary";

export default function MonomerPage() {
  const router = useRouter();
  const { reactionId } = router.query;

  const [data, setData] = useState<Reaction | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!reactionId) return;

    setLoading(true);
    setError(null);

    fetch(`${API_ENDPOINT}/reactions/${reactionId}?include_svg=true`)
      .then((res) => {
        if (!res.ok) throw new Error("Failed to fetch data");
        return res.json();
      })
      .then((json) => {
        setData(json);
      })
      .catch(() => {
        setError("Reaction not found");
        setData(null);
      })
      .finally(() => setLoading(false));
  }, [reactionId]);

  return (
    <>
      <link
        href="https://fonts.googleapis.com/css?family=Droid+Sans:400,700"
        rel="stylesheet"
      />
      <Breadcrumbs
        pages={["Home", "Reaction Search", String(reactionId || "")]}
      />

      {loading ? (
        <Center mt="xl">
          <Loader size="lg" />
        </Center>
      ) : error ? (
        <Center mt="xl">
          <Alert
            icon={<IconAlertCircle size="1rem" />}
            title="Error"
            color="red"
          >
            {error}
          </Alert>
        </Center>
      ) : (
        <Grid pt={30} gutter="xl" mb={50}>
          <Grid.Col span={3}>
            <div style={{ position: "sticky", top: "1rem" }}>
              <MonomerLogo
                svg={data?.monomer.svg || ""}
                monomerId={data?.monomer.monomer_id || ""}
              />
              <TableOfContents />
            </div>
          </Grid.Col>
          <GridCol span={9}>
            <ReactionSummary data={data} />
          </GridCol>
        </Grid>
      )}
    </>
  );
}
