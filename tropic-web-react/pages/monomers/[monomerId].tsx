import { useRouter } from "next/router";
import { useEffect, useState } from "react";
import { Grid, GridCol, Loader, Center, Alert } from "@mantine/core";
import { IconAlertCircle } from "@tabler/icons-react";
import { Breadcrumbs } from "@/components/Breadcrumbs";
import { MonomerSummary } from "@/components/MonomerSummary";
import { MonomerLogo } from "@/components/MonomerLogo";
import { TableOfContents } from "@/components/TableOfContents/TableOfContents";
import MonomerReactionTable from "@/components/MonomerReactionTable";
import { API_ENDPOINT } from "@/lib/constants";

export default function MonomerPage() {
  const router = useRouter();
  const { monomerId } = router.query;

  const [data, setData] = useState<any>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!monomerId) return;

    setLoading(true);
    setError(null);

    fetch(`${API_ENDPOINT}/monomer-summaries/${monomerId}`)
      .then((res) => {
        if (!res.ok) throw new Error("Failed to fetch data");
        return res.json();
      })
      .then((json) => {
        setData(json);
      })
      .catch(() => {
        setError("Monomer not found");
        setData(null);
      })
      .finally(() => setLoading(false));
  }, [monomerId]);

  return (
    <>
      <link
        href="https://fonts.googleapis.com/css?family=Droid+Sans:400,700"
        rel="stylesheet"
      />
      <Breadcrumbs
        pages={["Home", "Monomer Search", String(monomerId || "")]}
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
                svg={data.monomer.svg}
                monomerId={data.monomer.monomer_id}
              />
              <TableOfContents />
            </div>
          </Grid.Col>
          <GridCol span={9}>
            <MonomerSummary data={data} />
            <MonomerReactionTable data={data} />
          </GridCol>
        </Grid>
      )}
    </>
  );
}
