import ApiDocs from "@/components/ApiDocs";
import { Breadcrumbs } from "@/components/Breadcrumbs";
import { TableOfContents } from "@/components/TableOfContents/TableOfContents";
import { Grid } from "@mantine/core";

export default function ApiPage() {
  return (
    <>
      <Breadcrumbs pages={["Home", "API"]} />
      <Grid pt={30} gutter="xl" mb={50}>
        <Grid.Col span={3}>
          <div style={{ position: "sticky", top: "1rem" }}>
            <TableOfContents />
          </div>
        </Grid.Col>
        <Grid.Col span={9}>
          <ApiDocs />
        </Grid.Col>
      </Grid>
    </>
  );
}
