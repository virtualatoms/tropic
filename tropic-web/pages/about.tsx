import AboutSection from "@/components/AboutSection";
import { Breadcrumbs } from "@/components/Breadcrumbs";
import { TableOfContents } from "@/components/TableOfContents/TableOfContents";
import { Grid } from "@mantine/core";

export default function AboutPage() {
  return (
    <>
      <Breadcrumbs pages={["Home", "About"]} />
      <Grid pt={30} gutter="xl" mb={50}>
        <Grid.Col span={3} pt={25}>
          <div style={{ position: "sticky", top: "1rem" }}>
            <TableOfContents />
          </div>
        </Grid.Col>
        <Grid.Col span={9}>
          <AboutSection />
        </Grid.Col>
      </Grid>
    </>
  );
}
