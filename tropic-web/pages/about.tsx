import { Breadcrumbs } from "@/components/Breadcrumbs";
import { Title } from "@mantine/core";

export default function AboutPage() {
  return (
    <>
      <Breadcrumbs pages={["Home", "About"]} />
      <Title order={1}>About</Title>
    </>
  );
}
