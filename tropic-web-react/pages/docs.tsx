import { Breadcrumbs } from "../components/Breadcrumbs";
import { Title } from "@mantine/core";

export default function ApiPage() {
  return (
    <>
      <Breadcrumbs pages={["Home", "API"]} />
      <Title order={1}>API</Title>
    </>
  );
}
