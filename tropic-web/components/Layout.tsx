import { ReactNode } from "react";
import { Container, Box } from "@mantine/core";
import { Header } from "./Header";
import { Footer } from "./Footer";

export function Layout({ children }: { children: ReactNode }) {
  return (
    <Box
      style={{ display: "flex", flexDirection: "column", minHeight: "100vh" }}
    >
      <Box style={{ flex: 1 }}>
        <Container size="xl">
          <Header />
          {children}
        </Container>
      </Box>
      <Footer />
    </Box>
  );
}
