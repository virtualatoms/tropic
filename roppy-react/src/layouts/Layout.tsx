import { AppShell, Container } from '@mantine/core';
import { useState } from 'react';
import Header from '../components/Header';
import Sidebar from '../components/Sidebar';

interface LayoutProps {
  children: React.ReactNode;
}

export default function Layout({ children }: LayoutProps) {
  const [opened, setOpened] = useState(false);

  return (
    <AppShell
      header={<AppShell.Header height={60}><Header opened={opened} setOpened={setOpened} /></AppShell.Header>}
      navbar={<AppShell.Navbar width={{ base: 300 }}><Sidebar /></AppShell.Navbar>}
      padding="md"
    >
      <Container size="xl">
        {children}
      </Container>
    </AppShell>
  );
} 