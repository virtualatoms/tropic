import { Group, Burger, Title, ActionIcon, useMantineColorScheme } from '@mantine/core';
import { IconSun, IconMoonStars } from '@tabler/icons-react';
import { Link } from 'react-router-dom';

interface HeaderProps {
  opened: boolean;
  setOpened: (opened: boolean) => void;
}

export default function Header({ opened, setOpened }: HeaderProps) {
  const { colorScheme, toggleColorScheme } = useMantineColorScheme();

  return (
    <Group h="100%" px="md" justify="space-between">
      <Group>
        <Burger opened={opened} onClick={() => setOpened(!opened)} size="sm" />
        <Link to="/" style={{ textDecoration: 'none', color: 'inherit' }}>
          <Title order={3}>Roppy</Title>
        </Link>
      </Group>
      
      <Group>
        <ActionIcon
          variant="default"
          onClick={() => toggleColorScheme()}
          size="lg"
          aria-label="Toggle color scheme"
        >
          {colorScheme === 'dark' ? <IconSun size="1.2rem" /> : <IconMoonStars size="1.2rem" />}
        </ActionIcon>
      </Group>
    </Group>
  );
} 