import { NavLink, Stack } from '@mantine/core';
import { IconHome, IconApi, IconInfoCircle, IconApps } from '@tabler/icons-react';
import { Link, useLocation } from 'react-router-dom';

export default function Sidebar() {
  const location = useLocation();

  const links = [
    { icon: IconHome, label: 'Home', to: '/' },
    { icon: IconApps, label: 'Monomers', to: '/monomers' },
    { icon: IconApi, label: 'API Documentation', to: '/api' },
    { icon: IconInfoCircle, label: 'About', to: '/about' },
  ];

  return (
    <Stack>
      {links.map((link) => (
        <NavLink
          key={link.to}
          component={Link}
          to={link.to}
          label={link.label}
          leftSection={<link.icon size="1.2rem" stroke={1.5} />}
          active={location.pathname === link.to}
          variant="filled"
        />
      ))}
    </Stack>
  );
} 