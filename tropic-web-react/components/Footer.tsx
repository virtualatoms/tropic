import { Container, Group, Text, Box, Anchor } from "@mantine/core";
import { Icon } from "@iconify/react";
import Link from 'next/link';

export function Footer() {
	return (
		<Box bg="gray.0" mt={20}>
			<Container size="xl">
				<Group justify="space-between" py={20} px={15}>
					<Group gap="sm">
						<Icon
							icon="eos-icons:molecules-outlined"
							color="var(--mantine-color-dark-1)"
							width={24}
						/>
						<Text size="sm" c="dimmed" fw={600}>
							TROPIC
						</Text>
					</Group>
					<Group gap="sm">
						<Anchor href="/about" size="sm" c="dimmed" component={Link}>
							About
						</Anchor>
						<Anchor href="/api" size="sm" c="dimmed" component={Link}>
							API
						</Anchor>
						<Anchor
							href="https://github.com/virtualatoms/tropic"
							size="sm"
							c="dimmed"
							component={Link}
						>
							GitHub
						</Anchor>
					</Group>
				</Group>
			</Container>
		</Box>
	);
}
