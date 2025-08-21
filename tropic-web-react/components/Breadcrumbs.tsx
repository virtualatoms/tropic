import {
	Box,
	Breadcrumbs as MantineBreadcrumbs,
	Anchor,
	Text,
} from "@mantine/core";
import Link from "next/link";

const LINKS: Record<string, string> = {
	Home: "/",
	"Monomer Search": "/monomers",
	"Reaction Search": "/reactions",
	About: "/about",
	API: "/docs",
};

export function Breadcrumbs({ pages }: { pages: string[] }) {
	const items = pages.map((page, index) => {
		if (!(index === pages.length - 1) && LINKS[page]) {
			return (
				<Anchor component={Link} href={LINKS[page]} fw={700} key={page}>
					{page}
				</Anchor>
			);
		} else {
			return (
				<Text key={page} fw={700}>
					{page}
				</Text>
			);
		}
	});

	return (
		<Box py={10} px={20} mt={10} bg="var(--mantine-color-gray-0)">
			<MantineBreadcrumbs>{items}</MantineBreadcrumbs>
		</Box>
	);
}
