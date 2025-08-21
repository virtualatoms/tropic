import { Card, Badge, Stack, Box, Center, Image } from "@mantine/core";
import MoleculeImage from "./MoleculeImage/MoleculeImage";

export function MonomerLogo({
	svg,
	monomerId,
}: {
	svg: string;
	monomerId: string;
}) {
	return (
		<Box mb={20} px={0}>
			<Stack>
				<Card
					withBorder
					radius="lg"
					style={{
						borderColor: "var(--mantine-color-blue-5)",
						borderWidth: 3,
					}}
				>
					<Image src={`data:image/svg+xml;base64,${svg}`} alt={monomerId} />
					{/* <MoleculeImage smiles={smiles} /> */}
				</Card>
				<Center style={{ marginTop: -35, zIndex: 1 }}>
					<Badge size="lg" radius="lg">
						{monomerId}
					</Badge>
				</Center>
			</Stack>
		</Box>
	);
}
