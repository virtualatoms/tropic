import { Grid, GridCol, Title } from "@mantine/core";
import { MonomerInfoCard } from "./MonomerInfoCard";
import { MoleculeViewerCard } from "./MoleculeViewerCard";
import classes from "./TableOfContents/TableOfContents.module.css";

export function MonomerSummary({ data }: { data: any }) {
	return (
		<Grid gutter="xl" pb={10}>
			<Title order={1} id="summary" className={classes.visuallyHidden}>
				Summary
			</Title>
			<GridCol span={6}>
				<MoleculeViewerCard xyz={data.monomer.xyz} />
			</GridCol>
			<GridCol span={6}>
				<MonomerInfoCard data={data} />
			</GridCol>
		</Grid>
	);
}
