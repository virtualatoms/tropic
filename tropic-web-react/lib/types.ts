export interface MonomerFilters {
	smiles: string;
	ringSize: [number, number];
	molWeight: [number, number];
	hasExp: string;
	hasComp: string;
	functionalGroups: string[];
}
