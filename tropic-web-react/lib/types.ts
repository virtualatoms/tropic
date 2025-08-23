export interface MonomerFilters {
	smiles: string;
	ringSize: [number, number];
	molWeight: [number, number];
	hasExp: string;
	hasComp: string;
	functionalGroups: string[];
}


export interface ReactionData {
	is_experimental: boolean;
	delta_h: number | null;
	delta_s: number | null;
	ceiling_temperature: number | null;
	doi: string | null;
	formatted_reference: string | null;
}

export interface MonomerSummary {
	monomer: {
		monomer_id: string;
		smiles: string;
		ring_size: number;
	};
	data: ReactionData[];
	has_exp: boolean;
	has_comp: boolean;
}