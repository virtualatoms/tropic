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
    svg: string;
  };
  data: ReactionData[];
  has_exp: boolean;
  has_comp: boolean;
}

export interface Reaction {
  monomer: {
    monomer_id: string;
    smiles: string;
    ring_size: number;
  };
  type: string;
  parameters: {
    is_experimental: boolean;
    state_summary: string;
    initial_monomer_conc: number | null;
    bulk_monomer_conc: number | null;
    medium: string | null;
    method: string | null;
  };
  thermo: {
    delta_h: number | null;
    delta_s: number | null;
    ceiling_temperature: number | null;
  };
  product: {
    repeating_units: number | null;
  };
  metadata: {
    year: number;
    formatted_reference: string;
  };
}

export interface MonomerSummaryState {
  data: MonomerSummary[];
  isLoading: boolean;
}
