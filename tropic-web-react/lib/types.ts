export interface MonomerFilters {
  smiles: string;
  ringSize: [number, number];
  molWeight: [number, number];
  hasExp: string;
  hasComp: string;
  functionalGroups: string[];
}

export interface ReactionFilters {
  smiles: string;
  ringSize: [number, number];
  molecularWeight: [number, number];
  isExperimental: 'yes' | 'no' | 'both';
  initiator: string;
  medium: string[];
  reactionType: string[];
  method: string[];
  year: [number, number];
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

export interface Monomer {
  monomer_id: string;
  smiles: string; 
  inchi: string | null;
  molecular_weight: number | null;
  functional_groups: string[] | null;
  iupac_name: string | null;
  pubchem_cid: number | null;
  ring_size: number | null;
  xyz: string | null;
  svg: string;
}

export interface Product {
  smiles: string | null;
  repeating_units: number | null;
  deg_of_poly: number | null;
  dispersity: number | null;
  n_avg_molar_mass: number | null;
  m_avg_molar_mass: number | null;
}

export interface Parameters {
  is_experimental: boolean;
  temperature: number | null;
  pressure: number | null;
  monomer_state: string | null; 
  polymer_state: string | null;
  initiator_smiles: string | null;
  initial_monomer_conc: number | null;
  bulk_monomer_conc: number | null;
  medium: string | null;
  solvent: string | null;
  cosolvent: string | null;
  solvent_cosolvent_ratio: string | null;
  topology: string | null;
  method: string | null;
  method_calc: string | null;
  functional: string | null;
  basis_set: string | null;
  dispersion: string | null;
  forcefield: string | null;
  solvent_model: string | null;
  state_summary: string;
}

export interface Vanthoff {
  temperature: number[] | null;
  inverse_temperature: number[] | null;
  equilibrium_concentration: number[] | null;
  r_ln_equilibrium_concentration: number[] | null;
}

export interface ComputationalExtrapolation {
  repeating_units: number[] | null;
  inverse_repeating_units: number[] | null;
  delta_h: number[] | null;
  slope: number | null;
  intercept: number | null;
  std_err: number | null;
}

export interface Thermo {
  delta_h: number | null;
  delta_s: number | null;
  delta_h_std: number | null;
  delta_s_std: number | null;
  delta_g: number | null;
  ceiling_temperature: number | null;
  vanthoff: Vanthoff | null;
  extrapolation: ComputationalExtrapolation | null;
}

export interface Metadata {
  year: number | null;
  comment: string | null;
  doi: string | null;
  url: string | null;
  formatted_reference: string | null;
  flag: string | null;
}

export interface Reaction {
  reaction_id: string;
  type: string;
  monomer: Monomer;
  product: Product;
  parameters: Parameters;
  thermo: Thermo;
  metadata: Metadata;
}

export interface MonomerSummaryState {
  data: MonomerSummary[];
  isLoading: boolean;
}

export interface ReactionState {
  data: Reaction[];
  isLoading: boolean;
}
