"""Definition of TROPIC data models."""

from pydantic import BaseModel, Field

from tropic.core.validate import (
    Medium,
    Method,
    MethodCalc,
    ReactionType,
    Smiles,
    Solvent,
    State,
    Topology,
    get_func_groups,
    get_inchi,
    get_molecular_weight,
    get_ring_size,
    get_xyz,
)


class Monomer(BaseModel):
    """Model representing a monomer in the TROPIC database."""

    monomer_id: str = Field(
        "monomer-0",
        description="Unique display id for the monomer",
    )
    smiles: Smiles = Field(
        description="SMILES notation representing the molecular structure",
    )
    inchi: str | None = Field(
        description="Inchi identifier for the molecule",
        default_factory=lambda data: get_inchi(data["smiles"]),
    )
    molecular_weight: float | None = Field(
        description="Molecular weight of the molecule in g/mol",
        default_factory=lambda data: get_molecular_weight(data["smiles"]),
    )
    functional_groups: list[str] | None = Field(
        description="Primary functional group of molecule",
        default_factory=lambda data: get_func_groups(data["smiles"]),
    )
    iupac_name: str | None = Field(
        None,
        description="IUPAC name of the molecule",
    )
    pubchem_cid: int | None = Field(
        None,
        description="PubChem Compound ID (CID) for the molecule",
    )
    ring_size: int | None = Field(
        description="Size of the ring in the monomer structure",
        default_factory=lambda data: get_ring_size(data["smiles"]),
    )
    xyz: str | None = Field(
        description="XYZ coordinates for the molecule",
        default_factory=lambda data: get_xyz(data["smiles"]),
    )


class Product(BaseModel):
    """Model representing a polymer product in the TROPIC database."""

    smiles: str | None = Field(
        None,
        description="Big SMILES notation representing the polymer",
    )
    repeating_units: float | None = Field(
        None,
        description="Number of repeating monomer units in the polymer chain",
    )
    deg_of_poly: float | None = Field(
        None,
        description="average number of monomer units per polymer",
    )
    dispersity: float | None = Field(
        None,
        description="measure of the molar mass distribution of the polymer",
    )
    n_avg_molar_mass: float | None = Field(
        None,
        description="arithmetic mean of the molar masses of the polymers",
    )
    m_avg_molar_mass: float | None = Field(
        None,
        description="measure of the molar mass of the polymer",
    )


class Parameters(BaseModel):
    """Model representing experimental/computational parameters of a reaction."""

    is_experimental: bool = Field(
        description="Whether the reaction is experimental or computational",
    )
    temperature: float | None = Field(
        None,
        description="Temperature of the reaction in K",
    )
    pressure: float | None = Field(
        None,
        description="Pressure of the reaction (if not standard)",
    )
    monomer_state: State | None = Field(None, description="State of the monomer")
    polymer_state: State | None = Field(None, description="State of the polymer")
    initiator_smiles: str | None = Field(
        None,
        description="Initiator of the reaction",
    )
    initial_monomer_conc: float | None = Field(
        None,
        description="Initial concentration of the monomer",
    )
    bulk_monomer_conc: float | None = Field(
        None,
        description="Bulk concentration of the monomer",
    )
    medium: Medium | None = Field(
        None,
        description="Medium that the reaction is conducted within",
    )
    solvent: Solvent | None = Field(
        None,
        description="Solvent used in the reaction",
    )
    cosolvent: Solvent | None = Field(
        None,
        description="Co-solvent used in the reaction",
    )
    solvent_cosolvent_ratio: str | None = Field(
        None,
        description="Ratio of solvent to co-solvent used in the reaction",
    )
    topology: Topology | None = Field(
        None,
        description="Topology used to model the polymer structure",
    )
    method: Method | None = Field(None, description="Computational method used")
    method_calc: MethodCalc | None = Field(
        None,
        description="Approach used to obtain the thermodynamic parameters",
    )
    functional: str | None = Field(None, description="DFT functional")
    basis_set: str | None = Field(None, description="Basis set")
    dispersion: str | None = Field(None, description="Dispersion correction method")
    forcefield: str | None = Field(None, description="Molecule dynamics force field")
    solvent_model: str | None = Field(None, description="Solvent model")
    state_summary: str = Field(
        description="Formatted summary of monomer-polymer states",
        default_factory=lambda data: f"{data.get('monomer_state', 'x')}-{data.get('polymer_state', 'x')}",  # noqa: E501
    )


class Vanthoff(BaseModel):
    """Model representing Van't Hoff data for a reaction."""

    temperature: list[float] | None = Field(
        None,
        description="Temperatures at which the reaction was conducted (K)",
    )
    inverse_temperature: list[float] | None = Field(
        None,
        description="Inverse temperatures (1/T) (1/K)",
    )
    equilibrium_concentration: list[float] | None = Field(
        None,
        description="Equilibrium concentrations of the reactants/products (M)",
    )
    r_ln_equilibrium_concentration: list[float] | None = Field(
        None,
        description="R * ln(equilibrium concentration)",
    )


class ComputationalExtrapolation(BaseModel):
    """Model representing computational extrapolation data for a reaction."""

    repeating_units: list[float] | None = Field(
        None,
        description="List of repeating units",
    )
    inverse_repeating_units: list[float] | None = Field(
        None,
        description="List of 1 / repeating units",
    )
    delta_h: list[float] | None = Field(
        None,
        description="List of delta H values",
    )
    slope: float | None = Field(
        None,
        description="Slope of the extrapolation",
    )
    intercept: float | None = Field(
        None,
        description="Intercept of the extrapolation",
    )
    std_err: float | None = Field(
        None,
        description="Standard error of the extrapolation",
    )


class Thermo(BaseModel):
    """Model representing thermodynamic data for a reaction."""

    delta_h: float | None = Field(
        None,
        description="Enthalpy of polymerisation (kJ/mol)",
    )
    delta_s: float | None = Field(
        None,
        description="Entropy of polymerisation (J/molK)",
    )
    delta_h_std: float | None = Field(
        None,
        description="Uncertainty in enthalpy of polymerisation (kJ/mol)",
    )
    delta_s_std: float | None = Field(
        None,
        description="Uncertainty in entropy of polymerisation (J/molK)",
    )
    delta_g: float | None = Field(
        None,
        description="Free energy of polymerisation (kJ/mol)",
    )
    ceiling_temperature: float | None = Field(
        None,
        description="Ceiling temperature in C",
    )
    vanthoff: Vanthoff | None = Field(
        None,
        description="Van't Hoff data for the reaction",
    )
    extrapolation: ComputationalExtrapolation | None = Field(
        None,
        description="Computational extrapolation data for the reaction",
    )


class Metadata(BaseModel):
    """Model representing metadata for a reaction."""

    year: int | None = Field(description="Year of publication")
    comment: str | None = Field(description="Additional comments or notes")
    doi: str | None = Field(description="Digital Object Identifier for the source")
    url: str | None = Field(description="URL to related resource")
    formatted_reference: str | None = Field(
        None,
        description="Formatted reference string for the publication",
    )
    flag: str | None = Field(
        None,
        description="Any potential issues with the reaction data",
    )


class Reaction(BaseModel):
    """Model representing a reaction in the TROPIC database."""

    reaction_id: str = Field(
        "reaction-0",
        description="unique display id for the reaction",
    )
    type: ReactionType = Field(description="Type of reaction")
    monomer: Monomer = Field(description="Monomer of the reaction")
    product: Product = Field(description="Product of the reaction")
    parameters: Parameters = Field(description="Experimental/Computational parameters")
    thermo: Thermo = Field(description="Thermodynamic data/results")
    metadata: Metadata = Field(description="Reaction references and metadata")
