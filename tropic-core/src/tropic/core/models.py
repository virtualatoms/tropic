"""Definition of TROPIC data models."""

from pydantic import BaseModel, Field

from tropic.core.validate import (
    Medium,
    Method,
    PolymerisationType,
    Smiles,
    State,
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
    repeating_units: int | None = Field(
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
    """Model representing experimental/computational parameters of a polymerisation."""

    is_experimental: bool = Field(
        description="Whether the reaction is experimental or computational",
    )
    temperature: float | None = Field(
        None,
        description="Temperature of the polymerisation in C",
    )
    pressure: float | None = Field(
        None,
        description="Pressure of the polymerisation (if not standard)",
    )
    monomer_state: State | None = Field(None, description="State of the monomer")
    polymer_state: State | None = Field(None, description="State of the polymer")
    initiator_smiles: str | None = Field(
        None,
        description="Initiator of the polymerisation",
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
        description="Medium that the polymerisation is conducted within",
    )
    method: Method | None = Field(None, description="Computational method used")
    functional: str | None = Field(None, description="DFT functional")
    basis_set: str | None = Field(None, description="Basis set")
    dispersion: str | None = Field(None, description="Dispersion correction method")
    forcefield: str | None = Field(None, description="Molecule dynamics force field")
    solvent_model: str | None = Field(None, description="Solvent model")
    state_summary: str = Field(
        description="Formatted summary of monomer-polymer states",
        default_factory=lambda data: f"{data.get('monomer_state', 'x')}-{data.get('polymer_state', 'x')}",  # noqa: E501
    )


class Thermo(BaseModel):
    """Model representing thermodynamic data for a polymerisation."""

    delta_h: float | None = Field(
        None,
        description="Enthalpy of polymerisation (kJ/mol)",
    )
    delta_s: float | None = Field(
        None,
        description="Entropy of polymerisation (J/molK)",
    )
    delta_g: float | None = Field(
        None,
        description="Free energy of polymerisation (kJ/mol)",
    )
    ceiling_temperature: float | None = Field(
        None,
        description="Ceiling temperature in C",
    )


class Metadata(BaseModel):
    """Model representing metadata for a polymerisation."""

    year: int | None = Field(description="Year of publication")
    comment: str | None = Field(description="Additional comments or notes")
    doi: str | None = Field(description="Digital Object Identifier for the source")
    url: str | None = Field(description="URL to related resource")
    formatted_reference: str | None = Field(
        None,
        description="Formatted reference string for the publication",
    )


class Polymerisation(BaseModel):
    """Model representing a polymerisation in the TROPIC database."""

    polymerisation_id: str = Field(
        "poly-0",
        description="unique display id for the polymerisation",
    )
    type: PolymerisationType = Field(description="Type of polymerisation")
    monomer: Monomer = Field(description="Monomer of the polymerisation")
    product: Product = Field(description="Product of the polymerisation")
    parameters: Parameters = Field(description="Experimental/Computational parameters")
    thermo: Thermo = Field(description="Thermodynamic data/results")
    metadata: Metadata = Field(description="Polymerisation references and metadata")
