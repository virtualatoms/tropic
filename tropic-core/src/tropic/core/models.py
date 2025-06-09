from pydantic import BaseModel, Field

from tropic.core.validate import (
    Method,
    PolymerisationType,
    Smiles,
    Solvent,
    State,
    get_ring_size,
    get_xyz,
    get_molecular_weight,
    get_func_groups,
    get_inchi
)


class Monomer(BaseModel):
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
    smiles: str | None = Field(
        None,
        description="Big SMILES notation representing the polymer",
    )
    repeating_units: int | None = Field(
        None,
        description="Number of repeating monomer units in the computational polymer chain/ring",
    )
    deg_of_poly: float | None = Field(
        None, description="average number of monomer units per polymer"
    )
    dispersity: float | None = Field(
        None, description="measure of the molar mass distribution of the polymer"
    )
    n_avg_molar_mass: float | None = Field(
        None, description="arithmetic mean of the molar masses of the polymers"
    )
    m_avg_molar_mass: float | None = Field(
        None, description="measure of the molar mass of the polymer"
    )


class Parameters(BaseModel):
    is_experimental: bool = Field(
        description="Flag indicating whether the reaction is experimental (True) or computational (False)",
    )
    temperature: float | None = Field(
        None, description="Temperature of the polymerisation in C"
    )
    pressure: float | None = Field(
        None, description="Pressure of the polymerisation (if not standard)"
    )
    monomer_state: State | None = Field(None, description="State of the monomer")
    polymer_state: State | None = Field(None, description="State of the polymer")
    initiator_smiles: str | None = Field(None, description="Initiator of the polymerisation")
    initial_monomer_conc: float | None = Field(
        None, description="Initial concentration of the monomer"
    )
    bulk_monomer_conc: float | None = Field(
        None, description="Bulk concentration of the monomer"
    )
    solvent: Solvent | None = Field(
        None, description="Solvent that the polymerisation is conducted within"
    )
    method: Method | None = Field(None, description="Computational method used")
    functional: str | None = Field(None, description="DFT functional")
    basis_set: str | None = Field(None, description="Basis set")
    dispersion: str | None = Field(
        None, description="Dispersion correction method"
    )
    forcefield: str | None = Field(
        None, description="Molecule dynamics force field"
    )
    solvent_model: str | None = Field(None, description="Solvent model")
    state_summary: str = Field(
        description="Formatted summary of monomer-polymer states",
        default_factory=lambda data: f"{data['monomer_state'] or 'x'}-{data['polymer_state'] or 'x'}",
    )


class Thermo(BaseModel):
    delta_h: float | None = Field(None, description="Enthalpy of polymerisation (kJ/mol)")
    delta_s: float | None = Field(None, description="Entropy of polymerisation (J/molK)")
    delta_g: float | None = Field(None, description="Free energy of polymerisation (kJ/mol)")
    ceiling_temperature: float | None = Field(None, description="Ceiling temperature in C")


class Metadata(BaseModel):
    year: int | None = Field(description="Year of publication")
    comment: str | None = Field(description="Additional comments or notes")
    doi: str | None = Field(description="Digital Object Identifier for the source")
    url: str | None = Field(description="URL to related resource")
    formatted_reference: str | None = Field(
        None,
        description="Formatted reference string for the publication",
    )


class Polymerisation(BaseModel):
    polymerisation_id: str = Field(
        "poly-0", description="unique display id for the polymerisation"
    )
    type: PolymerisationType = Field(description="Type of polymerisation")
    monomer: Monomer = Field(description="Monomer of the polymerisation")
    product: Product = Field(description="Product of the polymerisation")
    parameters: Parameters = Field(description="Experimental/Computational parameters")
    thermo: Thermo = Field(description="Thermodynamic data/results")
    metadata: Metadata = Field(description="Polymerisation references and metadata")


class DataRow(BaseModel):
    type: PolymerisationType = Field(description="Type of polymerisation")
    polymerisation_id: str = Field(
        description="Unique display id for the polymerisation",
    )
    is_experimental: bool = Field(
        description="Flag indicating whether the reaction is experimental (True) or computational (False)",
    )
    state_summary: str = Field(
        description="Formatted summary of monomer-polymer states",
    )
    initial_monomer_conc: float | None = Field(
        description="Initial concentration of the monomer"
    )
    bulk_monomer_conc: float | None = Field(
        description="bulk concentration of the monomer"
    )
    solvent: str | None = Field(
        description="Solvent that the polymerisation is conducted within"
    )
    delta_h: float | None = Field(
        description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: float | None = Field(
        description="Change in entropy (ΔS) for the polymerization reaction"
    )
    repeating_units: int | None = Field(
        description="Number of repeating monomer units in the computational polymer chain/ring"
    )
    method: str | None = Field(
        description="Computational method used for the polymerisation"
    )
    ceiling_temperature: float | None = Field(description="Ceiling temperature in K")
    year: int | None = Field(description="Year of publication")
    doi: str | None = Field(description="DOI of publication")
    formatted_reference: str | None = Field(description="Formatted reference string")


class MonomerSummary(BaseModel):
    monomer_id: str = Field(description="unique display id for the monomer summary")
    monomer: Monomer = Field(description="corresponding monomer")
    data: list[DataRow] = Field(
        description="table of data where each row corresponds to a polymerisation (for display purposes)",
    )
    has_exp: bool = Field(
        description="Whether the molecule has experimental data",
        default_factory=lambda data: any(
            poly.is_experimental for poly in data.get("data", [])
        ),
    )
    has_calc: bool = Field(
        description="Whether the molecule has calculated data",
        default_factory=lambda data: any(
            not poly.is_experimental for poly in data.get("data", [])
        ),
    )
