from pydantic import BaseModel, Field, field_validator
from typing import Optional, Union, Literal, Annotated
from solvents import CANONICAL_SOLVENTS, SOLVENT_ALIASES

Solvent = Annotated[str, Field(description="Solvent name (canonical or alias)")]


def validate_solvent(value: Optional[str]) -> Optional[str]:
    """Validates and normalizes solvent names."""
    if value is None:
        return None

    if value in CANONICAL_SOLVENTS:
        return value

    if value in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[value]

    raise ValueError(
        f"Invalid solvent: {value}. Must be a valid ORCA solvent name or alias."
    )


class Monomer(BaseModel):
    monomer_id: int = Field(..., description="Unique identifier for the monomer")
    smiles: str = Field(
        ..., description="SMILES notation representing the monomer structure"
    )
    inchi: Optional[str] = Field(
        None, description="IUPAC International Chemical Identifier"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the monomer in g/mol"
    )
    functional_group: Optional[str] = Field(
        None, description="Primary functional group involved in polymerization"
    )
    ring_size: Optional[int] = Field(
        None,
        description="Size of the ring in the monomer structure",
    )
    iupac_name: Optional[str] = Field(None, description="IUPAC name of the monomer")
    common_name: Optional[str] = Field(None, description="Common name of the monomer")
    xyz: Optional[str] = Field(None, description="XYZ coordinates for the monomer")
    polymerisations: list[int] = Field(
        default_factory=list,
        description="List of polymerization ids involving the monomer (including ring-openings)",
    )

    # def __init__(self, monomer_id: int, smiles: str, polymerisations: list[int]):
    #     self.monomer_id = monomer_id
    #     self.smiles = smiles
    #     self.polymerisations = polymerisations
    #
    #     # TODO: complete rest of initialisation script to generate other fields


class Initiator(BaseModel):
    initiator_id: int = Field(..., description="Unique identifier for the initiator")
    smiles: str = Field(
        ..., description="SMILES notation representing the initiator structure"
    )
    inchi: Optional[str] = Field(
        None, description="IUPAC International Chemical Identifier"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the initiator in g/mol"
    )
    functional_group: Optional[str] = Field(
        None, description="Primary functional group of the initiator"
    )
    iupac_name: Optional[str] = Field(None, description="IUPAC name of the initiator")
    common_name: Optional[str] = Field(None, description="Common name of the initiator")
    xyz: Optional[str] = Field(None, description="XYZ coordinates for the initiator")
    polymerisations: list[int] = Field(
        default_factory=list,
        description="List of polymerization ids involving the initiator (including ring-openings)",
    )

    # def __init__(self, initiator_id: int, smiles: str, polymerisations: list[int]):
    #     self.initiator_id = initiator_id
    #     self.smiles = smiles
    #     self.polymerisations = polymerisations

    # TODO: complete rest of initialisation script to generate other fields


class Product(BaseModel):
    number_of_units: float = Field(
        ..., description="(average) Number of repeating units in the polymer"
    )
    repeating_unit: Optional[str] = Field(
        None, description="SMILES notation or chemical formula of the repeating unit"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Average molecular weight of the polymer in g/mol"
    )
    """ 
    TODO: 
    Add additional descriptors for product?
    E.g. dispersity, average molar mass, etc.

    Complete initialisation script to generate other fields
        """

    # def __init__(self, number_of_units: float):
    #     self.number_of_units = number_of_units

    # TODO: complete rest of initialisation script to generate other fields


class Parameters(BaseModel):
    temperature: Optional[float] = Field(
        None, description="Reaction temperature in Kelvin"
    )
    pressure: Optional[float] = Field(None, description="Reaction pressure in bar")
    solvent: Optional[Solvent] = Field(
        None, description="Solvent used in the reaction (must be ORCA-compatible)"
    )
    solvent_conc: Optional[float] = Field(
        None, description="Solvent concentration in mol/L"
    )
    monomer_state: Optional[Literal["a", "g", "l", "s", "c"]] = Field(
        None, description="Physical state of the monomer"
    )
    polymer_state: Optional[Literal["a", "g", "l", "s", "c"]] = Field(
        None, description="Physical state of the polymer"
    )
    method: Optional[Literal["dft", "ffmd", "aimd", "mlmd", "xtb", "ml"]] = Field(
        None, description="Computational method used for the calculation"
    )
    solvent_model: Optional[str] = Field(
        None, description="Computational solvent model used"
    )
    functional: Optional[str] = Field(
        None, description="DFT functional used (if applicable)"
    )
    basis_set: Optional[str] = Field(
        None, description="Basis set used in quantum calculations"
    )
    dispersion: Optional[str] = Field(None, description="Dispersion correction method")
    forcefield: Optional[str] = Field(
        None, description="Force field used for molecular dynamics"
    )

    _validate_solvent = field_validator("solvent")(validate_solvent)


class Thermo(BaseModel):
    delta_h: Optional[float] = Field(
        None, description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: Optional[float] = Field(
        None, description="Change in entropy (ΔS) for the polymerization reaction"
    )
    ceiling_temperature: Optional[float] = Field(
        None, description="Ceiling temperature in K"
    )


class Polymerisation(BaseModel):
    polymerisation_id: int = Field(
        ..., description="Unique identifier for the polymerization reaction"
    )
    product: Product = Field(..., description="Product of the polymerisation")
    is_experimental: bool = Field(
        ..., description="Flag indicating if data is experimental"
    )
    parameters: Parameters = Field(..., description="Computational parameters")
    thermo: Thermo = Field(..., description="Thermodynamic data/results")
    comment: Optional[str] = Field(None, description="Additional comments or notes")
    doi: Optional[str] = Field(
        None, description="Digital Object Identifier for the source"
    )
    url: Optional[str] = Field(None, description="URL to related resource")
    monomer_id: int = Field(..., description="Monomer id")
    initiator_id: Optional[int] = Field(None, description="Initiator id")


class PolymerisationData(BaseModel):
    is_exp: bool = Field(
        ..., description="Whether the reaction is experimental or computational"
    )
    state: Optional[str] = Field(
        None, description="State of the polymerization reaction"
    )
    number_of_units: Optional[int] = Field(
        None, description="Number of repeating units in the polymer"
    )
    initiator_smiles: Optional[str] = Field(
        None, description="SMILES notation representing the initiator"
    )
    delta_h: Optional[float] = Field(
        None, description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: Optional[float] = Field(
        None, description="Change in entropy (ΔS) for the polymerization reaction"
    )
    critical_temperature: Optional[float] = Field(
        None, description="Critical temperature in K"
    )
    equilibrium_temperature: Optional[float] = Field(
        None, description="Equilibrium temperature in K"
    )
    doi: Optional[str] = Field(
        None, description="Digital Object Identifier for the source"
    )
    solvent_conc: Optional[float] = Field(
        None, description="Solvent concentration in mol/L"
    )
    solvent_model: Optional[str] = Field(
        None, description="Method used to model the solvent"
    )
    polymerisation_id: str = Field(..., description="ID for the polymerisation")


class PolymerisationSummary(BaseModel):
    reaction_smiles: str = Field(
        ..., description="SMILES notation representing the polymerization reaction"
    )
    product: Product = Field(..., description="Polymer product information")
    is_ring_opening: bool = Field(
        ..., description="Whether the reaction is a ring-opening polymerization"
    )
    polyinfo_id: Optional[str] = Field(None, description="ID in PolyInfo database")
    polygenome_id: Optional[str] = Field(None, description="ID in polygenome database")
    data: list[PolymerisationData] = Field(
        default_factory=list,
        description="Data for the polymerization reaction",
    )
    polymerisation_summary_id: str = Field(
        ..., description="Unique identifier for the polymerization summary"
    )


class MonomerSummaryBrief(BaseModel):
    monomer_id: str = Field(..., description="Unique identifier for the monomer")
    smiles: str = Field(
        ..., description="SMILES notation representing the monomer structure"
    )
    ring_size: int = Field(
        ...,
        description="Size of the ring in the monomer structure",
    )
    has_exp: bool = Field(..., description="Whether the monomer has experimental data")
    has_calc: bool = Field(..., description="Whether the monomer has calculated data")
