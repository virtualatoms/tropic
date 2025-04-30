from pydantic import BaseModel, Field, field_validator
from typing import Optional, Union, Literal, Annotated
from roppy.core.solvents import CANONICAL_SOLVENTS, SOLVENT_ALIASES

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


class Conditions(BaseModel):
    temperature: Optional[float] = Field(
        None, description="Reaction temperature in Kelvin"
    )
    solvent: Optional[Solvent] = Field(
        None, description="Solvent used in the reaction (must be ORCA-compatible)"
    )
    monomer_state: Optional[Literal["gas", "liquid", "solution"]] = Field(
        None, description="Physical state of the monomer"
    )
    polymer_state: Optional[Literal["gas", "liquid", "solution"]] = Field(
        None, description="Physical state of the polymer"
    )
    pressure: Optional[float] = Field(None, description="Reaction pressure in bar")

    _validate_solvent = field_validator("solvent")(validate_solvent)


class Parameters(BaseModel):
    method: Optional[Literal["dft", "ffmd", "aimd", "mlmd", "xtb", "ml"]] = Field(
        None, description="Computational method used for the calculation"
    )
    functional: Optional[str] = Field(
        None, description="DFT functional used (if applicable)"
    )
    basis_set: Optional[str] = Field(
        None, description="Basis set used in quantum calculations"
    )
    dispersion: Optional[str] = Field(None, description="Dispersion correction method")
    solvent_model: Optional[str] = Field(
        None, description="Computational solvent model used"
    )
    solvent: Optional[Solvent] = Field(
        None, description="Solvent used in calculations (must be ORCA-compatible)"
    )
    forcefield: Optional[str] = Field(
        None, description="Force field used for molecular dynamics"
    )

    _validate_solvent = field_validator("solvent")(validate_solvent)


class Monomer(BaseModel):
    smiles: Optional[str] = Field(
        None, description="SMILES notation representing the monomer structure"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the monomer in g/mol"
    )
    functional_group: Optional[str] = Field(
        None, description="Primary functional group involved in polymerization"
    )


class Product(BaseModel):
    repeating_unit: Optional[str] = Field(
        None, description="SMILES notation or chemical formula of the repeating unit"
    )
    number_of_units: Optional[Union[int, Literal["many"]]] = Field(
        None, description="Number of repeating units in the polymer"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the polymer in g/mol"
    )


class Initiator(BaseModel):
    smiles: Optional[str] = Field(
        None, description="SMILES notation representing the initiator structure"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the initiator in g/mol"
    )
    functional_group: Optional[str] = Field(
        None, description="Primary functional group of the initiator"
    )


class Thermo(BaseModel):
    delta_h: Optional[float] = Field(
        None, description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: Optional[float] = Field(
        None, description="Change in entropy (ΔS) for the polymerization reaction"
    )


class Polymerisation(BaseModel):
    monomer: Monomer = Field(..., description="Monomer information")
    product: Product = Field(..., description="Polymer product information")
    initiator: Optional[Initiator] = Field(None, description="Initiator information")
    thermo: Thermo = Field(..., description="Thermodynamic data")
    conditions: Optional[Conditions] = Field(None, description="Reaction conditions")
    parameters: Optional[Parameters] = Field(
        None, description="Computational parameters"
    )
    experimental: Optional[bool] = Field(
        None, description="Flag indicating if data is experimental"
    )
    comment: Optional[str] = Field(None, description="Additional comments or notes")
    doi: Optional[str] = Field(
        None, description="Digital Object Identifier for the source"
    )
    link: Optional[str] = Field(None, description="URL to related resource")
    polymerisation_id: str = Field(
        ..., description="Unique identifier for the polymerization reaction"
    )


class ExperimentalData(BaseModel):
    state: str = Field(..., description="State of the polymerization reaction")
    solvent_conc: Optional[float] = Field(
        None, description="Solvent concentration in mol/L"
    )
    delta_h: float = Field(
        ..., description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: float = Field(
        ..., description="Change in entropy (ΔS) for the polymerization reaction"
    )
    critical_temperature: float = Field(float, description="Critical temperature in K")
    equilibrium_temperature: float = Field(
        float, description="Equilibrium temperature in K"
    )
    doi: str = Field(..., description="Digital Object Identifier for the source")
    polymerisation_id: str = Field(..., description="ID for the polymerisation")
    initiator_smiles: Optional[str] = Field(
        None, description="SMILES notation representing the initiator"
    )
    number_of_units: Optional[int] = Field(
        None, description="Number of repeating units in the polymer"
    )


class ComputationalData(BaseModel):
    state: str = Field(..., description="State of the polymerization reaction")
    solvent_model: Optional[str] = Field(
        None, description="Method used to model the solvent"
    )
    solvent_model: Optional[str] = Field(None, description="Solvent used")
    delta_h: float = Field(
        ..., description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: float = Field(
        ..., description="Change in entropy (ΔS) for the polymerization reaction"
    )
    critical_temperature: float = Field(float, description="Critical temperature in K")
    equilibrium_temperature: float = Field(
        float, description="Equilibrium temperature in K"
    )
    doi: str = Field(..., description="Digital Object Identifier for the source")
    polymerisation_id: str = Field(..., description="ID for the polymerisation")
    initiator_smiles: Optional[str] = Field(
        None, description="SMILES notation representing the initiator"
    )
    number_of_units: Optional[int] = Field(
        None, description="Number of repeating units in the polymer"
    )


class PolymerisationSummary(BaseModel):
    reaction_smiles: str = Field(
        ..., description="SMILES notation representing the polymerization reaction"
    )
    polyinfo_id: Optional[str] = Field(None, description="ID in PolyInfo database")
    polygenome_id: Optional[str] = Field(None, description="ID in polygenome database")
    experimental_data: list[ExperimentalData] = Field(
        default_factory=list,
        description="Experimental data for the polymerization reaction",
    )
    computational_data: list[ComputationalData] = Field(
        default_factory=list,
        description="Computational data for the polymerization reaction",
    )
    polymerisation_summary_id: str = Field(
        ..., description="Unique identifier for the polymerization summary"
    )


class MonomerSummary(BaseModel):
    monomer_id: str = Field(..., description="Unique identifier for the monomer")
    smiles: str = Field(
        ..., description="SMILES notation representing the monomer structure"
    )
    ring_size: int = Field(
        ...,
        description="Size of the ring in the monomer structure",
    )
    inchi: Optional[str] = Field(
        None, description="IUPAC International Chemical Identifier"
    )
    iupac_name: Optional[str] = Field(None, description="IUPAC name of the monomer")
    common_name: Optional[str] = Field(None, description="Common name of the monomer")
    xyz: Optional[str] = Field(None, description="XYZ coordinates for the monomer")
    has_exp: bool = Field(..., description="Whether the monomer has experimental data")
    has_calc: bool = Field(..., description="Whether the monomer has calculated data")
    polymerisation: PolymerisationSummary = Field(
        None,
        description="Polymerization summary involving the monomer",
    )
    ring_opening: list[PolymerisationSummary] = Field(
        default_factory=list,
        description="Ring opening summaries involving the monomer",
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
