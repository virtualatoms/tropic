from functools import cached_property
from pydantic import (
    BaseModel,
    Field,
    BeforeValidator,
    AfterValidator,
    computed_field,
)
from typing import Optional, Annotated, TypeVar
from roppy.core.constants import STATES, COMP_METHODS
from roppy.core.solvents import CANONICAL_SOLVENTS, SOLVENT_ALIASES
from beanie import Document, Indexed


def validate_solvent(solvent: Optional[str]) -> Optional[str]:
    """Validates and normalizes solvent names."""
    if solvent in CANONICAL_SOLVENTS:
        return solvent
    if solvent in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[solvent]
    return None


EmptyStringToNone = Annotated[
    Optional[TypeVar("T")],
    BeforeValidator(lambda s: s if s else None),
]
State = Annotated[Optional[str], BeforeValidator(lambda s: s if s in STATES else None)]
CompMethod = Annotated[
    Optional[str], BeforeValidator(lambda s: s if s in COMP_METHODS else None)
]
Solvent = Annotated[Optional[str], AfterValidator(validate_solvent)]


class Monomer(Document):

    smiles: str = Field(
        # need to perform post-validation normalisation to ensure bijectivity
        ...,
        description="SMILES notation representing the monomer structure",
    )

    @cached_property
    @computed_field(description="iupac identifier")
    def inchi(self) -> Optional[str]:
        return ""

    @cached_property
    @computed_field(description="Molecular weight of the monomer in g/mol")
    def molecular_weight(self) -> Optional[float]:
        return None

    @cached_property
    @computed_field(description="Primary functional group involved in polymerization")
    def functional_group(self) -> Optional[str]:
        return None

    @cached_property
    @computed_field(description="Size of the ring in the monomer structure")
    def ring_size(self) -> Optional[int]:
        return None

    @cached_property
    @computed_field(description="IUPAC name of the monomer")
    def iupac_name(self) -> Optional[str]:
        return None

    @cached_property
    @computed_field(description="Common name of the monomer")
    def common_name(self) -> Optional[str]:
        return None

    @cached_property
    @computed_field(description="XYZ coordinates for the monomer")
    def xyz(self) -> Optional[str]:
        return None

    class Settings:
        name = "monomers"


class Initiator(BaseModel):
    smiles: str = Field(
        ..., description="SMILES notation representing the initiator structure"
    )

    # TODO: Add same fields as Monomer?


class Product(BaseModel):
    number_of_units: EmptyStringToNone[float] = (
        Field(..., description="(average) Number of repeating units in the polymer"),
    )
    repeating_unit: EmptyStringToNone[str] = (
        Field(
            ..., description="SMILES notation or chemical formula of the repeating unit"
        ),
    )

    @cached_property
    @computed_field(description="Average molecular weight of the polymer in g/mol")
    def molecular_weight(self) -> Optional[float]:
        return None

    """ 
    TODO: 
    Add additional descriptors for product?
    E.g. dispersity, average molar mass, etc.

        """


class Parameters(BaseModel):
    temperature: EmptyStringToNone[float] = Field(
        None, description="Temperature of the polymerisation in K"
    )
    pressure: EmptyStringToNone[float] = Field(
        None, description="Pressure of the polymerisation in K"
    )
    solvent: Solvent = Field(None, description="Solvent medium of the polymerisation")
    solvent_conc: EmptyStringToNone[float] = Field(
        None, description="Concentration of the polymerisation solvent"
    )
    solvent_model: EmptyStringToNone[str] = Field(
        None, description="Computational solvent model used (if applicable)"
    )
    monomer_state: State = Field(None, description="State of the monomer")
    polymer_state: State = Field(None, description="State of the polymer")
    method: CompMethod = Field(None, description="Computational method used")
    functional: EmptyStringToNone[str] = Field(
        None, description="DFT functional used (if applicable)"
    )
    basis_set: EmptyStringToNone[str] = Field(
        None, description="Basis set used in quantum calculations"
    )
    dispersion: EmptyStringToNone[str] = Field(
        None, description="Dispersion correction method"
    )
    forcefield: EmptyStringToNone[str] = Field(
        None, description="Force field used for molecular dynamics"
    )


class Thermo(BaseModel):
    delta_h: EmptyStringToNone[float] = (
        Field(
            None, description="Change in enthalpy (ΔH) for the polymerization reaction"
        ),
    )
    delta_s: EmptyStringToNone[float] = (
        Field(
            None, description="Change in entropy (ΔS) for the polymerization reaction"
        ),
    )
    ceiling_temperature: EmptyStringToNone[float] = (
        Field(None, description="Ceiling temperature in K"),
    )

    # compute ceiling temperature from H and S if present?
    def model_post_init(self, _):
        pass


class Metadata(BaseModel):
    comment: EmptyStringToNone[str] = (
        Field(..., description="Additional comments or notes"),
    )
    doi: EmptyStringToNone[str] = (
        Field(..., description="Digital Object Identifier for the source"),
    )
    url: EmptyStringToNone[str] = (Field(..., description="URL to related resource"),)


class Polymerisation(Document):

    display_id: Optional[int] = Field(
        None, description="unique display id for polymerisation"
    )

    monomer: Monomer = Field(..., description="Monomer of the polymerisation")
    initiator: Optional[Initiator] = Field(
        ..., description="Initiator of the polymerisation"
    )
    product: Product = Field(..., description="Product of the polymerisation")
    parameters: Parameters = Field(..., description="Computational parameters")
    thermo: Thermo = Field(..., description="Thermodynamic data/results")
    metadata: Metadata = Field(
        ..., description="Polymerisation references and metadata"
    )

    # polyinfo_id: Optional[str] = Field(None, description="ID in PolyInfo database")
    # polygenome_id: Optional[str] = Field(None, description="ID in polygenome database")

    @cached_property
    @computed_field(
        description="Smiles notation representing the polymerization reaction"
    )
    def reaction_smiles(self) -> Optional[str]:
        return None

    @cached_property
    @computed_field(description="Flag indicating if data is experimental")
    def is_experimental(self) -> bool:
        return False

    def model_post_init(self, _) -> None:
        if isinstance(self.initiator, Initiator):
            if not self.initiator.smiles:
                self.initiator = None

    class Settings:
        name = "polymerisations"


class DataRow(BaseModel):
    monomer_state: Optional[str] = Field(..., description="State of the monomer")
    polymer_state: Optional[str] = Field(..., description="State of the polymer")
    solvent: Optional[str] = Field(
        ..., description="Solvent medium of the polymerisation"
    )
    solvent_conc: Optional[str] = Field(
        ..., description="Concentration of the polymerisation solvent"
    )
    delta_h: Optional[float] = Field(
        ..., description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: Optional[float] = Field(
        ..., description="Change in entropy (ΔS) for the polymerization reaction"
    )
    ceiling_temperature: Optional[float] = Field(
        ..., description="Ceiling temperature in K"
    )

    @cached_property
    @computed_field(description="Whether the monomer has experimental data")
    def state_summary(self) -> Optional[str]:
        if not (self.monomer_state or self.polymer_state):
            return None
        if self.monomer_state and self.polymer_state:
            return f"{self.monomer_state}-{self.polymer_state}"
        else:
            if self.monomer_state:
                return f"{self.monomer_state}-x"
            else:
                return f"x-{self.polymer_state}"


class MonomerSummary(Document):
    display_id: int = Field(
        ..., description="unique display id for the monomer summary"
    )
    monomer: Monomer = Field(..., description="corresponding monomer")
    polymerisations: list[Polymerisation]

    @cached_property
    @computed_field(description="Whether the monomer has experimental data")
    def has_experimental(self) -> bool:
        return False

    @cached_property
    @computed_field(description="Whether the monomer has computational data")
    def has_computational(self) -> bool:
        return False


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
