from functools import cached_property
from pydantic import (
    BaseModel,
    Field,
    BeforeValidator,
    AfterValidator,
    computed_field,
)
from typing import Optional, Annotated, TypeVar, Any
from roppy.core.constants import STATES, COMP_METHODS, COMP_FIELDS
from roppy.core.solvents import CANONICAL_SOLVENTS, SOLVENT_ALIASES
from roppy.core.efgs import get_dec_fgs
from beanie import Document, Indexed
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from rdkit.Chem.rdinchi import MolToInchi
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles


def validate_solvent(solvent: Optional[str]) -> Optional[str]:
    """Validates and normalizes solvent names."""
    if solvent in CANONICAL_SOLVENTS:
        return solvent
    if solvent in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[solvent]
    return None


def validate_smiles(smiles: str) -> Optional[str]:
    try:
        return StandardizeSmiles(smiles)
    except:
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
Smiles = Annotated[str, AfterValidator(validate_smiles)]


class Monomer(Document):

    smiles: Smiles = Field(
        ...,
        description="SMILES notation representing the monomer structure",
    )

    # @computed_field(description="internal rdkit.Chem.rdchem.mol object")
    # def _mol(self) -> Optional[Mol]:
    #     return MolFromSmiles(self.smiles)

    @computed_field(description="iupac identifier")
    @property
    def inchi(self) -> Optional[str]:
        try:
            return MolToInchi(MolFromSmiles(self.smiles))[0]
        except:
            return None

    @computed_field(description="Molecular weight of the monomer in g/mol")
    @property
    def molecular_weight(self) -> Optional[float]:
        try:
            return CalcExactMolWt(MolFromSmiles(self.smiles))
        except:
            return None

    @computed_field(description="Primary functional group involved in polymerization")
    @property
    def functional_group(self) -> Optional[list[str]]:
        try:
            return get_dec_fgs(MolFromSmiles(self.smiles))[2]
        except:
            return None

    @computed_field(description="Size of the ring in the monomer structure")
    @property
    def ring_size(self) -> Optional[int]:
        try:
            mol = MolFromSmiles(self.smiles)
            ring_info = mol.GetRingInfo()
            n_atoms = mol.GetNumAtoms()
            max_ring_size = 0
            for i in range(n_atoms):
                ring_size = ring_info.MinAtomRingSize(i)
                if ring_size > max_ring_size:
                    max_ring_size = ring_size
            if max_ring_size == 0:
                return None
            return max_ring_size

        except:
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
    repeating_unit: EmptyStringToNone[str] = (
        Field(
            None,
            description="SMILES notation or chemical formula of the repeating unit",
        ),
    )
    n_avg_molar_mass: EmptyStringToNone[float] = Field(
        None, description="arithmetic mean of the molar masses of the polymers"
    )
    m_avg_molar_mass: EmptyStringToNone[float] = Field(
        None, description="measure of the molar mass of the polymer"
    )
    dispersity: EmptyStringToNone[float] = Field(
        None, description="measure of the molar mass distribution of the polymer"
    )
    deg_of_poly: EmptyStringToNone[float] = Field(
        None, description="average number of monomer units per polymer"
    )

    @computed_field(description="Molecular weight of the polymer repeat unit in g/mol")
    @property
    def molecular_weight(self) -> Optional[float]:
        try:
            return CalcExactMolWt(MolFromSmiles(self.repeating_unit))
        except:
            return None

    # compute dispersity and DP from molar masses if present?
    def model_post_init(self, _):
        if self.dispersity is None:
            if isinstance(self.n_avg_molar_mass, float) and isinstance(
                self.m_avg_molar_mass, float
            ):
                self.dispersity = self.m_avg_molar_mass / self.n_avg_molar_mass
        if self.deg_of_poly is None:
            if isinstance(self.n_avg_molar_mass, float) and isinstance(
                self.molecular_weight, float
            ):
                self.deg_of_poly = self.n_avg_molar_mass / self.molecular_weight


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
    monomer_conc: EmptyStringToNone[float] = Field(
        None, description="Initial concentration of the monomer"
    )
    monomer_state: State = Field(None, description="State of the monomer")
    polymer_state: State = Field(None, description="State of the polymer")
    solvent_model: EmptyStringToNone[str] = Field(
        None, description="Computational solvent model used (if applicable)"
    )
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

    @computed_field(description="Formatted summary of monomer-polymer states")
    @property
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
        if self.ceiling_temperature is None:
            if isinstance(self.delta_h, float) and isinstance(self.delta_s, float):
                self.ceiling_temperature = 1000 * self.delta_h / self.delta_s


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

    @computed_field(
        description="Smiles notation representing the polymerization reaction"
    )
    @property
    def reaction_smiles(self) -> Optional[str]:
        return None

    @computed_field(description="Flag indicating if data is experimental")
    @property
    def is_experimental(self) -> bool:
        for comp_field in COMP_FIELDS:
            if getattr(self.parameters, comp_field):
                return False
        return True

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


class MonomerSummary(Document):
    display_id: int = Field(
        ..., description="unique display id for the monomer summary"
    )
    monomer: Monomer = Field(..., description="corresponding monomer")
    # polymerisations: list[Polymerisation] = Field(
    #     ..., description="list of polymerisations for the corresponding monomer"
    # )
    exp_polymerisations: list[Polymerisation] = Field(
        ...,
        description="list of experimental polymerisations for the corresponding monomer",
    )
    comp_polymerisations: list[Polymerisation] = Field(
        ...,
        description="list of computational polymerisations for the corresponding monomer",
    )

    @computed_field(description="Whether the monomer has experimental data")
    @property
    def has_experimental(self) -> bool:
        # for poly in self.polymerisations:
        #     if poly.is_experimental:
        #         return True
        # return False
        return not (len(self.exp_polymerisations) == 0)

    @computed_field(description="Whether the monomer has computational data")
    @property
    def has_computational(self) -> bool:
        # n_experimental = sum(
        #     [1 if poly.is_experimental else 0 for poly in self.polymerisations]
        # )
        # n_poly = len(self.polymerisations)
        # if (n_poly - n_experimental) == 0:
        #     return False
        # return True
        return not (len(self.comp_polymerisations) == 0)


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
