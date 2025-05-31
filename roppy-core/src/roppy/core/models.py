from functools import cached_property
from pydantic import (
    BaseModel,
    Field,
    computed_field,
)
from typing import Optional
from roppy.core.validate import (
    EmptyStringToNone,
    Smiles,
    PolymerisationType,
    Solvent,
    State,
    Method,
)
from roppy.core.efgs import get_dec_fgs
from beanie import Document, Indexed
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from rdkit.Chem.rdinchi import MolToInchi


class Molecule(Document):

    smiles: Smiles = Field(
        ...,
        description="SMILES notation representing the molecular structure",
    )

    @computed_field(description="Inchi identifier")
    @property
    def inchi(self) -> Optional[str]:
        try:
            return MolToInchi(MolFromSmiles(self.smiles))[0]
        except:
            return None

    @computed_field(description="Molecular weight of the molecule in g/mol")
    @property
    def molecular_weight(self) -> Optional[float]:
        try:
            return CalcExactMolWt(MolFromSmiles(self.smiles))
        except:
            return None

    @computed_field(description="Primary functional group of molecule")
    @property
    def functional_groups(self) -> Optional[list[str]]:
        try:
            return get_dec_fgs(MolFromSmiles(self.smiles))[2]
        except:
            return None

    # TODO: Remaining fields for molecule parent class
    @cached_property
    @computed_field(description="IUPAC name of the molecule")
    def iupac_name(self) -> Optional[str]:
        return None

    @cached_property
    @computed_field(description="Common name of the molecule")
    def common_name(self) -> Optional[str]:
        return None

    @cached_property
    @computed_field(description="XYZ coordinates for the molecule")
    def xyz(self) -> Optional[str]:
        return None


class Monomer(Molecule):

    @computed_field(description="Size of the ring in the monomer structure")
    @property
    def ring_size(self) -> Optional[int]:
        try:
            mol = MolFromSmiles(self.smiles)
            ring_info = mol.GetRingInfo()
            n_atoms = mol.GetNumAtoms()
            max_ring_size = max([ring_info.MinAtomRingSize(i) for i in range(n_atoms)])
            if max_ring_size == 0:
                return None
            return max_ring_size

        except:
            return None

    class Settings:
        name = "monomers"


class Initiator(Molecule):

    class Settings:
        name = "initiators"


class Product(BaseModel):
    deg_of_poly: EmptyStringToNone[float] = Field(
        None, description="average number of monomer units per polymer"
    )
    dispersity: EmptyStringToNone[float] = Field(
        None, description="measure of the molar mass distribution of the polymer"
    )
    n_avg_molar_mass: EmptyStringToNone[float] = Field(
        None, description="arithmetic mean of the molar masses of the polymers"
    )
    m_avg_molar_mass: EmptyStringToNone[float] = Field(
        None, description="measure of the molar mass of the polymer"
    )

    # compute dispersity and DP from molar masses if present?
    # def model_post_init(self, _):
    #     if self.dispersity is None:
    #         if isinstance(self.n_avg_molar_mass, float) and isinstance(
    #             self.m_avg_molar_mass, float
    #         ):
    #             self.dispersity = self.m_avg_molar_mass / self.n_avg_molar_mass
    #     if self.deg_of_poly is None:
    #         if isinstance(self.n_avg_molar_mass, float) and isinstance(
    #             self.molecular_weight, float
    #         ):
    #             self.deg_of_poly = self.n_avg_molar_mass / self.molecular_weight


class Parameters(BaseModel):
    is_experimental: bool = Field(
        ...,
        description="Flag indicating whether the reaction is experimental (True) or computational (False)",
    )
    temperature: EmptyStringToNone[float] = Field(
        None, description="Temperature of the polymerisation in C"
    )
    pressure: EmptyStringToNone[float] = Field(
        None, description="Pressure of the polymerisation (if not standard)"
    )
    monomer_state: State = Field(None, description="State of the monomer")
    polymer_state: State = Field(None, description="State of the polymer")
    monomer_conc: EmptyStringToNone[float] = Field(
        None, description="Initial concentration of the monomer"
    )
    solvent: Solvent = Field(
        None, description="Solvent that the polymerisation is conducted within"
    )
    solvent_conc: EmptyStringToNone[float] = Field(
        None, description="Concentration of the polymerisation solvent"
    )
    method: Method = Field(None, description="Computational method used")
    functional: EmptyStringToNone[str] = Field(None, description="DFT functional")
    basis_set: EmptyStringToNone[str] = Field(None, description="Basis set")
    dispersion: EmptyStringToNone[str] = Field(
        None, description="Dispersion correction method"
    )
    forcefield: EmptyStringToNone[str] = Field(
        None, description="Molecule dynamics force field"
    )
    solvent_model: EmptyStringToNone[str] = Field(None, description="Solvent model")

    @computed_field(description="Formatted summary of monomer-polymer states")
    @property
    def state_summary(self) -> str:
        get_state = lambda s: s if s else "x"
        return f"{get_state(self.monomer_state)}-{get_state(self.polymer_state)}"


class Thermo(BaseModel):
    delta_h: EmptyStringToNone[float] = (
        Field(None, description="Enthalpy of polymerisation (kJ/mol)"),
    )
    delta_s: EmptyStringToNone[float] = (
        Field(None, description="Entropy of polymerisation (J/mol)"),
    )
    ceiling_temperature: EmptyStringToNone[float] = (
        Field(None, description="Ceiling temperature in C"),
    )

    # compute ceiling temperature from H and S if present?
    # def model_post_init(self, _):
    #     if self.ceiling_temperature is None:
    #         if isinstance(self.delta_h, float) and isinstance(self.delta_s, float):
    #             self.ceiling_temperature = 1000 * self.delta_h / self.delta_s


class Metadata(BaseModel):
    # TODO: change year to be a datetime object?
    year: EmptyStringToNone[str] = (Field(..., description="Year of publication"),)
    comment: EmptyStringToNone[str] = (
        Field(..., description="Additional comments or notes"),
    )
    doi: EmptyStringToNone[str] = (
        Field(..., description="Digital Object Identifier for the source"),
    )
    url: EmptyStringToNone[str] = (Field(..., description="URL to related resource"),)


class Polymerisation(Document):

    polymerisation_id: Optional[int] = Field(
        None, description="unique display id for the polymerisation"
    )
    type: PolymerisationType = Field(..., description="Type of polymerisation")
    monomer: Monomer = Field(..., description="Monomer of the polymerisation")
    initiator: Initiator = Field(..., description="Initiator of the polymerisation")
    product: Product = Field(..., description="Product of the polymerisation")
    parameters: Parameters = Field(
        ..., description="Experimental/Computational parameters"
    )
    thermo: Thermo = Field(..., description="Thermodynamic data/results")
    metadata: Metadata = Field(
        ..., description="Polymerisation references and metadata"
    )

    class Settings:
        name = "polymerisations"


class DataRow(BaseModel):
    type: PolymerisationType = Field(..., description="Type of polymerisation")
    is_experimental: bool = Field(
        ...,
        description="Flag indicating whether the reaction is experimental (True) or computational (False)",
    )
    monomer_state: Optional[str] = Field(..., description="State of the monomer")
    polymer_state: Optional[str] = Field(..., description="State of the polymer")
    monomer_conc: Optional[float] = Field(
        ..., description="Initial concentration of the monomer"
    )
    solvent: Optional[str] = Field(
        ..., description="Solvent that the polymerisation is conducted within"
    )
    solvent_conc: Optional[float] = Field(
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
    year: Optional[str] = Field(..., description="Year of publication")


class MoleculeSummary(Document):
    polymerisations: list[Polymerisation] = Field(
        ..., description="list of polymerisations for the corresponding molecule"
    )
    data: list[DataRow] = Field(
        ...,
        description="table of data where each row corresponds to a polymerisation (for display purposes)",
    )

    @computed_field(description="Whether the molecule has experimental data")
    @property
    def has_experimental(self) -> bool:
        for poly in self.polymerisations:
            if poly.parameters.is_experimental:
                return True
        return False


class MonomerSummary(MoleculeSummary):
    monomer_id: int = Field(
        ..., description="unique display id for the monomer summary"
    )
    monomer: Monomer = Field(..., description="corresponding monomer")

    class Settings:
        name = "monomerSummaries"


class InitiatorSummary(MoleculeSummary):
    initiator_id: int = Field(
        ..., description="unique display id for the initiator summary"
    )
    initiator: Initiator = Field(..., description="corresponding initiator")

    class Settings:
        name = "initiatorSummaries"


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
