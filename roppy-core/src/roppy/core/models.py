from datetime import datetime
from typing import Optional

from pydantic import BaseModel, Field
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdinchi import MolToInchi
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

from roppy.core.efgs import get_dec_fgs
from roppy.core.validate import (
    EmptyStringToNone,
    Method,
    PolymerisationType,
    Smiles,
    Solvent,
    State,
    get_ring_size,
    get_xyz,
)


class Molecule(BaseModel):
    smiles: Smiles = Field(
        description="SMILES notation representing the molecular structure",
    )

    inchi: Optional[str] = Field(
        description="Inchi identifier for the molecule",
        default_factory=lambda data: MolToInchi(MolFromSmiles(data["smiles"]))[0],
    )

    molecular_weight: Optional[float] = Field(
        description="Molecular weight of the molecule in g/mol",
        default_factory=lambda data: CalcExactMolWt(MolFromSmiles(data["smiles"])),
    )

    functional_groups: Optional[list[str]] = Field(
        description="Primary functional group of molecule",
        default_factory=lambda data: get_dec_fgs(MolFromSmiles(data["smiles"]))[2],
    )

    iupac_name: Optional[str] = Field(
        None,
        description="IUPAC name of the molecule",
    )

    common_name: Optional[str] = Field(
        None,
        description="Common name of the molecule",
    )


class Monomer(Molecule):
    ring_size: Optional[int] = Field(
        description="Size of the ring in the monomer structure",
        default_factory=lambda data: get_ring_size(data["smiles"]),
    )

    xyz: Optional[str] = Field(
        description="XYZ coordinates for the molecule",
        default_factory=lambda data: get_xyz(data["smiles"]),
    )


class Product(BaseModel):
    smiles: EmptyStringToNone[str] = Field(
        None,
        description="Big SMILES notation representing the polymer",
    )
    repeating_units: EmptyStringToNone[float] = Field(
        None,
        description="Number of repeating monomer units in the computational polymer chain/ring",
    )
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


class Parameters(BaseModel):
    is_experimental: bool = Field(
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
    initial_monomer_conc: EmptyStringToNone[float] = Field(
        None, description="Initial concentration of the monomer"
    )
    bulk_monomer_conc: EmptyStringToNone[float] = Field(
        None, description="Bulk concentration of the monomer"
    )
    solvent: Solvent = Field(
        None, description="Solvent that the polymerisation is conducted within"
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
    state_summary: str = Field(
        description="Formatted summary of monomer-polymer states",
        default_factory=lambda data: f"{data['monomer_state'] or 'x'}-{data['polymer_state'] or 'x'}",
    )


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


class Metadata(BaseModel):
    date: EmptyStringToNone[datetime] = Field(description="Date of publication")
    comment: EmptyStringToNone[str] = (
        Field(description="Additional comments or notes"),
    )
    doi: EmptyStringToNone[str] = (
        Field(description="Digital Object Identifier for the source"),
    )
    url: EmptyStringToNone[str] = (Field(description="URL to related resource"),)
    formatted_reference: EmptyStringToNone[str] = Field(
        None,
        description="Formatted reference string for the publication",
    )


class Polymerisation(BaseModel):
    polymerisation_id: str = Field(
        "poly-0", description="unique display id for the polymerisation"
    )
    type: PolymerisationType = Field(description="Type of polymerisation")
    monomer: Monomer = Field(description="Monomer of the polymerisation")
    initiator: Molecule = Field(description="Initiator of the polymerisation")
    product: Product = Field(description="Product of the polymerisation")
    parameters: Parameters = Field(description="Experimental/Computational parameters")
    thermo: Thermo = Field(..., description="Thermodynamic data/results")
    metadata: Metadata = Field(description="Polymerisation references and metadata")


class DataRow(BaseModel):
    type: PolymerisationType = Field(description="Type of polymerisation")
    is_experimental: bool = Field(
        description="Flag indicating whether the reaction is experimental (True) or computational (False)",
    )
    monomer_state: Optional[str] = Field(description="State of the monomer")
    polymer_state: Optional[str] = Field(description="State of the polymer")
    initial_monomer_conc: Optional[float] = Field(
        description="Initial concentration of the monomer"
    )
    bulk_monomer_conc: Optional[float] = Field(
        description="bulk concentration of the monomer"
    )
    solvent: Optional[str] = Field(
        description="Solvent that the polymerisation is conducted within"
    )
    delta_h: Optional[float] = Field(
        description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: Optional[float] = Field(
        description="Change in entropy (ΔS) for the polymerization reaction"
    )
    ceiling_temperature: Optional[float] = Field(description="Ceiling temperature in K")
    date: Optional[datetime] = Field(description="Year of publication")
    doi: Optional[str] = Field(description="DOI of publication")
    formatted_reference: Optional[str] = Field(description="Formatted reference string")


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
