from pydantic import BaseModel, Field, field_validator
from typing import Optional, Union, Literal, Annotated
from beanie import Document, Link

# Dictionary mapping aliases to canonical names
SOLVENT_ALIASES = {
    "dioxane": "1,4-dioxane",
    "bromooctane": "1-bromooctane",
    "butanol": "1-butanol",
    "chlorohexane": "1-chlorohexane",
    "decanol": "1-decanol",
    "heptanol": "1-heptanol",
    "hexanol": "1-hexanol",
    "nonanol": "1-nonanol",
    "octanol": "1-octanol",
    "pentanol": "1-pentanol",
    "propanol": "1-propanol",
    "secbutanol": "2-butanol",
    "isobutanol": "2-methyl-1-propanol",
    "isopropanol": "2-propanol",
    "aceticacid": "acetic acid",
    "mecn": "acetonitrile",
    "ch3cn": "acetonitrile",
    "benzylalcohol": "benzyl alcohol",
    "butylacetate": "butyl ethanoate",
    "butylbenzene": "n-butylbenzene",
    "secbutylbenzene": "sec-butylbenzene",
    "tbutylbenzene": "tert-butylbenzene",
    "carbondisulfide": "carbon disulfide",
    "cs2": "carbon disulfide",
    "ccl4": "carbon tetrachloride",
    "chcl3": "chloroform",
    "mcresol": "m-cresol",
    "decane": "n-decane",
    "odichlorobenzene": "o-dichlorobenzene",
    "ch2cl2": "dichloromethane",
    "dcm": "dichloromethane",
    "diethylether": "diethyl ether",
    "diisopropylether": "diisopropyl ether",
    "dimethylacetamide": "n,n-dimethylacetamide",
    "dimethylformamide": "n,n-dimethylformamide",
    "dmf": "n,n-dimethylformamide",
    "dmso": "dimethylsulfoxide",
    "dodecane": "n-dodecane",
    "ethylacetate": "ethyl acetate",
    "ethanoate": "ethyl acetate",
    "ethoxybenzene": "ethyl phenyl ether",
    "heptane": "n-heptane",
    "hexadecane": "n-hexadecane",
    "hexane": "n-hexane",
    "isopropyltoluene": "p-isopropyltoluene",
    "methylformamide": "n-methylformamide",
    "phno2": "nitrobenzene",
    "meno2": "nitromethane",
    "onitrotoluene": "o-nitrotoluene",
    "nonane": "n-nonane",
    "octane": "n-octane",
    "pentadecane": "n-pentadecane",
    "wetoctanol": "octanol(wet)",
    "woctanol": "octanol(wet)",
    "pentane": "n-pentane",
    "hexafluorobenzene": "perfluorobenzene",
    "c2cl4": "tetrachloroethene",
    "thf": "tetrahydrofuran",
    "sulfolane": "tetrahydrothiophene-s,s-dioxide",
    "undecane": "n-undecane",
    "h2o": "water",
}

# Set of canonical solvent names
CANONICAL_SOLVENTS = {
    "1,1,1-trichloroethane",
    "1,1,2-trichloroethane",
    "1,2,4-trimethylbenzene",
    "1,2-dibromoethane",
    "1,2-dichloroethane",
    "1,2-ethanediol",
    "1,4-dioxane",
    "1-bromo-2-methylpropane",
    "1-bromooctane",
    "1-bromopentane",
    "1-bromopropane",
    "1-butanol",
    "1-chlorohexane",
    "1-chloropentane",
    "1-chloropropane",
    "1-decanol",
    "1-fluorooctane",
    "1-heptanol",
    "1-hexanol",
    "1-hexene",
    "1-hexyne",
    "1-iodobutane",
    "1-iodohexadecane",
    "1-iodopentane",
    "1-iodopropane",
    "1-nitropropane",
    "1-nonanol",
    "1-octanol",
    "1-pentanol",
    "1-pentene",
    "1-propanol",
    "2,2,2-trifluoroethanol",
    "2,2,4-trimethylpentane",
    "2,4-dimethylpentane",
    "2,4-dimethylpyridine",
    "2,6-dimethylpyridine",
    "2-bromopropane",
    "2-butanol",
    "2-chlorobutane",
    "2-heptanone",
    "2-hexanone",
    "2-methoxyethanol",
    "2-methyl-1-propanol",
    "2-methyl-2-propanol",
    "2-methylpentane",
    "2-methylpyridine",
    "2-nitropropane",
    "2-octanone",
    "2-pentanone",
    "2-propanol",
    "2-propen-1-ol",
    "e-2-pentene",
    "3-methylpyridine",
    "3-pentanone",
    "4-heptanone",
    "4-methyl-2-pentanone",
    "4-methylpyridine",
    "5-nonanone",
    "acetic acid",
    "acetone",
    "acetonitrile",
    "acetophenone",
    "ammonia",
    "aniline",
    "anisole",
    "benzaldehyde",
    "benzene",
    "benzonitrile",
    "benzyl alcohol",
    "bromobenzene",
    "bromoethane",
    "bromoform",
    "butanal",
    "butanoic acid",
    "butanone",
    "butanonitrile",
    "butyl ethanoate",
    "butylamine",
    "n-butylbenzene",
    "sec-butylbenzene",
    "tert-butylbenzene",
    "carbon disulfide",
    "carbon tetrachloride",
    "chlorobenzene",
    "chloroform",
    "a-chlorotoluene",
    "o-chlorotoluene",
    "m-cresol",
    "o-cresol",
    "cyclohexane",
    "cyclohexanone",
    "cyclopentane",
    "cyclopentanol",
    "cyclopentanone",
    "decalin",
    "cis-decalin",
    "n-decane",
    "dibromomethane",
    "dibutylether",
    "o-dichlorobenzene",
    "e-1,2-dichloroethene",
    "z-1,2-dichloroethene",
    "dichloromethane",
    "diethyl ether",
    "diethyl sulfide",
    "diethylamine",
    "diiodomethane",
    "diisopropyl ether",
    "cis-1,2-dimethylcyclohexane",
    "dimethyl disulfide",
    "n,n-dimethylacetamide",
    "n,n-dimethylformamide",
    "dimethylsulfoxide",
    "diphenylether",
    "dipropylamine",
    "n-dodecane",
    "ethanethiol",
    "ethanol",
    "ethyl acetate",
    "ethyl methanoate",
    "ethyl phenyl ether",
    "ethylbenzene",
    "fluorobenzene",
    "formamide",
    "formic acid",
    "furan",
    "n-heptane",
    "n-hexadecane",
    "n-hexane",
    "hexanoic acid",
    "iodobenzene",
    "iodoethane",
    "iodomethane",
    "isopropylbenzene",
    "p-isopropyltoluene",
    "mesitylene",
    "methanol",
    "methyl benzoate",
    "methyl butanoate",
    "methyl ethanoate",
    "methyl methanoate",
    "methyl propanoate",
    "n-methylaniline",
    "methylcyclohexane",
    "n-methylformamide",
    "nitrobenzene",
    "nitroethane",
    "nitromethane",
    "o-nitrotoluene",
    "n-nonane",
    "n-octane",
    "n-pentadecane",
    "octanol(wet)",
    "pentanal",
    "n-pentane",
    "pentanoic acid",
    "pentyl ethanoate",
    "pentylamine",
    "perfluorobenzene",
    "phenol",
    "propanal",
    "propanoic acid",
    "propanonitrile",
    "propyl ethanoate",
    "propylamine",
    "pyridine",
    "tetrachloroethene",
    "tetrahydrofuran",
    "tetrahydrothiophene-s,s-dioxide",
    "tetralin",
    "thiophene",
    "thiophenol",
    "toluene",
    "trans-decalin",
    "tributylphosphate",
    "trichloroethene",
    "triethylamine",
    "n-undecane",
    "water",
    "xylene",
    "m-xylene",
    "o-xylene",
    "p-xylene",
}

# Combined type for all possible solvent inputs
SolventInput = Annotated[str, Field(description="Solvent name (canonical or alias)")]


def validate_solvent(value: Optional[str]) -> Optional[str]:
    """Validates and normalizes solvent names."""
    if value is None:
        return None

    # Check if it's already a canonical name
    if value in CANONICAL_SOLVENTS:
        return value

    # Check if it's an alias
    if value in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[value]

    raise ValueError(
        f"Invalid solvent: {value}. Must be a valid ORCA solvent name or alias."
    )


class Conditions(Document):
    temperature: Optional[float] = Field(
        None, description="Reaction temperature in Kelvin"
    )
    solvent: Optional[SolventInput] = Field(
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


class Parameters(Document):
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
    solvent: Optional[SolventInput] = Field(
        None, description="Solvent used in calculations (must be ORCA-compatible)"
    )
    forcefield: Optional[str] = Field(
        None, description="Force field used for molecular dynamics"
    )

    _validate_solvent = field_validator("solvent")(validate_solvent)


class Monomer(Document):
    smiles: Optional[str] = Field(
        None, description="SMILES notation representing the monomer structure"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the monomer in g/mol"
    )
    functional_group: Optional[str] = Field(
        None, description="Primary functional group involved in polymerization"
    )


class Product(Document):
    repeating_unit: Optional[str] = Field(
        None, description="SMILES notation or chemical formula of the repeating unit"
    )
    number_of_units: Optional[Union[int, Literal["many"]]] = Field(
        None, description="Number of repeating units in the polymer"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the polymer in g/mol"
    )


class Initiator(Document):
    smiles: Optional[str] = Field(
        None, description="SMILES notation representing the initiator structure"
    )
    molecular_weight: Optional[float] = Field(
        None, description="Molecular weight of the initiator in g/mol"
    )
    functional_group: Optional[str] = Field(
        None, description="Primary functional group of the initiator"
    )


class Thermo(Document):
    delta_h: Optional[float] = Field(
        None, description="Change in enthalpy (ΔH) for the polymerization reaction"
    )
    delta_s: Optional[float] = Field(
        None, description="Change in entropy (ΔS) for the polymerization reaction"
    )


class Polymerization(Document):
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

    class Settings:
        name = "polymerizations"
        indexes = [
            "polymerisation_id",
            "monomer.smiles",
        ]


class MonomerSummary(Document):
    monomer_id: str = Field(..., description="Unique identifier for the monomer")
    smiles: str = Field(
        ..., description="SMILES notation representing the monomer structure"
    )
    ring_size: Union[int, Literal["many"]] = Field(
        ...,
        description="Size of the ring in the monomer structure, 'many' for complex structures",
    )
    average_delta_h: float = Field(
        ..., description="Average enthalpy change across all polymerization reactions"
    )
    average_delta_s: float = Field(
        ..., description="Average entropy change across all polymerization reactions"
    )
    polymerisation_id: list[Link[Polymerization]] = Field(
        default_factory=list,
        description="List of IDs referencing related polymerization reactions",
    )

    class Settings:
        name = "monomers_summary"
        indexes = [
            "monomer_id",
            "smiles",
        ]
