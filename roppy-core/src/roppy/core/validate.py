from typing import Optional, Annotated, TypeVar
from pydantic import BeforeValidator, AfterValidator
from roppy.core.constants import STATES, COMP_METHODS, POLY_TYPES
from roppy.core.solvents import CANONICAL_SOLVENTS, SOLVENT_ALIASES
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
Smiles = Annotated[str, AfterValidator(validate_smiles)]
Solvent = Annotated[Optional[str], BeforeValidator(validate_solvent)]
PolymerisationType = Annotated[
    Optional[str], BeforeValidator(lambda s: s if s in POLY_TYPES else None)
]
State = Annotated[Optional[str], BeforeValidator(lambda s: s if s in STATES else None)]
CompMethod = Annotated[
    Optional[str], BeforeValidator(lambda s: s if s in COMP_METHODS else None)
]
