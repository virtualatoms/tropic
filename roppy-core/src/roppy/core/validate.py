from typing import Optional, Annotated, TypeVar
from pydantic import BeforeValidator, AfterValidator
from roppy.core.constants import STATES, METHODS, POLY_TYPES
from roppy.core.solvents import CANONICAL_SOLVENTS, SOLVENT_ALIASES
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles
from rdkit.Chem import MolFromSmiles


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


def get_ring_size(smiles: str) -> Optional[int]:
    """Calculates the size of the largest ring in a molecule given its SMILES."""
    try:
        mol = MolFromSmiles(smiles)
        ring_info = mol.GetRingInfo()
        n_atoms = mol.GetNumAtoms()
        max_ring_size = max([ring_info.MinAtomRingSize(i) for i in range(n_atoms)])
        if max_ring_size == 0:
            return None
        return max_ring_size
    except:
        return None


EmptyStringToNone = Annotated[
    Optional[TypeVar("T")],
    BeforeValidator(lambda s: s if s else None),
]
Smiles = Annotated[Optional[str], AfterValidator(validate_smiles)]
Solvent = Annotated[Optional[str], BeforeValidator(validate_solvent)]
PolymerisationType = Annotated[
    Optional[str], BeforeValidator(lambda s: s if s in POLY_TYPES else None)
]
State = Annotated[Optional[str], BeforeValidator(lambda s: s if s in STATES else None)]
Method = Annotated[
    Optional[str], BeforeValidator(lambda s: s if s in METHODS else None)
]
