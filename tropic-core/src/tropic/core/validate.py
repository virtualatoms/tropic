from typing import Annotated, Literal, TypeAlias

from pydantic import AfterValidator
from rdkit.Chem import AddHs, AllChem, MolFromSmiles, MolToXYZBlock
from rdkit.Chem.rdinchi import MolToInchi
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from rdkit.Contrib.efgs.efgs import get_dec_fgs


def validate_solvent(solvent: str) -> str:
    """Validates and normalizes solvent names."""
    if solvent in CANONICAL_SOLVENTS:
        return solvent
    if solvent in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[solvent]
    raise ValueError(f"Invalid solvent name: {solvent}. Must be one of {CANONICAL_SOLVENTS}.")

def get_mol(smiles: str) -> str | None:
    """Returns a molecule object from a SMILES string."""
    mol = MolFromSmiles(smiles)
    if mol is None:
        return ValueError("Invalid SMILES string")
    return mol


def get_ring_size(smiles: str) -> int | None:
    """Calculates the size of the largest ring in a molecule given its SMILES."""
    mol = get_mol(smiles)
    ring_info = mol.GetRingInfo()
    n_atoms = mol.GetNumAtoms()
    max_ring_size = max([ring_info.MinAtomRingSize(i) for i in range(n_atoms)])
    if max_ring_size == 0:
        return None
    return max_ring_size


def get_xyz(smiles: str) -> str:
    """Generates XYZ coordinates for a molecule given its SMILES."""
    mol = AddHs(get_mol(smiles))
    AllChem.EmbedMolecule(mol)
    return MolToXYZBlock(mol)

def get_inchi(smiles: str) -> str:
    """Generates InChI for a molecule given its SMILES."""
    return MolToInchi(get_mol(smiles))

def get_molecular_weight(smiles: str) -> float:
    """Calculates the molecular weight of a molecule given its SMILES."""
    return CalcExactMolWt(get_mol(smiles))


def get_func_groups(smiles: str) -> list[str]:
    """Extracts functional groups from a molecule given its SMILES."""
    return get_dec_fgs(get_mol(smiles))[2]

Smiles = Annotated[str, AfterValidator(StandardizeSmiles)]
Solvent = Annotated[str, AfterValidator(validate_solvent)]

# ROP: Ring-Opening Polymerisation (monomer -> polymer, no end groups)
# ROR: Ring-Opening Reaction (monomer -> polymer, with end groups)
# RER: Ring-Expansion Reaction (monomer -> cyclic chain)
# RCE: Ring-Chain Equilibrium (all rings -> polymer)
PolymerisationType = Literal["ROR", "RER", "ROP", "RCE"]

# g: gaseous state
# l: liquid state
# s: in solution
# c: crystalline or partially crystalline state
# a: condensed, glassy (amorphous) state
State: TypeAlias = Literal["g", "l", "s", "c", "a"]

# vant_hoff: Van't Hoff Curve
# DSC: Differential Scanning Calorimetry
Method: TypeAlias = Literal[
    "dft",
    "ffmd",
    "aimd",
    "mlmd",
    "xtb",
    "ml",
    "vant_hoff",
    "DSC",
    "NMR",
    "calorimetry",
]

# TODO: choose accepted values for computational fields
SOLVENT_MODELS = []
FUNCTIONALS = []
BASIS_SETS = []
DISPERSIONS = []
FORCEFIELDS = []

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