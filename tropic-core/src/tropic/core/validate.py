"""Validation functions for chemical properties and structures."""

from typing import Annotated, Literal, TypeAlias

from pydantic import AfterValidator
from rdkit.Chem import AddHs, AllChem, MolFromSmiles, MolToXYZBlock
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles
from rdkit.Chem.rdchem import Atom, Mol, RingInfo
from rdkit.Chem.rdinchi import MolToInchi
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from rdkit.Chem.rdmolfiles import MolFromSmarts

MONOMER_GROUPS: dict[str, str] = {
    "CC": "[O;R]-[C;R](=[O;!R])-[O;R]",
    "CtC": "[O;R]-[C;R](=[O;!R])-[S;R]",
    "CdtC": "[S;R]-[C;R](=[O;!R])-[S;R]",
    "CtnC": "[O;R]-[C;R](=[S;!R])-[O;R]",
    "CX": "[O;R]-[C;R](=[S;!R])-[S;R]",
    "CtX": "[S;R]-[C;R](=[S;!R])-[S;R]",
    "L": "[C,c;R]-[C;R](=[O;!R])-[O;R]",
    "tL": "[C,c;R]-[C;R](=[O;!R])-[S;R]",
    "tnL": "[C,c;R]-[C;R](=[S;!R])-[O;R]",
    "dtL": "[C,c;R]-[C;R](=[S;!R])-[S;R]",
    "Lm": "[C,c;R]-[C;R](=[O;!R])-[N;R]",
    "oA": "[O;R]-[C;R](=[O;!R])-[N;R]",
}


def validate_solvent(solvent: str) -> str:
    """Validate and normalize solvent names."""
    if solvent in CANONICAL_SOLVENTS:
        return solvent
    if solvent in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[solvent]
    return solvent
    # For now, we return the solvent as is.
    # raise ValueError(f"Invalid solvent name: {solvent}.")


def get_mol(smiles: str) -> Mol:
    """Get an rdkit molecule object from a SMILES string."""
    mol: Mol | None = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Failed to create RDKit Mol from SMILES string")
    return mol


def get_func_group(smiles: str) -> str:
    """Extract ROP functional group for monomer."""
    mol: Mol = get_mol(smiles)
    functional_groups: dict[str, Mol] = {
        fg_name: MolFromSmarts(fg_smarts)
        for fg_name, fg_smarts in MONOMER_GROUPS.items()
    }
    for fg_name, fg_mol in functional_groups.items():
        if mol.HasSubstructMatch(fg_mol):
            return fg_name
    return "other"


def get_ring_size(smiles: str) -> int | None:
    """Get the size of the functional group ring in the monomer."""
    mol: Mol = get_mol(smiles)
    ring_info: RingInfo = mol.GetRingInfo()
    functional_groups: list[Mol] = [
        MolFromSmarts(fg_smarts) for fg_smarts in MONOMER_GROUPS.values()
    ]
    for fg_mol in functional_groups:
        if mol.HasSubstructMatch(fg_mol):
            substruct_atom_idx: tuple[int] = mol.GetSubstructMatch(fg_mol)
            for atom_id in substruct_atom_idx:
                substruct_atom: Atom = mol.GetAtomWithIdx(atom_id)
                if (substruct_atom.GetAtomicNum() != 6) and (substruct_atom.IsInRing()):
                    return ring_info.MinAtomRingSize(atom_id)
    return None


def get_xyz(smiles: str) -> str:
    """Generate XYZ coordinates for a molecule."""
    mol = AddHs(get_mol(smiles))
    AllChem.EmbedMolecule(mol)
    return MolToXYZBlock(mol)


def get_inchi(smiles: str) -> str:
    """Get InChI key for a molecule."""
    return MolToInchi(get_mol(smiles))[0]


def get_molecular_weight(smiles: str) -> float:
    """Calculate the molecular weight of a molecule."""
    return CalcExactMolWt(get_mol(smiles))


Smiles = Annotated[str, AfterValidator(StandardizeSmiles)]
Medium = Annotated[str, AfterValidator(validate_solvent)]

# ROP: Ring-Opening Reaction (monomer -> polymer, no end groups)
# ROR: Ring-Opening Reaction (monomer -> polymer, with end groups)
# RER: Ring-Expansion Reaction (monomer -> cyclic chain)
# RCE: Ring-Chain Equilibrium (all rings -> polymer)
ReactionType = Literal["ROR", "RER", "ROP", "RCE"]

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
# SOLVENT_MODELS = []
# FUNCTIONALS = []
# BASIS_SETS = []
# DISPERSIONS = []
# FORCEFIELDS = []

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
