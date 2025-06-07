# ROP: Ring-Opening Polymerisation (monomer -> polymer, no end groups)
# ROR: Ring-Opening Reaction (monomer -> polymer, with end groups)
# RER: Ring-Expansion Reaction (monomer -> cyclic chain)
# RCE: Ring-Chain Equilibrium (all rings -> polymer)
POLY_TYPES = ["ROR", "RER", "ROP", "RCE"]

# g: gaseous state
# l: liquid state
# s: in solution
# c: crystalline or partially crystalline state
# a: condensed, glassy (amorphous) state
STATES = ["g", "l", "s", "c", "a"]

# vant_hoff: Van't Hoff Curve
# DSC: Differential Scanning Calorimetry
METHODS = [
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
