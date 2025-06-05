# ROP: Ring-Opening Polymerisation (monomer -> polymer)
# RCE: Ring-Chain Equilibrium (all rings -> polymer)
POLY_TYPES = ["ROR", "RER", "ROP", "RCE"]

# g: gas
# l: liquid
# s: solution
# c: crystalline solid
# a: amorphous solid
STATES = ["g", "l", "s", "c", "a"]

# vant_hof: Van't Hoff Curve
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
    "calorimetry",
]

# TODO: choose accepted values for computational fields
SOLVENT_MODELS = []
FUNCTIONALS = []
BASIS_SETS = []
DISPERSIONS = []
FORCEFIELDS = []
