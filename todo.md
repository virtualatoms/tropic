# Input template
 - Change state fields to a (amorphous state, bulk polymer melt), g (glass), l (liquid monomer), s (solution), c (semicrystalline)
 - Include units for input form fields (e.g. kJmol-1 for enthalpy, Jmol-1K-1 for entropy)
 - Merge temperature fields into single ceiling temperature
 - Include fields for average molar mass and dispersity
 - Two different templates depending on whether data is computational/experimental


# Database
 - Remove field for product smiles
 - Merge temperature fields into single ceiling temperature
 - Include fields for average molar mass and dispersity
 - Move Monomer.molecular_weight and Monomer.functional_group into MonomerSummaryA
 - Change Monomer class into monomer_smiles (str)
 - Merge PolymerisationData and PolymerisationSummary into one Polymerisation class
 - Rename Results to Thermo
