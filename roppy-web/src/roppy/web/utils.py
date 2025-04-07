import base64
from rdkit import Chem
from rdkit.Chem import Draw, AllChem


def smiles_to_image(smiles, size=(150, 100)):
    drawer = Draw.MolDraw2DSVG(*size)
    dopts = drawer.drawOptions()
    dopts.bondLineWidth = 1.0  # default is 2.
    dopts.clearBackground = False  # default is True
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            svg_base64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
            return svg_base64
    except Exception as e:
        print(f"Error converting SMILES to image: {e}")
        return None
    return None

def reaction_to_image(smarts, size=(630, 275)):
    drawer = Draw.MolDraw2DSVG(*size)
    dopts = drawer.drawOptions()
    dopts.bondLineWidth = 1.5  # default is 2.
    dopts.padding = 0
    dopts.clearBackground = False  # default is True

    try:
        rxn = Chem.AllChem.ReactionFromSmarts(smarts, useSmiles=True)
        if rxn:
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            svg_base64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
            return svg_base64
    except Exception as e:
        print(f"Error converting reaction to image: {e}")
        return None
    return None
