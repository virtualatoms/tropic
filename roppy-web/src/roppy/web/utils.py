import base64
from rdkit import Chem
from rdkit.Chem import Draw


def smiles_to_image(smiles, size=(150, 100)):
    drawer = Draw.MolDraw2DSVG(*size)
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            svg = svg.replace("<rect style='opacity:1.0", "<rect style='opacity: 0")
            svg = svg.replace("stroke-width:2.0px", f"stroke-width:1.2px")
            svg_base64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
            return svg_base64
    except Exception as e:
        print(f"Error converting SMILES to image: {e}")
        return None
    return None
