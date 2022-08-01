from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64


def smiles_to_base64_img(SMILES):
    mol = Chem.MolFromSmiles(SMILES)
    image = Draw.MolToImage(mol)
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    encoded_image = base64.b64encode(buffer.getvalue())
    return encoded_image.decode("utf-8")
