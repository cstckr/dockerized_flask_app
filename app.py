from flask import Flask, render_template, request, redirect, url_for
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


app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/output", methods=["POST"])
def output():
    name = request.form.get("name")
    if name:
        try:
            return render_template("output.html", name=name,
                                   image=smiles_to_base64_img(name))
        except ValueError:
            return render_template("invalid_input.html", name=name)
    else:
        return redirect(url_for("index"))


if __name__ == "__main__":
    app.run(host="0.0.0.0")
