from flask import Blueprint, render_template, request, redirect, url_for
from flask_login import login_required, current_user
from extra.functions import smiles_to_base64_img

index = Blueprint("index", __name__)


@index.route("/main", methods=["GET", "POST"])
@login_required
def main():
    return render_template("main.html", username=current_user.username)


@index.route("/output", methods=["GET", "POST"])
@login_required
def output():
    name = request.form.get("name")
    if name:
        try:
            return render_template("output.html", name=name,
                                   image=smiles_to_base64_img(name))
        except ValueError:
            return render_template("invalid_input.html", name=name)
    else:
        return redirect(url_for("index.main"))
