from flask import Flask, render_template, request, redirect, url_for
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
from flask_login import UserMixin, login_user, LoginManager, login_required, logout_user, current_user
from flask_sqlalchemy import SQLAlchemy
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, PasswordField
from wtforms.validators import InputRequired
from flask_bcrypt import Bcrypt
from sqlalchemy.ext.automap import automap_base


def smiles_to_base64_img(SMILES):
    mol = Chem.MolFromSmiles(SMILES)
    image = Draw.MolToImage(mol)
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    encoded_image = base64.b64encode(buffer.getvalue())
    return encoded_image.decode("utf-8")


app = Flask(__name__)
db = SQLAlchemy(app)
flask_bcrypt = Bcrypt(app)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///user_database1.db"
app.config["SECRET_KEY"] = "notasafekey"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False


Base = automap_base()
Base.prepare(db.engine, reflect=True)

User = Base.classes.user


class LoginForm(FlaskForm):
    username = StringField(validators=[InputRequired()],
                           render_kw={"placeholder": "Username"})
    password = PasswordField(validators=[InputRequired()],
                             render_kw={"placeholder": "Password"})
    submit = SubmitField('Login')


@app.route("/", methods=["GET", "POST"])
def login():
    form = LoginForm()
    user = db.session.query(User).filter_by(
        username=form.username.data).first()
    if user:
        if flask_bcrypt.check_password_hash(user.password, form.password.data):
            login_user(user)
            return redirect(url_for('main'))
    return render_template("login.html", form=form)


@app.route("/main")
# @login_required
def main():
    return render_template("main.html")


@app.route("/output", methods=["POST"])
# @login_required
def output():
    name = request.form.get("name")
    if name:
        try:
            return render_template("output.html", name=name,
                                   image=smiles_to_base64_img(name))
        except ValueError:
            return render_template("invalid_input.html", name=name)
    else:
        return redirect(url_for("main"))


if __name__ == "__main__":
    app.run(host="0.0.0.0")
    # app.run(host="0.0.0.0")
