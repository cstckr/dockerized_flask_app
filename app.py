from extra.forms import LoginForm, ChangePasswordForm
from extra.users import User, Base
from extra.functions import smiles_to_base64_img
from flask import (Flask, render_template, request, redirect, url_for,
                   session, flash)
from flask_login import (login_user, LoginManager, login_required,
                         logout_user, current_user)
from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt
from datetime import timedelta


app = Flask(__name__)
db = SQLAlchemy(app)
flask_bcrypt = Bcrypt(app)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///database1.db"
app.config["SECRET_KEY"] = "notasafekey"
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False


login_manager = LoginManager()
login_manager.init_app(app)
login_manager.session_protection = "strong"
login_manager.login_view = "login"


@login_manager.user_loader
def load_user(user_id):
    return db.session.query(User).get(int(user_id))


Base.prepare(db.engine, reflect=True)


@app.before_first_request
def make_session_permanent():
    session.permanent = True
    app.permanent_session_lifetime = timedelta(seconds=60)


@app.route("/", methods=["GET", "POST"])
def login():
    form = LoginForm()

    user = db.session.query(User).filter_by(
        username=form.username.data).first()
    if user:
        if flask_bcrypt.check_password_hash(user.password, form.password.data):
            login_user(user)
            return redirect(url_for("main"))
        else:
            flash("Invalid username or password")
    return render_template("login.html", form=form)


@app.route("/logout", methods=["GET", "POST"])
@login_required
def logout():
    logout_user()
    return redirect(url_for("login"))


@app.route("/change_password", methods=["GET", "POST"])
@login_required
def change_password():
    form = ChangePasswordForm()
    if form.validate_on_submit():
        condition1 = flask_bcrypt.check_password_hash(current_user.password,
                                                      form.old_password.data)
        if condition1:
            new_password_hash = flask_bcrypt.generate_password_hash(
                form.new_password1.data).decode("utf-8")
            db.session.query(User).filter_by(
                username=current_user.username).update(
                    dict(password=new_password_hash))
            db.session.commit()
            return redirect(url_for("main"))
        else:
            flash("Old password is incorrect")
    return render_template("change_password.html", form=form)


@app.route("/main", methods=["GET", "POST"])
@login_required
def main():
    return render_template("main.html", username=current_user.username)


@app.route("/output", methods=["GET", "POST"])
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
        return redirect(url_for("main"))


if __name__ == "__main__":
    app.run(host="0.0.0.0")
