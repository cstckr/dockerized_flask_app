from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, PasswordField
from wtforms.validators import InputRequired


class LoginForm(FlaskForm):
    username = StringField(validators=[InputRequired()],
                           render_kw={"placeholder": "Username"})
    password = PasswordField(validators=[InputRequired()],
                             render_kw={"placeholder": "Password"})
    submit = SubmitField("Login")


class ChangePasswordForm(FlaskForm):
    old_password = PasswordField(validators=[InputRequired()],
                                 render_kw={"placeholder": "old password"})
    new_password1 = PasswordField(validators=[InputRequired()],
                                  render_kw={"placeholder": "new password"})
    new_password2 = PasswordField(validators=[InputRequired()],
                                  render_kw={"placeholder":
                                             "repeat new password"})
    submit = SubmitField("Change")
