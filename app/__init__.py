from flask import Flask
from flask_login import LoginManager
from extra.models import User, Base
from extra.extensions import db
from blueprints.authentication import authentication
from blueprints.index import index


def create_app():
    app = Flask(__name__)

    app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///database1.db"
    app.config["SECRET_KEY"] = "notasafekey"
    app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

    db.app = app
    db.init_app(app)
    Base.prepare(db.engine, reflect=True)

    login_manager = LoginManager()
    login_manager.init_app(app)
    login_manager.session_protection = "strong"
    login_manager.login_view = "authentication.login"

    @login_manager.user_loader
    def load_user(user_id):
        return db.session.query(User).get(int(user_id))

    app.register_blueprint(authentication)
    app.register_blueprint(index)

    return app
