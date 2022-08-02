from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt
# from flask_login import UserMixin
# from sqlalchemy.ext.automap import automap_base


db = SQLAlchemy()
flask_bcrypt = Bcrypt()

# Base = automap_base()


# class User(Base, UserMixin):
#     __tablename__ = 'user'

#     def get_id(self):
#         return self.user_id
