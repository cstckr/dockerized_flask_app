import sqlite3
from flask_bcrypt import Bcrypt


def create_database(database):
    try:
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        cursor.execute("""CREATE TABLE IF NOT EXISTS user
                        (user_id INTEGER PRIMARY KEY,
                        username TEXT UNIQUE,
                        password TEXT)""")
        connection.commit()
        connection.close()
    except Exception as error:
        print(error)


def add_user(database, username, password):
    try:
        flask_bcrypt = Bcrypt()
        pw_hash = flask_bcrypt.generate_password_hash(password).decode("utf-8")

        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        cursor.execute(f"""INSERT INTO user VALUES (NULL, "{username}", 
                       "{pw_hash}")""")
        # cursor.execute("""SELECT * FROM users """)
        # print(cursor.fetchall())
        connection.commit()
        connection.close()
    except Exception as error:
        print(error)


create_database("user_database1.db")
add_user("user_database1.db", "guest", "topsecret")
