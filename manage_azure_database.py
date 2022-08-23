import pyodbc
from credentials.database import server, database, username, password
from flask_bcrypt import Bcrypt

driver = "{ODBC Driver 17 for SQL Server}"

connection_string = "DRIVER=" + driver + ";SERVER=tcp:" + server + \
    ";PORT=1433;DATABASE=" + database + ";UID=" + username + ";PWD=" + password


table_name = "users"


def sql_query(sql_statement):
    with pyodbc.connect(connection_string) as connection:
        with connection.cursor() as cursor:
            cursor.execute(sql_statement)
            cursor.commit()


def create_sql_str_insert_user(username, password):
    flask_bcrypt = Bcrypt()
    pw_hash = flask_bcrypt.generate_password_hash(password).decode("utf-8")
    sql_insert_user = f"""
        INSERT INTO {table_name}(username, password)
        VALUES ('{username}', '{pw_hash}')
        """
    return sql_insert_user


sql_str_drop_table = f"DROP TABLE IF EXISTS {table_name}"

sql_str_create_table = f"""
    IF OBJECT_ID(N'{table_name}', N'U') IS NULL
    BEGIN
    CREATE TABLE {table_name}
        (
        user_id INT IDENTITY PRIMARY KEY,
        username NVARCHAR(128) NOT NULL UNIQUE,
        password NVARCHAR(75) NOT NULL,
        )
    END;
    """


sql_query(sql_str_drop_table)
sql_query(sql_str_create_table)
sql_query(create_sql_str_insert_user("guest3", "123456"))
