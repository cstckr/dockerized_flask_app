import sqlite3

connection = sqlite3.connect('user_database.db')
cursor = connection.cursor()

cursor.execute("""INSERT INTO users VALUES (NULL, "guest", "topsecret")""")

# cursor.execute("""SELECT * FROM users """)
# print(cursor.fetchall())

# cursor.execute("""DROP TABLE IF EXISTS users """)
connection.commit()
connection.close()
