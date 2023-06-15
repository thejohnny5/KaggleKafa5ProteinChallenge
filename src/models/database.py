import sqlite3

# Create a connection to the database

connection = sqlite3.connect('models.db')
cursor = connection.cursor()

# Create the necessary tables
cursor.execute('''
    CREATE TABLE IF NOT EXISTS model (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        model_data BLOB,
        model_md5 TEXT,
        model_params TEXT,
        model_type TEXT,
        created_on DATETIME,
        CONSTRAINT unique_model_md5 UNIQUE (model_md5)

    )
''')

# Create additional tables or perform other necessary operations

# Commit the changes and close the connection
connection.commit()
connection.close()