import json
import pandas as pd
import sqlite3
import joblib
from datetime import datetime
"""
Connect to database
Check if entry exists
Save Model using joblib
Save model params as metadata
"""
class InvalidModel(Exception): pass

class DBHandler:
    def __init__(self, db_local: str) -> None:
        self.db_local = db_local
        self.db_conn = None
        self.db_cursor = None

    def __enter__(self):
        # This ensure, whenever an object is created using "with"
        # this magic method is called, where you can create the connection.
        self.db_conn = sqlite3.connect(**self.db_local)
        self.db_cursor = self.db_conn.cursor()
        return self

    def __exit__(self, exception_type, exception_val, trace) -> None:
        # once the with block is over, the __exit__ method would be called
        # with that, you close the connnection
        try:
           self.db_cursor.close()
           self.db_conn.close()
        except AttributeError: # isn't closable
           AttributeError("Not Closable")
           return True # exception handled successfully

    def get_row(self, sql: str, data = None) -> dict:
        self.db_cursor.execute(sql)
        self.resultset = self.db_cursor.fetchall()
        return self.resultset
    
    def write_row(self, sql: str, data: tuple) -> str:
        self.db_cursor.execute(sql, data)
        self.db_conn.commit()
        return self.db_cursor.lastrowid


class ModelHandler:

    def __init__(self, model: object, name: str) -> None:
        self.model = model
        self.classname = model.__class__.__name__
        self.name = name
        self._save()
        self.model_params = self._get_params()
    
    def _save(self) -> None:
        """Saves model as joblib serialized model"""
        joblib.dump(self.model, self.name + ".joblib")
    def _get_params(self) -> dict:
        try:
            return self.model.get_params()
        except InvalidModel:
            InvalidModel("Model does not support get_params() method.")
    def write_to_db(self, db: str) -> None:
        sql_line = """
        INSERT INTO model(model_data, model_md5, model_params, model_type, created_on)
        VALUES(?,?,?,?,?)
        """
        line = ()

        with DBHandler(db) as f:
            f.write_row()
    
#joblib.dump(model, name + ".joblib")
#parameters = model.get_params()
#time = datetime.now()
