import pickle
import hashlib
from datetime import datetime
from sqlalchemy import create_engine, Column, Integer, String, DateTime, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import pandas as pd
Base = declarative_base()

class Model(Base):
    __tablename__ = 'model'

    id = Column(Integer, primary_key=True)
    model_data = Column(LargeBinary)
    model_md5 = Column(String)
    model_params = Column(String)
    model_type = Column(String)
    created_on = Column(DateTime)

class InvalidModel(Exception):
    pass

class ModelHandler:
    def __init__(self) -> None:
        self.model: object = None
        self.classname: str = None
        self.name: str = None
        self.model_params: str = None

    def set_model(self, model: object, name: str) -> None:
        self.classname = model.__class__.__name__
        self.model = model
        self.name = name
        self.model_params = self._get_params()

    def _save(self) -> bytes:
        """Saves model as pickled object"""
        return pickle.dumps(self.model)

    def _get_params(self) -> dict:
        try:
            return self.model.get_params()
        except InvalidModel:
            raise InvalidModel("Model does not support get_params() method.")

    def _hash(self) -> str:
        return hashlib.md5(self._save()).hexdigest()

    def write_to_db(self, db: str) -> None:
        if self.model is None:
            raise Exception("Must Write Model First")

        model_data = self._save()
        model_md5 = self._hash()
        model_params = str(self.model_params)
        model_type = self.classname
        created_on = datetime.now()

        engine = create_engine('sqlite:///' + db)  # SQLite connection string
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)
        session = Session()

        model_entry = Model(
            model_data=model_data,
            model_md5=model_md5,
            model_params=model_params,
            model_type=model_type,
            created_on=created_on
        )

        session.add(model_entry)
        session.commit()
        session.close()

    def load_model(self, db: str, id: int) -> object:
        engine = create_engine('sqlite:///' + db)  # SQLite connection string
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)
        session = Session()

        model_entry = session.query(Model).filter_by(id=id).first()
        if model_entry:
            loaded_model = pickle.loads(model_entry.model_data)
        else:
            loaded_model = None

        session.close()
        return loaded_model
    
    def load_top_n_rows(self, db: str, n: int) -> pd.DataFrame:
        engine = create_engine('sqlite:///' + db)  # SQLite connection string
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)
        session = Session()

        query = session.query(Model).limit(n)
        df = pd.read_sql(query.statement, query.session.bind)

        session.close()
        return df