from sqlalchemy import Column, Integer, String, JSON
from sqlalchemy.ext.declarative import declarative_base
import json

Base = declarative_base()


class InteractionBaseModel(Base):
    __abstract__ = True

    row_id = Column(Integer, primary_key=True)
    uniparc_id = Column(String)
    interactions = Column(JSON)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.row_id = kwargs.get('row_id')
        self.uniparc_id = kwargs.get('uniparc_id')
        self.interactions = kwargs.get('interactions')

    def to_dict(self):
        return {
            'row_id': self.row_id,
            'interactions': self.interactions
        }

    def __repr__(self):
        return f"<InteractionModel(row_id='{self.row_id}', uniparc_id='{self.uniparc_id}', interactions={self.interactions})>"


class InteractionModel(InteractionBaseModel):
    __tablename__ = 'interaction'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class Interaction300T10Model(InteractionBaseModel):
    __tablename__ = 'interaction_300_t10'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class Interaction300Model(InteractionBaseModel):
    __tablename__ = 'interaction_300'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class InteractionT10Model(InteractionBaseModel):
    __tablename__ = 'interaction_t10'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)