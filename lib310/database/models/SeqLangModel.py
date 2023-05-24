from typing import Any
from sqlalchemy import Column, Integer, String, Float, ARRAY
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class SeqLangBaseModel(Base):
    __abstract__ = True

    row_id = Column(Integer, primary_key=True)
    uniparc_id = Column(String)
    sequence = Column(String)
    len = Column(Integer)
    stage = Column(String)
    token_size = Column(Integer)
    token_ids = Column(ARRAY(Integer))
    token_names = Column(String)
    protein_names = Column(String)
    taxonomies = Column(String)
    organisms = Column(String)
    rand = Column(Float)

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.row_id = kwargs.get('row_id')
        self.uniparc_id = kwargs.get('uniparc_id')
        self.sequence = kwargs.get('sequence')
        self.len = kwargs.get('len')
        self.stage = kwargs.get('stage')
        self.token_size = kwargs.get('token_size')
        self.token_ids = kwargs.get('token_ids')
        self.token_names = kwargs.get('token_names')
        self.protein_names = kwargs.get('protein_names')
        self.taxonomies = kwargs.get('taxonomies')
        self.organisms = kwargs.get('organisms')
        self.rand = kwargs.get('rand')

    def to_dict(self) -> dict:
        d = {
            'row_id': self.row_id,
            'uniparc_id': self.uniparc_id,
            'sequence': self.sequence,
            'stage': self.stage,
            'token_ids': self.token_ids,
            'token_names': self.token_names,
            'protein_names': self.protein_names,
            'taxonomies': self.taxonomies,
            'organisms': self.organisms,
        }
        return d

    def __repr__(self) -> str:
        return f"<SeqLangModel(uniparc_id='{self.uniparc_id}', sequence='{self.sequence}', len={self.len}, stage='{self.stage}', token_size={self.token_size}, token_ids={self.token_ids}, token_names='{self.token_names}', protein_names='{self.protein_names}', taxonomies='{self.taxonomies}', organisms='{self.organisms}', rand={self.rand})>"


class SeqLangModel(SeqLangBaseModel):
    __tablename__ = 'seq_lang'

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)


class SeqLang300T10Model(SeqLangBaseModel):
    __tablename__ = 'seq_lang_300_t10'

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)


class SeqLang300Model(SeqLangBaseModel):
    __tablename__ = 'seq_lang_300'

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)


class SeqLangT10Model(SeqLangBaseModel):
    __tablename__ = 'seq_lang_t10'

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)