import pandas as pd

from ._base import ProteinDataTable


def read_csv(*args, **kwargs):
    df = pd.read_csv(*args, **kwargs)
    return ProteinDataTable.from_pandas(df)
