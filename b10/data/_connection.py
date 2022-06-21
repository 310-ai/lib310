import os
from typing import List

from .client import Client

import numpy as np


def set_gcloud_key_path(path: str):
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = path


class DatabaseConnection(object):
    def __new__(cls, *args, **kwargs):
        if "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ:
            raise Exception('You must set BigQuery credentials first! try using `b10.dataAPI.set_db_key(path) method.`')
        else:
            return cls()

    def __init__(self):
        self.client = Client()
        print(f'Successfully established a connection with b10 Database!')

    def query(self, **query_params):
        table_name = query_params.get('table_name')
        query_string = f"SELECT * FROM `{table_name}` WHERE Entry IN {str(tuple(query_params.get('name', [])))}"
        query_results = self.client.query(query_string).result().to_dataframe()

        return query_results
