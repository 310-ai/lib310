from typing import Optional
from ._connection import DatabaseConnection, set_gcloud_key_path

db_connection = None


def fetch(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)

    query_results = db_connection.query(**kwargs)
    return query_results
