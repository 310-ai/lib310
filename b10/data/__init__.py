from typing import Optional
from ._connection import DatabaseConnection

db_connection = None


def fetch(limit: Optional[int] = 1000, **query_params):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection()

    query_results = db_connection.query(**query_params, limit=limit)
    return query_results
