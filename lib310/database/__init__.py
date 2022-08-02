from typing import Optional
from ._connection import DatabaseConnection, set_gcloud_key_path
from ._visualize import *

db_connection: DatabaseConnection = None
db_info = 'sandbox.info'


def fetch(*args, **kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)

    query_results = db_connection.query(*args, **kwargs)
    return query_results


def get_table_info(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)
    elif db_connection.table_name != kwargs.get('table_name', db_connection.table_name):
        db_connection = DatabaseConnection(**kwargs)

    return db_connection.get_table_info()


def list_datasets(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)
    return list(db_connection.client.list_datasets())


def list_tables(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)

    dataset_ids = kwargs.get('dataset_ids')

    if dataset_ids is None:
        dataset_ids = [d.dataset_id for d in list(db_connection.client.list_datasets())]

    if not isinstance(dataset_ids, list):
        dataset_ids = [dataset_ids]
        return list(db_connection.client.list_tables(dataset_ids[0]))

    result = {}
    for dataset_id in dataset_ids:
        result[dataset_id] = list(db_connection.client.list_tables(dataset_id))
    return result


def summary(**kwargs):
    global db_connection
    db_connector(**kwargs)
    client = db_connection.client
    df = client.query(f"SELECT dataset_name, table_name, num_rows, num_cols, size FROM `{db_info}`").to_dataframe()

    print_it = kwargs.get('print')
    if print_it is None or print_it is True:
        df['size'] = df['size'].apply(size_to_abbrevation)
        print(df.to_string(index = False))
    return df


def visualize(**kwargs):
    name = kwargs.get('dataset')
    global db_connection
    db_connector(**kwargs)
    client = db_connection.client
    if name is None:
        df = client.query(f"""
                SELECT dataset, dataset_name, table, table_name, num_rows, num_cols, size, treemap
                FROM `{db_info}`
                WHERE treemap.show = true AND num_rows > 0""").to_dataframe()

        visualize_all(df)
        return
    df = client.query(f"SELECT dataset_name, table_name, num_rows, num_cols, size FROM `{db_info}`").to_dataframe()
    tables = df.where(df['dataset_name'] == name.upper()).dropna()
    if len(tables) <= 0:
        return
    visualize_dataset(tables, name)


def db_connector(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)
    return db_connection
