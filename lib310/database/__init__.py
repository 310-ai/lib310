import time
from typing import Optional
from ._connection import DatabaseConnection, set_gcloud_key_path
from ._visualize import *
from ._constants import FileFormat
from datetime import datetime, timedelta
from .gcs_dataset import GCSDataset
import logging


db_connection: DatabaseConnection = None
__db_info = 'system.info'
__db_gcs_cache = 'system.gcs_cache'

log = logging.getLogger("lib310")


def fetch(*args, **kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection()

    query_results = db_connection.query(*args, **kwargs)
    return query_results


def get_table_info(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection()
    elif db_connection.table_name != kwargs.get('table_name', db_connection.table_name):
        db_connection = DatabaseConnection()

    return db_connection.get_table_info()


def list_datasets(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection()
    return list(db_connection.client.list_datasets())


def list_tables(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection()

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
    db_connector()
    client = db_connection.client
    df = client.query(f"SELECT dataset_name, table_name, num_rows, num_cols, size FROM `{__db_info}`").to_dataframe()

    print_it = kwargs.get('print')
    if print_it is None or print_it is True:
        df['size'] = df['size'].apply(size_to_abbrevation)
        print(df.to_string(index = False))
    return df


def visualize(**kwargs):
    name = kwargs.get('dataset')
    global db_connection
    db_connector()
    client = db_connection.client
    if name is None:
        df = client.query(f"""
                SELECT dataset, dataset_name, table, table_name, num_rows, num_cols, size, treemap
                FROM `{__db_info}`
                WHERE treemap.show = true AND num_rows > 0""").to_dataframe()

        visualize_all(df)
        return
    df = client.query(f"SELECT dataset_name, table_name, num_rows, num_cols, size FROM `{__db_info}`").to_dataframe()
    tables = df.where(df['dataset_name'] == name.upper()).dropna()
    if len(tables) <= 0:
        return
    visualize_dataset(tables, name)


def db_connector():
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection()
    return db_connection


def cache_query(query, name, destination_format=FileFormat.CSV):
    db = db_connector()
    table = db.client.build_temp_table()
    log.debug(f'Created table {table.table_id} in {db.client.CACHE_DATASET}')

    db.client.query_to_cached_dataset(query=query, destination=table)
    log.debug(f'Insert query {query} to {table.table_id}')

    destination_format = destination_format.upper()
    try:
        destination_format = FileFormat.to_format(destination_format)
    except KeyError:
        raise TypeError('The destination format is not supported')

    res = db.client.export_to_gcs(table, name, destination_format)
    log.debug(f'Export {table.table_id} to GCS(/{table.table_id}/{name}_*.{destination_format})')

    try:
        info = db.client.insert_rows_json(__db_gcs_cache, [{
            'name': name,
            'folder': table.table_id,
            'uri': res.destination_uris[0],
            'length': res.destination_uri_file_counts[0],
            'created_at': datetime.now().isoformat(),
            'expired_at': (datetime.now() + timedelta(days=7)).isoformat()
        }])
        if len(info) > 0 and len(info[0]['errors']) != 0:
            log.error(info[0]['errors'])
    except Exception as e:
        log.error(e)


    db.client.delete_table(table)
    log.debug(f'Deleted table {table.table_id} from {db.client.CACHE_DATASET}')

    return res.destination_uris, res.destination_uri_file_counts
