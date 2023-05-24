from ._constants import FileFormat, CacheResponseType
from ._functions import get_random_string
from .gcs_dataset import GCSDataset
import hashlib
from datetime import datetime, timedelta
import logging as log

__db_gcs_cache = 'system.cached_queries'


def cache_query(query,
                name: str = None,
                destination_format: str or FileFormat = FileFormat.CSV,
                days: int = 7,
                ignore_hit: bool = False,
                response_type: CacheResponseType = CacheResponseType.CACHE_INFO,
                target_column='',
                db_connection=None):
    """
    Cache a query result in Google Cloud Storage (GCS)
    :param query: query to run on bigquery and cache in GCS
    :param name: name of the file to store in GCS
    :param destination_format: format of the file to store in GCS
    :param days: days to keep the file in GCS
    :param ignore_hit: ignore cache hit and run the query again
    :param response_type: return response as cache info or as a dataset
    :param target_column: target column to use for the dataset
    :param db_connection: database connection
    :return: cache info or dataset
    """

    if response_type == CacheResponseType.DATASET and FileFormat.to_format(destination_format) == FileFormat.AVRO:
        raise ValueError("Cannot return dataset for AVRO format")
    if str is None:
        name = get_random_string(10)

    hashed_query = hashlib.sha1(query.encode('utf-8')).hexdigest()
    if days <= 0:
        days = 7
    db = db_connection

    # check for the hit
    row = None
    if not ignore_hit:
        cached = db.client.query('SELECT * from `system.cached_queries` WHERE status_code != -1').to_dataframe()
        if len(cached) > 0:
            cached = cached.where(cached['hash'] == hashed_query).dropna().sort_values(by=['created_at'],
                                                                                       ascending=False)
            if len(cached) > 0:
                row = cached.iloc[0].to_dict()
    if row is not None:
        # hit happened
        row['hit'] = True
        if response_type == CacheResponseType.CACHE_INFO:
            return row
        return GCSDataset(row['uri'], target_col_name=target_column,
                          file_format=FileFormat.to_format(destination_format), size=row['total_rows'])

    table = db.client.build_temp_table()
    log.debug(f'Created table {table.table_id} in {db.client.CACHE_DATASET}')

    qresult = db.client.query_to_cached_dataset(query=query, destination=table)
    log.debug(f'Insert query {query} to {table.table_id}')

    destination_format = destination_format.upper()
    try:
        destination_format = FileFormat.to_format(destination_format)
    except KeyError:
        raise TypeError('The destination format is not supported')

    res = db.client.export_to_gcs(table, name, destination_format)
    log.debug(f'Export {table.table_id} to GCS(/{table.table_id}/{name}_*.{destination_format})')

    try:
        row = {
            'name': name,
            'folder': table.table_id,
            'uri': res.destination_uris[0],
            'length': res.destination_uri_file_counts[0],
            'created_at': datetime.now().isoformat(),
            'expired_at': (datetime.now() + timedelta(days=days)).isoformat(),
            'query': query,
            'hash': hashed_query,
            'total_rows': qresult.total_rows,
            'status_code': 1
        }
        info = db.client.insert_rows_json(__db_gcs_cache, [row])
        if len(info) > 0 and len(info[0]['errors']) != 0:
            log.error(info[0]['errors'])
    except Exception as e:
        log.error(e)

    db.client.delete_table(table)
    log.debug(f'Deleted table {table.table_id} from {db.client.CACHE_DATASET}')

    row['hit'] = False

    if response_type == CacheResponseType.CACHE_INFO:
        return row
    return GCSDataset(row['uri'], target_col_name=target_column, file_format=FileFormat.to_format(destination_format), size=row['total_rows'])
