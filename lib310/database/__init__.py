from typing import Optional
from ._connection import DatabaseConnection, set_gcloud_key_path
from ._visualize import *
import io
import json

db_connection: DatabaseConnection = None


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
    import pandas as pd
    all_datasets = list_tables()
    result = []
    for dataset_id, tables in all_datasets.items():
        for table in tables:
            row = {}
            temp = db_connection.client.get_table(table.full_table_id.replace(':', '.'))
            row['dataset'] = dataset_id.split('_')[-1].upper()
            row['table'] = temp.table_id
            row['num_rows'] = temp.num_rows
            row['num_cols'] = len(temp.schema)
            size = temp.num_bytes
            for x in ['B', 'KB', 'MB', 'GB', 'TB']:
                if size < 1024.0:
                    row['size'] = "%3.1f %s" % (size, x)
                    break
                size /= 1024.0
                row['size'] = size
            result.append(row)

    print_it = kwargs.get('print')
    df = pd.DataFrame(result)
    if print_it is None or print_it is True:
        print(df.to_string(index = False))
    return df


def visualize(name=None):
    df = summary(print=False)
    if name is None:
        ds = df.where(df['num_rows'] > 0)\
                .where(df['dataset'] != 'MID')\
                .where(df['dataset'] != 'SANDBOX').dropna()
        visualize_all(ds)
        return
    if name.lower() == 'datasets':
        ds = df.where(df['num_rows'] > 0).groupby('dataset').sum().dropna()
        visualize_all_datasets(ds)
        return
    tables = df.where(df['dataset'] == name.upper()).dropna()
    if len(tables) <= 0:
        return
    visualize_dataset(tables, name)


def gather_info(**kwargs):
    global db_connection
    if db_connection is None:
        db_connection = DatabaseConnection(**kwargs)
    client = db_connection.client

    with open('dataset_info.json', 'r') as f:
        dict_info = json.load(f)
    datasets = list(client.list_datasets())
    result = []
    for dataset in datasets:
        if dataset.dataset_id == 'sandbox':
            continue

        tables = list(db_connection.client.list_tables(dataset.dataset_id))
        for t in tables:
            table = client.get_table(t.full_table_id.replace(':', '.'))
            schema = io.StringIO("")
            client.schema_to_json(table.schema, schema)
            item = {
                'dataset': dataset.dataset_id,
                'dataset_name': extract_value(f'{dataset.dataset_id}.name', dict_info, default=dataset.dataset_id),
                'table': t.table_id,
                'full_table_id': t.full_table_id,
                'table_name': extract_value(f'{dataset.dataset_id}.tables.{t.table_id}.name', dict_info, default=t.table_id),
                'num_rows': table.num_rows,
                'num_cols': len(table.schema),
                'size': table.num_bytes,
                'schema': json.dumps(schema.getvalue()),
                'treemap': {
                    'show': extract_value(f'{dataset.dataset_id}.treemap.show', dict_info, default=True),
                    'color': extract_value(f'{dataset.dataset_id}.treemap.children.color', dict_info, default='#FF0000'),
                    'fontsize': extract_value(f'{dataset.dataset_id}.treemap.children.size', dict_info, default=18),
                    'textcolor': extract_value(f'{dataset.dataset_id}.treemap.children.size', dict_info, default='#FFFFFF'),
                },
                'last_modified': table.modified.timestamp() * 1000,
            }
            item['treemap']['show'] = extract_value(f'{dataset.dataset_id}.treemap.tables.{t.table_id}.show', dict_info, default=item['treemap']['show'])
            item['treemap']['color'] = extract_value(f'{dataset.dataset_id}.treemap.tables.{t.table_id}.color', dict_info, default=item['treemap']['color'])
            item['treemap']['fontsize'] = extract_value(f'{dataset.dataset_id}.treemap.tables.{t.table_id}.fontsize', dict_info, default=item['treemap']['fontsize'])
            item['treemap']['textcolor'] = extract_value(f'{dataset.dataset_id}.treemap.tables.{t.table_id}.textcolor', dict_info, default=item['treemap']['textcolor'])

            result.append(item)
    result = [json.dumps(record) for record in result]

    with open('info.json', 'w') as obj:
        for i in result:
            obj.write(i + '\n')



def extract_value(key, dictionary, default=None):
    parts = key.split('.')
    value = dictionary
    for i in range(len(parts)):
        if parts[i] not in value:
            return default
        value = value[parts[i]]
    return value
