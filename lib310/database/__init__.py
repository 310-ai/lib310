import matplotlib.pyplot as plt
import seaborn as sns
from .charts.bubble_chart import BubbleChart
from .charts.constants import Component
from typing import Optional
from ._connection import DatabaseConnection, set_gcloud_key_path
import math

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


def visualize(component):
    if component.lower() != Component.DATASETS.value:
        return
    df = summary()
    ds = df.where(df['num_rows'] > 0).groupby('dataset').sum().dropna()
    bubble_chart = BubbleChart(area=[math.log(x) for x in ds['num_rows']],
                               bubble_spacing=0.1)
    bubble_chart.collapse()

    fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"))
    bubble_chart.plot(ax, list(ds.index), sns.color_palette("Spectral", len(list(ds.index))).as_hex())
    ax.axis("off")
    ax.relim()
    ax.autoscale_view()
    ax.set_title('Dataset (log)')

    plt.show()

