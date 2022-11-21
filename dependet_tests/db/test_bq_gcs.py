import unittest
import dotenv
import lib310
import time
import dask.dataframe as dd
from google.cloud import storage
import gcsfs
import google.auth
import os
import pandas as pd


from google.cloud import bigquery


class TestBqGcs(unittest.TestCase):

    def setUp(self):
        dotenv.load_dotenv('.env')
        self.client = bigquery.Client()


    def test_export_to_gcs(self):
        bucket_name = 'bigquery_1'
        project = "pfsdb3"
        dataset_id = "system"
        table_id = "info_test"

        destination_uri = "gs://{}/{}".format(bucket_name, "onefile.csv")
        dataset_ref = bigquery.DatasetReference(project, dataset_id)
        table_ref = dataset_ref.table(table_id)
        extract_job = self.client.extract_table(
            table_ref,
            destination_uri,
            # Location must match that of the source table.
            # location="US",
        )  # API request
        extract_job.result()  # Waits for job to complete.

        print(
            "Exported {}:{}.{} to {}".format(project, dataset_id, table_id, destination_uri)
        )

    def test_simple_query(self):
        res = lib310.db.fetch('SELECT * FROM `pfsdb3.system.info`')
        print(res)
    def test_create_table(self):
        lib310.db.cache_query('SELECT * EXCEPT(treemap) FROM `pfsdb3.system.info`', 'salam')

    def test_read_from_gcs(self):
        # # Instantiates a client
        storage_client = storage.Client()
        print(storage_client.SCOPE)
        #
        # # The name for the new bucket
        # bucket_name = "bigquery_1"
        #
        # x = storage_client.get_bucket(bucket_name).get_blob('lang_4/tokens/tokens_000000000000')
        # print(x.download_as_string())
        # for z in x:
        #     print(z)

        # z = pd.read_csv('gs://bigquery_1/lang_4/tokens/tokens_000000000000')
        # print(z)
        # x = dd.read_csv('gs://bigquery_1/lang_4/tokens/tokens_000000000000', storage_options={'token': '../google-service.json'})
        # print(x)
        print(os.getenv('GOOGLE_APPLICATION_CREDENTIALS'))
        # credentials, x = google.auth.load_credentials_from_file(os.getenv('GOOGLE_APPLICATION_CREDENTIALS'), scopes=storage_client.SCOPE)
        credentials, x = google.auth.default(scopes=storage_client.SCOPE)
        # print(x)
        # print(credentials)
        print(credentials.scopes)
        fs = gcsfs.GCSFileSystem(project='pfsdb3', token=credentials)
        print(fs.scopes)
        print(fs.credentials)
        fs.credentials.maybe_refresh()
        print(fs.credentials.token)
        print(fs.credentials.token.valid)
        # print(fs.credentials.project)
        # print(fs.credentials.access)
        # print(fs.credentials.scope)
        # print(fs.credentials.method)
        # print(fs.protocol)
        # print(fs.buckets)
        # exit(1)
        #
        ddf = dd.read_csv('gs://bigquery_1/lang_4/tokens/tokens_000000000000', storage_options={'token': fs.credentials.token})
        print(ddf)
        ddf = ddf.repartition(partition_size="1MB", force=True)
        print(ddf.index.pipe(len))
        # print('sssss', len(ddf.partitions[0]))
        # print('aaaaa', ddf.map_partitions(len).compute())
        parts = ddf.map_partitions(len).compute()
        idx = 4000
        n = 0
        division_start = 0
        for part in parts:
            if idx < division_start + part:
                break
            n += 1
            division_start += part
        df = ddf.get_partition(n).compute()
        print(division_start, idx)
        print(df.iloc[idx - division_start, df.columns.get_indexer(df.columns[df.columns != 'num'])])


        # print(ddf.columns.get_indexer(ddf.columns[ddf.columns != 'num']))
        # item = ddf.iloc[5, 1]
        # print(item)
