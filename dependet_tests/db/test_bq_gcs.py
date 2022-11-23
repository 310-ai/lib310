import unittest
import dotenv
import lib310


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

    def test_cache_query_csv(self):
        dst, n = lib310.db.cache_query('SELECT * EXCEPT(treemap) FROM `pfsdb3.system.info`', 'mycsv')
        print(dst)
        print(n)

    def test_cache_query_json(self):
        dst, n = lib310.db.cache_query('SELECT * FROM `pfsdb3.system.info`', 'myjson', 'json')
        print(dst)
        print(n)

    def test_read_from_gcs(self):
        ddf = lib310.db.GCSDataset('gs://bigquery_1/lang_4/tokens/tokens_000000000000', 'num')
        print(len(ddf))
        print(ddf[4000])

    def test_read_from_gcs_json(self):
        ddf = lib310.db.GCSDataset('gs://bigquery_1/QJZDCWPD/myjson_*.json', 'num', 'json')
        print(len(ddf))

