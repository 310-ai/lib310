import unittest
import dotenv
import lib310
import time

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
