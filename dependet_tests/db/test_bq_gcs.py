import unittest
import dotenv
import lib310
from torch.utils.data import DataLoader
from lib310.database import FileFormat, CacheResponseType
from torch.utils.data import DataLoader


from google.cloud import bigquery


def salam(data):
    print(data)
    return data
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
        row = lib310.db.cache_query('SELECT * EXCEPT(treemap) FROM `pfsdb3.system.info`', 'mycsv')
        print(row)

    def test_cache_query_json(self):
        row = lib310.db.cache_query('SELECT * FROM `pfsdb3.system.info`', 'myjson', 'json')
        print(row)

    def test_read_from_gcs(self):
        ddf = lib310.db.GCSDataset('gs://bigquery_1/lang_4/tokens/tokens_000000000000', 'num')
        print(len(ddf))
        print(ddf[4000])

    def test_read_from_gcs_json(self):
        ddf = lib310.db.GCSDataset('gs://bigquery_1/QJZDCWPD/myjson_*.json', 'num', 'json')
        print(len(ddf))

    def test_direct_dataset(self):
        ddf = lib310.db.cache_query('SELECT * FROM `pfsdb3.system.info`', 'myjson', 'json', response_type=lib310.db.CacheResponseType.DATASET)
        print(len(ddf))

    def test_simple_dataloader(self):
        row = lib310.db.cache_query('SELECT * EXCEPT(treemap) FROM `pfsdb3.system.info`', 'mycsv')
        print(row['uri'])
        ddf = lib310.db.GCSDataset(row['uri'], 'size')
        data = DataLoader(ddf, batch_size=2, shuffle=False)
        for batch in data:
            print(batch)

    def test_icm_creating(self):
        ds = lib310.db.cache_query(
            query='''
               SELECT parc.Entry, terms.go_ids
               FROM `1_go.sequence_term_expanded` terms
               JOIN `0_uniprot.uniparc` parc
               ON terms.sequence = parc.Sequence
           ''',
            name='ICR',
            destination_format=FileFormat.JSON,
            response_type=CacheResponseType.DATASET,
        )
        print(ds)
        data = DataLoader(ds, batch_size=16, shuffle=False)
        for batch in data:
            print(batch)


