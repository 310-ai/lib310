import unittest
import lib310
import dotenv


class TestDB(unittest.TestCase):

    def setUp(self):
        dotenv.load_dotenv('.env')

    def test_bigquery_select(self):
        seqs = lib310.db.fetch(query="SELECT * FROM `pfsdb3.0_go.gaf` LIMIT 10")
        self.assertEqual(len(seqs), 10)

    def test_list_datasets(self):
        datasets = lib310.db.list_datasets()
        self.assertGreater(len(datasets), 0)
        datasets_ids = [d.dataset_id for d in datasets]
        self.assertIn('0_go', datasets_ids)

    def test_list_tables(self):
        tables = lib310.db.list_tables()
        self.assertGreater(len(tables), 0)
        tables = lib310.db.list_tables(dataset_ids='0_go')
        self.assertGreater(len(tables), 0)
