import unittest
import lib310
import dotenv


class TestDB(unittest.TestCase):

    def setUp(self):
        dotenv.load_dotenv('.env')

    def test_bigquery_select(self):
        seqs = lib310.db.fetch(query="SELECT * FROM `pfsdb3.0_go.gaf` LIMIT 10")
        self.assertEqual(len(seqs), 10)

    def test_get_datasets(self):
        datasets = lib310.db.get_datasets()
