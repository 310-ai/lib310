import unittest
import lib310
import dotenv


class TestDB(unittest.TestCase):

    def setUp(self):
        dotenv.load_dotenv('.env')

    def test_bigquery_select(self):
        seqs1 = lib310.db.fetch(query="SELECT GO_ID FROM `pfsdb3.0_go.gaf` LIMIT 10")
        seqs2 = lib310.db.fetch("SELECT GO_ID FROM `pfsdb3.0_go.gaf` LIMIT 10")
        self.assertEqual(len(seqs1), 10)
        self.assertEqual(len(seqs2), 10)


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

    def test_summary(self):
        tables = lib310.db.summary()

    def test_db_visualize(self):
        # tables = lib310.db.summary()
        # lib310.db.visualize()
        lib310.db.visualize(dataset='go')
        lib310.db.visualize(dataset='interpro')
