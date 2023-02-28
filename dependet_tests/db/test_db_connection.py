import unittest
import lib310
import dotenv
import os
from lib310.database import MLDL
import time



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

    def test_mldl(test):
        mldl = MLDL(os.getenv('MONGO_URL'))
        x = mldl.get_batch(10, 500, 20, 'TRAIN')
        mldl.terminate()
        print(len(x))

    def test_mldl_multiple_time(test):
        mldl = MLDL(os.getenv('MONGO_URL'))
        for i in range(50):
            start = time.perf_counter()
            x = mldl.get_batch(10, 500, 20, 'TRAIN')
            end = time.perf_counter()
            print(f"Time: {end - start} s")
        print()
        for i in range(50):
            start = time.perf_counter()
            x = mldl.get_batch(100, 500, 20, 'TRAIN')
            end = time.perf_counter()
            print(f"Time: {end - start} s")
        mldl.terminate()

    def test_mldl_multiple_time_sla(test):
        mldl = MLDL(os.getenv('MONGO_URL'))
        d = []
        for i in range(51):
            start = time.perf_counter()
            x = mldl.get_batch(200, 1000, 20, 'TRAIN')
            end = time.perf_counter()
            time.sleep(0.05)
            if i > 0:
                d.append(end - start)
            print(f"Time: {end - start} s")
        mldl.terminate()
        print(sum(d)/len(d))
