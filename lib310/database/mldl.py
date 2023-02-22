from pymongo import MongoClient, ASCENDING
import random
import time
import pandas as pd
from threading import Thread, Event
from queue import Queue
import concurrent.futures
# from . import cache_query, FileFormat
from .cache_bgquery import cache_query, FileFormat
import dask.dataframe as dd
from dotenv import load_dotenv

load_dotenv()

TRAIN = 'TRAIN'
TEST = 'TEST'
MAIN = 'MAIN'
INTERACTION = 'INTERACTION'


class MLDL:
    def __init__(self, mongo_url, cache_size=10):
        random.seed(time.time())
        self.client = MongoClient(mongo_url)
        self.db = self.client['A310']
        self.collections = {
            TRAIN: {
                MAIN: self.db['seq_lang5_train'],
                INTERACTION: self.db['interactions_lang5_train'],
            },
            TEST: {
                MAIN: self.db['seq_lang5_test'],
                INTERACTION: self.db['interactions_lang5_test'],
            }
        }

        self.cache = {TRAIN: {}, TEST: {}}
        self.sample_cache = []
        self.max_queue_size = cache_size
        self.queue_key = None
        self.queue = Queue(self.max_queue_size)
        self.filler_thread = None
        self.kill_event = Event()
        # self.kill_event_cache = Event()
        # self.filler_thread_cache = None
        # self.sample_queue = Queue()
        # self.sample_max_worker = sample_max_worker

        # Cache the data with regard to length of sequence
        tmp = self.db['length_freq_lang5_train'].find().sort('len', ASCENDING)
        for item in tmp:
            self.cache[TRAIN][item['len']] = {'freq': item['freq'], 'maxi': item['maxi'], 'mini': item['mini']}
        tmp = self.db['length_freq_lang5_test'].find().sort('len', ASCENDING)
        for item in tmp:
            self.cache[TEST][item['len']] = {'freq': item['freq'], 'maxi': item['maxi'], 'mini': item['mini']}

    def get_batch(self, num, max_length, min_num_feature=0, stage=TRAIN):
        """
        Get a batch of data from the database with regard to the length of the sequence and the number of features
        It also start the background thread to fill the caches and queue
        :param num: number of samples
        :param max_length: maximum length of the sequence
        :param min_num_feature: minimum number of features
        :param stage: TRAIN or TEST
        :return: a pandas dataframe
        """
        stage = stage.upper()
        if stage == TRAIN:
            col = self.collections[TRAIN]
        elif stage == TEST:
            col = self.collections[TEST]
        else:
            raise ValueError('Stage must be either TRAIN or TEST')

        if max_length not in self.cache[stage.upper()].keys():
            arr = [i for i in self.cache[stage.upper()].keys() if i <= max_length]
            if len(arr) == 0:
                return pd.DataFrame()
            max_length = max(arr)

        hit = self.queue_key == f'{num}_{max_length}_{min_num_feature}_{stage}'
        if hit:
            # HIT
            return pd.DataFrame(self.queue.get())

        # MISS
        # if self.filler_thread_cache is not None:
        #     self.kill_event_cache.set()
        #     self.filler_thread_cache.join()
        #     self.kill_event_cache.clear()
        #     self.filler_thread_cache = None
        #     self.sample_cache = []

        maxi = self.cache[stage.upper()][max_length]['maxi']

        if min_num_feature > 11:
            table = '1_uniprot.mongo_lang5_train'
            if stage == TEST:
                table = '1_uniprot.mongo_lang5_test'

            res = cache_query(
                query=f"SELECT row_id FROM `{table}` WHERE len <= {max_length} and token_size >= {min_num_feature}",
                name=f'{max_length}_{min_num_feature}_{stage}',
                destination_format=FileFormat.CSV,
            )
            # print(res)
            ddf = dd.read_csv(res['uri'])
            df = ddf.compute()
            self.sample_cache = df['row_id'].tolist()
            # print(len(self.sample_cache))

            # factor = 10
            # if num < 1000:
            #     factor = 10000 // num
            # start = time.perf_counter()
            # c = col[MAIN].find({'len': {'$lt': max_length + 1}, 'token_size': {'$gt': min_num_feature - 1}},
            #                    {'row_id': 1, '_id': 0, 'rand': 1}, batch_size=1000).sort("rand", -1).limit(factor * num)
            # print(f'Query time: {time.perf_counter() - start}')
            # start = time.perf_counter()
            # tmp = list(c)
            # print(f'Convert time: {time.perf_counter() - start}')
            # print(tmp[0], tmp[-1])
            # # print(len(tmp))
            # start = time.perf_counter()
            # self.sample_cache = list(map(lambda x: x['row_id'], tmp))
            # print(f'Convert time 2: {time.perf_counter() - start}')

            # self.filler_thread_cache = Thread(target=self.background_sample_cache,
            #                                   args=(max_length, num, min_num_feature, col[MAIN], tmp[-1]['rand']))
            # self.filler_thread_cache.start()
        else:
            self.sample_cache = range(1, maxi)

        if self.filler_thread is not None:
            self.kill_event.set()
            self.queue.get()
            self.filler_thread.join()
            self.kill_event.clear()
            self.queue = Queue(self.max_queue_size)

        self.queue_key = f'{num}_{max_length}_{min_num_feature}_{stage}'

        self.filler_thread = Thread(target=self.filler, args=(num, col))
        self.filler_thread.start()

        return self.fetch_sample(num, col)

    def filler(self, num, col):
        """
        This function is used to fill the queue with data then get batch read data from the queue
        :param num: number of samples
        :param col: collections of the data
        """
        while True:
            df = self.fetch_sample(num, col)
            self.queue.put(df)
            if self.kill_event.is_set():
                break

    def fetch_sample(self, num, collections):
        """
        This function is used to fetch the data from the database
        Run two threads to fetch the data from the database
        :param num: number of samples
        :param collections: collections of the data
        :return: dataframe of the data
        """
        index_list = random.sample(self.sample_cache, min(num, len(self.sample_cache)))
        with concurrent.futures.ThreadPoolExecutor() as executor:
            main_thread = executor.submit(self.fetch_sample_main, index_list, collections[MAIN])
            interaction_thread = executor.submit(self.fetch_sample_interaction, index_list, collections[INTERACTION])

            main_df = main_thread.result()
            interaction_df = interaction_thread.result()

            return pd.merge(main_df, interaction_df, on='row_id', how='left')

    def fetch_sample_main(self, index_list, col):
        """
        This function is used to fetch the main data from the database
        :param index_list: list of row_ids
        :param col: main collection
        :return: dataframe of main data
        """
        cursor = col.find({'row_id': {'$in': index_list}}, {'_id': 0, 'row_id': 1, 'sequence': 1, 'token_ids': 1})
        return pd.DataFrame(list(cursor))

    def fetch_sample_interaction(self, index_list, col):
        """
        This function is used to fetch the interaction data from the database
        :param index_list: list of row_ids
        :param col: interaction collection
        :return: dataframe of interaction data
        """
        cursor = col.find({'row_id': {'$in': index_list}}, {'_id': 0, 'row_id': 1, 'interactions': 1})
        return pd.DataFrame(list(cursor))

    # def background_sample_cache(self, max_length, num, min_num_feature, col, last_rand):
    #     """
    #     This function is used to increase the sample cache (suitable row_ids) in the background
    #     :param max_length: max length of sequence
    #     :param num: number of samples
    #     :param min_num_feature: min number of features
    #     :param col: main collection
    #     :param last_rand: last rand value
    #     :return: extend the sample cache
    #     """
    #     max_rand = last_rand
    #     num = max(num, 1000)
    #     while True:
    #         print('filling the cache')
    #         c = col.find({'len': {'$lt': max_length + 1}, 'token_size': {'$gt': min_num_feature - 1}, 'rand': {'$lt': max_rand}},
    #                     {'row_id': 1, '_id': 0, 'rand': 1}, batch_size=1000).sort("rand", -1).limit(num)
    #         tmp = list(c)
    #         max_rand = tmp[-1]['rand']
    #         self.sample_cache += list(map(lambda x: x['row_id'], tmp))
    #         if len(tmp) < num or self.kill_event_cache.is_set():
    #             break
    #
    #     # print('starting the background sample cache')
    #     max_workers = self.sample_max_worker
    #     kill_event = Event()
    #     i = 10
    #     if num < 1000:
    #         i = 10000 // num
    #     listener = concurrent.futures.ThreadPoolExecutor()
    #     # print('starting the listener')
    #     listener_thread = listener.submit(self.background_sample_cache_listener, kill_event, max_workers)

    #     # print('starting the executor')

    #     executor = BoundedThreadPoolExecutor(max_workers=max_workers)
    #     while i > 0 and not kill_event.is_set():
    #         # print(f'getting the sample cache {i}')
    #         executor.submit(self.background_sample_cache_filler, max_length, num, min_num_feature, col, i)
    #         i += 1

    #     # print('shutting down the executor')
    #     listener_thread.result()
    #     executor.shutdown(wait=True, cancel_futures=True)

    #     # while i > 0:
    #     #     c = col.find({'len': {'$lt': max_length + 1}, 'token_size': {'$gt': min_num_feature - 1}},
    #     #                  {'row_id': 1, '_id': 0}).sort("rand").skip(i * num).limit(num)
    #     #     tmp = list(c)
    #     #     self.sample_cache += list(map(lambda x: x['row_id'], tmp))
    #     #     i += 1
    #     #     if len(tmp) < num or self.kill_event_cache.is_set():
    #     #         i = 0
    #     #         break

    # def background_sample_cache_filler(self, max_length, num, min_num_feature, col, i, min_rand, max_rand):
    #     """
    #     This function is used to increase the sample cache (suitable row_ids) in the background
    #     :param max_length: max length of sequence
    #     :param num: number of samples
    #     :param min_num_feature: min number of features
    #     :param col: main collection
    #     :param i: index of the skip
    #     :param last_rand: last rand value
    #     :return: extend the sample cache
    #     """
    #     min_rand = min_rand
    #     c = col.find({'len': {'$lt': max_length + 1}, 'token_size': {'$gt': min_num_feature - 1}, 'rand': {'$gt': min_rand}'},
    #                  {'row_id': 1, '_id': 0}).sort("rand").skip(i * num).limit(num)
    #     tmp = list(c)
    #     tmp = list(map(lambda x: x['row_id'], tmp))
    #     if self.kill_event_cache.is_set():
    #         self.sample_queue.put([])
    #         return

    #     self.sample_queue.put(tmp)
    #     if len(tmp) < num:
    #         self.sample_queue.put([])

    # def background_sample_cache_listener(self, kill_event, max_workers):
    #     i = 0
    #     while True:
    #         try:
    #             tmp = self.sample_queue.get(timeout=3 * 60)
    #             # print(f'get sample cache with len {len(tmp)}')
    #             if self.kill_event_cache.is_set():
    #                 self.sample_queue = Queue()
    #                 kill_event.set()
    #                 return
    #             if len(tmp) == 0:
    #                 # print('sample cache is empty')
    #                 kill_event.set()
    #                 i += 1

    #             self.sample_cache += tmp

    #             if i >= max_workers:
    #                 # print('workers are all done')
    #                 return

    #         except Exception as error:
    #             # print('time out on get sample cache')
    #             kill_event.set()
    #             self.sample_queue = Queue()
    #             return

    def terminate(self):
        """
        Terminate the filler threads and clear the queue
        """
        # Terminate the filler thread
        if self.filler_thread is not None:
            self.kill_event.set()
            self.queue.get()
            self.filler_thread.join()
            self.kill_event.clear()
            self.queue = Queue(self.max_queue_size)

        # Terminate the filler thread cache
        # if self.filler_thread_cache is not None:
        #     self.kill_event_cache.set()
        #     self.filler_thread_cache.join()
        #     self.kill_event_cache.clear()
        #     self.filler_thread_cache = None
        #     self.sample_cache = []
