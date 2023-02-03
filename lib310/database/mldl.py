from pymongo import MongoClient, ASCENDING
import random
import time
import pandas as pd
from threading import Thread, Event
from queue import Queue


class MLDL:
    def __init__(self, mongo_url, cache_size=10):
        random.seed(time.time())
        self.client = MongoClient(mongo_url)
        self.db = self.client['A310']
        self.col_train = self.db['seq_lang5_train']
        self.col_test = self.db['seq_lang5_test']
        self.cache = {'TRAIN': {}, 'TEST': {}}
        self.sample_cache = []
        self.max_queue_size = cache_size
        self.queue_key = None
        self.queue = Queue(self.max_queue_size)
        self.filler_thread = None
        self.kill_event = Event()
        self.kill_event_cache = Event()
        self.filler_thread_cache = None

        tmp = self.db['length_freq_lang5_train'].find().sort('len', ASCENDING)
        for item in tmp:
            self.cache['TRAIN'][item['len']] = {'freq': item['freq'], 'maxi': item['maxi'], 'mini': item['mini']}
        tmp = self.db['length_freq_lang5_test'].find().sort('len', ASCENDING)
        for item in tmp:
            self.cache['TEST'][item['len']] = {'freq': item['freq'], 'maxi': item['maxi'], 'mini': item['mini']}

    def get_batch(self, num, max_length, min_num_feature=0, stage='TRAIN'):
        stage = stage.upper()
        if stage == 'TRAIN':
            col = self.col_train
        elif stage == 'TEST':
            col = self.col_test
        else:
            raise ValueError('Stage must be either TRAIN or TEST')
        if max_length not in self.cache[stage.upper()].keys():
            arr = [i for i in self.cache[stage.upper()].keys() if i <= max_length]
            if (len(arr) == 0):
                return pd.DataFrame()
            max_length = max(arr)

        hit = self.queue_key == f'{num}_{max_length}_{min_num_feature}_{stage}'
        if hit:
            # HIT
            return pd.DataFrame(self.queue.get())

        # MISS
        maxi = self.cache[stage.upper()][max_length]['maxi']

        if self.filler_thread_cache is not None:
            self.kill_event_cache.set()
            self.filler_thread_cache.join()
            self.kill_event_cache.clear()
            self.filler_thread_cache = None

        if min_num_feature > 11:
            c = col.find({'len': {'$lt': max_length + 1}, 'token_ids': {'$size': {'$gt': min_num_feature - 1}}},
                         {'row_id': 1, '_id': 0}).sort("rand").limit(5 * num)
            self.sample_cache = list(map(lambda x: x['row_id'], list(c)))
            self.filler_thread_cache = Thread(target=self.background_sample_cache,
                                              args=(max_length, num, min_num_feature, col))
            self.filler_thread_cache.start()
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

        index_list = random.sample(self.sample_cache, min(num, len(self.sample_cache)))
        cursor = col.find({'row_id': {'$in': index_list}}, {'_id': 1, 'sequence': 1, 'token_ids': 1})
        return pd.DataFrame(list(cursor))

    def filler(self, num, col):
        while True:
            index_list = random.sample(self.sample_cache, min(num, len(self.sample_cache)))
            cursor = col.find({'row_id': {'$in': index_list}}, {'_id': 1, 'sequence': 1, 'token_ids': 1})
            df = pd.DataFrame(list(cursor))
            self.queue.put(df)
            if self.kill_event.is_set():
                break

    def background_sample_cache(self, max_length, num, min_num_feature, col):
        i = 5
        while i > 0:
            c = col.find({'len': {'$lt': max_length + 1}, 'token_ids': {'$size': {'$gt': min_num_feature - 1}}},
                         {'row_id': 1, '_id': 0}).sort("rand").skip(i * num).limit(num)
            tmp = list(c)
            self.sample_cache += list(map(lambda x: x['row_id'], tmp))
            i += 1
            if len(tmp) < num or self.kill_event_cache.is_set():
                i = 0
                break

    def terminate(self):
        self.kill_event.set()
        self.queue.get()
        self.filler_thread.join()
        self.kill_event.clear()
        self.kill_event_cache.set()
        self.filler_thread_cache.join()
        self.kill_event_cache.clear()
        self.filler_thread_cache = None