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
        self.max_queue_size = cache_size
        self.queue_key = None
        self.queue = Queue(self.max_queue_size)
        self.filler_thread = None
        self.kill_event = Event()

        tmp = self.db['length_freq_lang5_train'].find().sort('len', ASCENDING)
        for item in tmp:
            self.cache['TRAIN'][item['len']] = {'freq': item['freq'], 'maxi': item['maxi'], 'mini': item['mini']}
        tmp = self.db['length_freq_lang5_test'].find().sort('len', ASCENDING)
        for item in tmp:
            self.cache['TEST'][item['len']] = {'freq': item['freq'], 'maxi': item['maxi'], 'mini': item['mini']}
        
    
    def get_batch(self, num, max_length, stage):
        stage = stage.upper()
        if stage == 'TRAIN':
            col = self.col_train
        elif stage == 'TEST':
            col = self.col_test
        else:
            raise ValueError('Stage must be either TRAIN or TEST')
        if max_length not in self.cache[stage.upper()].keys():
            arr = [i for i in self.cache[stage.upper()].keys() if i <= max_length]
            if(len(arr) == 0):
                return pd.DataFrame()
            max_length = max(arr)
        
        hit = self.queue_key == f'{num}_{max_length}_{stage}'
        if hit:
            return pd.DataFrame(self.queue.get())

        if self.filler_thread is not None:
            self.kill_event.set()
            self.queue.get()
            self.filler_thread.join()
            self.kill_event.clear()
            self.queue = Queue(self.max_queue_size)

        self.queue_key = f'{num}_{max_length}_{stage}'

        maxi = self.cache[stage.upper()][max_length]['maxi']
        self.filler_thread = Thread(target=self.filler, args=(self.queue, num, max_length, maxi, col, self.kill_event))
        self.filler_thread.start()

        index_list = random.sample(range(1, maxi), num)
        cursor = col.find({'len': {'$lt': max_length + 1}, 'row_id': {'$in': index_list}}).limit(num)
        return pd.DataFrame(list(cursor))

    @staticmethod
    def filler(queue, num, max_length, maxi, col, event):
        while True:
            index_list = random.sample(range(1, maxi), num)
            cursor = col.find({'len': {'$lt': max_length + 1}, 'row_id': {'$in': index_list}}).limit(num)
            df = pd.DataFrame(list(cursor))
            queue.put(df)
            if event.is_set():
                break
    def terminate(self):
        self.kill_event.set()
        self.queue.get()
        self.filler_thread.join()
        self.kill_event.clear()
