from google.cloud import storage
import google.auth
import gcsfs
from torch.utils.data import Dataset
import dask.dataframe as dd
from ._constants import FileFormat
import threading



class GCSDataset(Dataset):
    def __init__(self, file_uri, target_col_name, file_format=FileFormat.CSV, size=-1):
        """
        :param file_uri: the uri of file even with wildcards in gcs
        :return:
        """
        if not size:
            size = -1
        storage_client = storage.Client()
        credentials, project = google.auth.default(scopes=storage_client.SCOPE)
        fs = gcsfs.GCSFileSystem(project=project, token=credentials)
        fs.credentials.maybe_refresh()
        file_format = FileFormat.to_format(file_format)
        if file_format == FileFormat.CSV:
            self.ddf = dd.read_csv(file_uri, storage_options={'token': fs.credentials.token})
        elif file_format == FileFormat.JSON:
            self.ddf = dd.read_json(file_uri, storage_options={'token': fs.credentials.token})
        elif file_format == FileFormat.PARQUET:
            self.ddf = dd.read_parquet(file_uri, storage_options={'token': fs.credentials.token})
        else:
            raise TypeError('The file_format is not supported. supported("csv", "json", "parquet")')
        self.target_name = target_col_name
        self.file_uri = file_uri
        self.lock = threading.Lock()


        if size < 0:
            self.n = 0
            self.parts = []
            for i in range(self.ddf.npartitions):
                df = self.ddf.partitions[i].compute()
                self.n = self.n + len(df)
                self.parts.append({
                    'fetched': True,
                    'start': self.n - len(df),
                    'end': self.n - 1
                })
            self.last_partition_loaded = self.ddf.get_partition(0).compute()
            self.last_partition_loaded_index = 0
        else:
            self.last_partition_loaded = self.ddf.get_partition(0).compute()
            self.last_partition_loaded_index = 0
            first_partition_size = len(self.last_partition_loaded)
            self.n = size
            self.parts = [
                {
                    'fetched': (i == 0),
                    'start': (i * first_partition_size),
                    'end': min(size, ((i + 1) * first_partition_size) - 1)
                } for i in range(self.ddf.npartitions)]

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        try:
            _, part_idx, df = self.find_correct_dask_partition(idx)
            if part_idx >= len(df):
                print('idx: ', idx)
                return None
            item = df.iloc[part_idx, df.columns.get_indexer(df.columns[df.columns != self.target_name])]
            if self.target_name is not None and self.target_name != '' and self.target_name in df.columns:
                target = df.iloc[part_idx, df.columns.get_indexer(df.columns[df.columns == self.target_name])]
                target = target.to_list()
                return dict(item), dict(target)
            return dict(item)
        except Exception as e:
            print(e)
            return None

    def find_correct_dask_partition(self, idx):
        self.lock.acquire()
        for i, part in enumerate(self.parts):
            if not part['fetched']:
                prv = self.parts[i - 1]
                if self.last_partition_loaded_index != i:
                    self.last_partition_loaded = self.ddf.partitions[i].compute()
                    self.last_partition_loaded_index = i
                tmp = {
                    'fetched': True,
                    'start': prv['end'] + 1,
                    'end': prv['end'] + len(self.last_partition_loaded)
                }
                self.parts[i] = tmp

            if part['start'] <= idx <= part['end']:
                if self.last_partition_loaded_index != i:
                    lp = self.ddf.partitions[i].compute()
                    li = i
                    self.last_partition_loaded = lp
                    self.last_partition_loaded_index = li
                else:
                    li = self.last_partition_loaded_index
                    lp = self.last_partition_loaded
                part_idx = idx - part['start']
                result = (li, part_idx, lp)
                self.lock.release()
                return result

        part_idx = idx - self.parts[-1]['start']
        li = self.last_partition_loaded_index
        lp = self.last_partition_loaded
        result = (li, part_idx, lp)
        self.lock.release()
        return result
