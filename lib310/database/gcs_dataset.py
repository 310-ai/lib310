import torch
from google.cloud import storage
import google.auth
import gcsfs
from torch.utils.data import Dataset
import dask.dataframe as dd


class GCSDataset(Dataset):
    def __int__(self, file_uri, target_col_name):
        """
        :param file_uri: the uri of file even with wildcards in gcs
        :return:
        """
        storage_client = storage.Client()
        credentials, project = google.auth.default(scopes=storage_client.SCOPE)
        fs = gcsfs.GCSFileSystem(project=project, token=credentials)
        fs.credentials.maybe_refresh()
        self.ddf = dd.read_csv('gs://bigquery_1/lang_4/tokens/tokens_000000000000', storage_options={'token': fs.credentials.token})
        self.target_name = target_col_name

    def __len__(self):
        return self.ddf.index.pipe(len)

    def __getitem__(self, idx):
        parts = self.ddf.map_partitions(len).compute()
        n = 0
        division_start = 0
        for part in parts:
            if idx < division_start + part:
                break
            n += 1
            division_start += part
        part_idx = idx - division_start
        df = self.ddf.get_partition(n).compute()
        item = df.iloc[part_idx, df.columns.get_indexer(df.columns[df.columns != self.target_name])]
        target = df.iloc[part_idx, df.columns.get_indexer(df.columns[df.columns == self.target_name])]
        return item, target
