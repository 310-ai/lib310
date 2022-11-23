import torch
from google.cloud import storage
import google.auth
import gcsfs
from torch.utils.data import Dataset
import dask.dataframe as dd
from ._constants import FileFormat


class GCSDataset(Dataset):
    def __init__(self, file_uri, target_col_name, file_format=FileFormat.CSV):
        """
        :param file_uri: the uri of file even with wildcards in gcs
        :return:
        """
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
