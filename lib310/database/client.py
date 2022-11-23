from google.cloud import bigquery
from google.cloud.bigquery import retry as r
from google.cloud.bigquery import enums
from google.api_core.exceptions import AlreadyExists, Conflict

from . import exceptions
from ._constants import FileFormat
from datetime import datetime
from . import _functions as fn


class Client(bigquery.Client):
    PROJECT = 'pfsdb3'
    CACHE_DATASET = 'cached'
    BUCKET_NAME = 'bigquery_1'

    def __init__(self):
        super(Client, self).__init__()

    def query(self,
              query: str,
              job_config=None,
              job_id: str = None,
              job_id_prefix: str = None,
              location: str = None,
              project: str = None,
              retry=r.DEFAULT_RETRY,
              timeout=r.DEFAULT_TIMEOUT,
              job_retry=r.DEFAULT_JOB_RETRY,
              api_method=enums.QueryApiMethod.QUERY):

        dry_config = bigquery.QueryJobConfig(dry_run=True, use_query_cache=False)
        dry_run = super(Client, self).query(query, dry_config, job_id, job_id_prefix, location, project, retry, timeout, job_retry, api_method)
        vol = int(dry_run.total_bytes_processed)
        # if vol > self.single_query_limit:
        #     raise exceptions.VolumeLimitException()
        self.__write_usage_log(vol)
        return super(Client, self).query(query, job_config, job_id, job_id_prefix, location, project, retry, timeout, job_retry, api_method)

    def build_temp_table(self):
        """
        making a table in the cached dataset
        :return: table
        """
        random_string = fn.get_random_string(8)
        ok = False
        table = None
        while not ok:
            try:
                table = self.create_table(f'{self.PROJECT}.{self.CACHE_DATASET}.{random_string}')
                ok = True
            except Conflict:
                random_string += fn.get_random_string(1)
        return table

    def query_to_cached_dataset(self, query: str, destination: bigquery.Table):

        config = bigquery.QueryJobConfig(allow_large_results=True)
        config.destination = destination
        config.write_disposition = bigquery.WriteDisposition().WRITE_TRUNCATE

        query_result = self.query(query, job_config=config, api_method=enums.QueryApiMethod.INSERT)
        return query_result.result()

    def export_to_gcs(self, table, name, destination_format=FileFormat.CSV):
        try:
            destination_uri = f"gs://{self.BUCKET_NAME}/{table.table_id}/{name}_*.{FileFormat.to_extension(destination_format)}"
        except KeyError:
            destination_uri = f"gs://{self.BUCKET_NAME}/{table.table_id}/{name}_*"
        dataset_ref = bigquery.DatasetReference(self.PROJECT, self.CACHE_DATASET)
        table_ref = dataset_ref.table(table.table_id)
        try:
            job_config = bigquery.job.ExtractJobConfig()
            job_config.destination_format = FileFormat.to_format(destination_format)
        except KeyError:
            raise TypeError('The destination format is not supported')
        extract_job = self.extract_table(
            table_ref,
            destination_uri,
            job_config=job_config
            # Location must match that of the source table.
            # location="US",
        )  # API request
        res = extract_job.result()  # Waits for job to complete.
        return res

    def __write_usage_log(self, vol):
        try:
            with open(self.__usage_log_file, 'a', encoding='utf-8') as f:
                f.write(f"{datetime.now()}: {vol} bytes\n")
                return True
        except Exception as e:
            return False

    @property
    def __usage_log_file(self):
        return '.usage_log'

    @property
    def single_query_limit(self):
        return 3 * 10 ** 12