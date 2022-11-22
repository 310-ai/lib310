import os


def set_gcloud_key_path(path: str):
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = path


class DatabaseConnection(object):
    def __init__(self):
        from .client import Client  # This import is here for lazy load

        if "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ.keys():
            raise Exception(
                'You must set BigQuery credentials first! try using `lib310.db.set_gcloud_key_path(path) method.`')

        self.client = Client()
        # self.table_name = kwargs.get('table_name', 'pfsdb3.0_uniprot.uniref')
        # self.table = self.client.get_table(self.table_name) if self.table_name else None
        # self.table_schema = {schema.name: {'type': schema.field_type, } for schema in
        #                      self.table.schema} if self.table else {}
        # self.verbose = kwargs.pop('verbose', False)
        #
        # if self.verbose:
        #     print(f'Successfully established a connection with lib310 Database!')

    def get_table_info(self):
        for column_name, column_info in self.table_schema.items():
            print(f'{column_name}:{column_info["type"]}')

    def query(self, *args, **query_params):
        from google.cloud import bigquery # This import is here for lazy load
        from ..data._base import ProteinDataTable # This import is here for lazy load
        """
        :param args: can have a string which is query
        :param query_params: if query not mentioned we read other query parameters other wise we only read query
        :param - query
        :param - table_name :default=pfsdb3.0_go.gaf
        :param - limit      :default=100
        :return:
        """

        has_query = False
        if 'query' in query_params.keys():
            has_query = True
            query_results = self.client.query(query_params.get('query')).result().to_dataframe()
        if not has_query and len(args) > 0:
            has_query = True
            query_results = self.client.query(args[0]).result().to_dataframe()

        if not has_query:
            table_name = query_params.pop('table_name', 'pfsdb3.0_go.gaf')
            if table_name != self.table_name:
                self.table_name = table_name
                self.table = self.client.get_table(table_name)
                self.table_schema = {schema.name: {'type': schema.field_type, } for schema in
                                     self.table.schema} if self.table else {}
                if self.verbose:
                    print(f'Successfully connected to {table_name} with {self.table.num_rows} rows!')

            query_string = f"SELECT * FROM `{table_name}`"

            force_limit = query_params.pop('force_limit', True)
            limit = query_params.pop('limit', 100)

            query_parameters = []
            for key, query_value in query_params.items():
                if key in self.table_schema.keys():
                    if not isinstance(query_value, list):
                        query_value = [query_value]
                    query_parameters.append(
                        bigquery.ArrayQueryParameter(f'{key}s', self.table_schema[key]['type'], query_value))

            job_config = bigquery.QueryJobConfig(
                query_parameters=query_parameters,
            )

            constraints = 0
            for key, query_value in query_params.items():
                if key in self.table_schema.keys():
                    if constraints == 0:
                        query_string += ' WHERE '

                    query_string += f'{key} IN UNNEST(@{key}s)'

                    if constraints >= 0:
                        query_string += ' AND '
                    constraints += 1

            query_string = query_string[:-5]

            if force_limit:
                query_string += f' LIMIT {limit}'

            print(query_string)
            query_results = self.client.query(query_string, job_config=job_config).result().to_dataframe()

        query_results = ProteinDataTable.from_pandas(query_results)
        return query_results
