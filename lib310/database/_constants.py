from google.cloud.bigquery import DestinationFormat
class CONSTATNS:
    SEQUENCE_KEY = 'Sequence'


class FileFormat(object):
    CSV = DestinationFormat.CSV
    JSON = DestinationFormat.NEWLINE_DELIMITED_JSON
    AVRO = DestinationFormat.AVRO
    PARQUET = DestinationFormat.PARQUET

    @staticmethod
    def to_extension(value):
        extensions = {
            FileFormat.CSV: 'csv',
            FileFormat.JSON: 'json',
            FileFormat.AVRO: 'avro',
            FileFormat.PARQUET: 'parquet'
        }
        return extensions[value]

    @staticmethod
    def to_format(value):
        value = value.upper()
        formats = {
            'CSV': FileFormat.CSV,
            'JSON': FileFormat.JSON,
            'AVRO': FileFormat.AVRO,
            'PARQUET': FileFormat.PARQUET,
            'JSONNL': FileFormat.JSON,
            FileFormat.JSON: FileFormat.JSON
        }
        return formats[value]


class CacheResponseType(object):
    CACHE_INFO = 'cache_info'
    DATASET = 'dataset'