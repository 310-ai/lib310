def test_bigquery_select():
    from lib310.data.client import Client
    import dotenv
    dotenv.load_dotenv('.env')
    c = Client()
    query_job = c.query('SELECT * FROM `pfsdb3.0_uniprot.uniref` LIMIT 10')
    results = query_job.result()
    assert results.total_rows == 10