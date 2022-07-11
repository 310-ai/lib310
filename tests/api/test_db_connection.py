def test_bigquery_select():
    import lib310
    import dotenv
    dotenv.load_dotenv('.env')
    seqs = lib310.db.fetch(query="SELECT * FROM `pfsdb3.0_go.gaf` LIMIT 10")
    print(seqs)
    pass