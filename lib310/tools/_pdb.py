import requests

BASE_URL = "https://files.rcsb.org/download"


def download(pdb_id, file_type, out_dir):
    if '_' in pdb_id:
        pdb_id = pdb_id[:-2]
    url = "{BASE_URL}/{pdb_id}.{type}.gz".format(pdb_id=pdb_id, BASE_URL=BASE_URL, type=file_type)
    res = requests.get(url)
    outfile = "{out_dir}/{pdb_id}.{type}.gz".format(pdb_id=pdb_id, type=file_type, out_dir=out_dir)
    open(outfile, "wb").write(res.content)
    import gzip
    import shutil
    zipfile = outfile
    output_file = zipfile[:-3]
    with gzip.open(zipfile, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    import os
    os.remove(outfile)


def download_list_pdb(pdb_list, type_list, out_dir):
    import os
    os.makedirs(out_dir, exist_ok=True)
    for pdb in pdb_list:
        for tp in type_list:
            download(pdb, tp, out_dir)


def get_pdb_parser(pdb_name, pdb_file):
    from Bio.PDB import PDBParser
    return PDBParser().get_structure(pdb_name, pdb_file)


# download_list_pdb(["7niu_A", "6s7p_A"], ["pdb", "pdb1", "xml"],"./data")


def download_from_alphafold(source_blob_name, destination_file_name):
    from google.cloud import storage
    query = """with file_rows AS (
      with file_cols AS (
        SELECT
          CONCAT(entryID, '-model_v', latestVersion, '.cif') as m,
          CONCAT(entryID, '-predicted_aligned_error_v', latestVersion, '.json') as p
        FROM bigquery-public-data.deepmind_alphafold.metadata
        WHERE organismScientificName = "Homo sapiens"
          AND (fractionPlddtVeryHigh + fractionPlddtConfident) > 0.5
      )
      SELECT * FROM file_cols UNPIVOT (files for filetype in (m, p))
    )
    SELECT CONCAT('gs://public-datasets-deepmind-alphafold/', files) as files
    FROM file_rows LIMIT 10"""
    bucket_name = "public-datasets-deepmind-alphafold"
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)
    print(
        "Downloaded storage object {} from bucket {} to local file {}.".format(
            source_blob_name, bucket_name, destination_file_name
        )
    )
# download_blob("public-datasets-deepmind-alphafold","AF-J3QLM0-F1-model_v3.cif","a.cif")
