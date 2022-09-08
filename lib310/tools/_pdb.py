import requests

BASE_URL = "https://files.rcsb.org/download"


def download(pdb_id, file_type, out_dir):
    if '_' in pdb_id:
        pdb_id = pdb_id[:-2]
    url = "{BASE_URL}/{pdb_id}.{type}.gz".format(pdb_id=pdb_id, BASE_URL=BASE_URL, type=file_type)
    res = requests.get(url)
    outfile = "{out_dir}/{pdb_id}.{type}.gz".format(pdb_id=pdb_id, type=file_type, out_dir=out_dir)
    open(outfile, "wb").write(res.content)


def download_list_pdb(pdb_list, type_list, out_dir):
    import os
    os.makedirs(out_dir, exist_ok=True)
    for pdb in pdb_list:
        for tp in type_list:
            download(pdb, tp, out_dir)

def get_pdb_parser(pdb_name,pdb_file):
    from Bio.PDB import PDBParser
    return PDBParser().get_structure(pdb_name, pdb_file)



# download_list_pdb(["7niu_A", "6s7p_A"], ["pdb", "pdb1", "xml"],"./data")
