import requests

BASE_URL = "https://files.rcsb.org/download"


def download(pdb_id, type,outdir):
    if pdb_id.find('_'):
        pdb_id = pdb_id[:-2]
    url = "{BASE_URL}/{pdb_id}.{type}.gz".format(pdb_id=pdb_id, BASE_URL=BASE_URL, type=type)
    res = requests.get(url)
    outfile = "{outdir}/{pdb_id}.{type}.gz".format(pdb_id=pdb_id, type=type, outdir=outdir)
    open(outfile, "wb").write(res.content)


def download_list_pdb(pdb_list, type_list,outdir):
    import os
    os.mkdir(outdir)
    for pdb in pdb_list:
        for tp in type_list:
            download(pdb, tp,outdir)


# download_list_pdb(["7niu_A", "6s7p_A"], ["pdb", "pdb1", "xml"],"./data")
