from copy import deepcopy
import scanpy as sc


def cluster(model_results, method: str = 'leiden', copy: bool = False):
    adata = sc.AnnData(X=model_results.embeddings)

    if method.lower() == 'leiden':
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)

    if copy:
        model_results = deepcopy(model_results)
        model_results.obs['leiden'] = adata.obs['leiden'].values
        return model_results

    else:
        model_results.obs['leiden'] = adata.obs['leiden'].values
