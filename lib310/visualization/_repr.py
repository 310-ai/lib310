from typing import Optional

import scanpy as sc


def _construct_adata(model_results, color: Optional[str] = None, **kwargs):
    adata = sc.AnnData(X=model_results.embeddings)

    if color:
        adata.obs['color'] = model_results.obs[color].values
    else:
        pass

    return adata


def umap(model_results, color: Optional[str] = None, **kwargs):
    adata = _construct_adata(model_results, color)

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    sc.pl.umap(adata, color='color', frameon=False, **kwargs)


def pca(model_results, color: Optional[str] = None, **kwargs):
    adata = _construct_adata(model_results, color)
    sc.pp.pca(adata, n_comps=50)
    sc.pl.pca(adata, color='color', frameon=False, **kwargs)


def tsne(model_results, color: Optional[str] = None, **kwargs):
    adata = _construct_adata(model_results, color)

    sc.pp.neighbors(adata)
    sc.tl.tsne(adata)

    sc.pl.tsne(adata, color='color', frameon=False, **kwargs)
