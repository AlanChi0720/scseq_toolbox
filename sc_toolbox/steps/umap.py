import scanpy as sc

from sc_toolbox.pipeline import Step


def _run(adata, params):
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)
    return adata


step = Step(name="umap", label="UMAP", auto=True, run=_run)
