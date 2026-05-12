import scanpy as sc

from sc_toolbox.pipeline import Step


def _run(adata, params):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=float(params.get("target_sum", 1e4)))
    sc.pp.log1p(adata)
    return adata


step = Step(name="normalize", label="Normalization (log1p)", auto=True, run=_run)
