import scanpy as sc

from sc_toolbox.pipeline import Step


def _run(adata, params):
    n_pcs = int(adata.uns.get("sc_toolbox", {}).get("n_pcs", 20))
    sc.pp.neighbors(
        adata,
        n_neighbors=int(params.get("n_neighbors", 15)),
        n_pcs=n_pcs,
    )
    return adata


step = Step(name="neighbors", label="Neighbors graph", auto=True, run=_run)
