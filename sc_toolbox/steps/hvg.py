import scanpy as sc

from sc_toolbox.pipeline import Step


def _run(adata, params):
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=int(params.get("n_top_genes", 2000)),
        flavor=params.get("flavor", "seurat"),
    )
    return adata


step = Step(name="hvg", label="Highly variable genes", auto=True, run=_run)
