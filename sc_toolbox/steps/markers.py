import scanpy as sc

from sc_toolbox.pipeline import Step


def _run(adata, params):
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method=params.get("method", "wilcoxon"),
        n_genes=int(params.get("n_genes", 50)),
        use_raw=False,
    )
    return adata


step = Step(
    name="markers",
    label="Marker genes (rank_genes_groups)",
    auto=True,
    run=_run,
)
