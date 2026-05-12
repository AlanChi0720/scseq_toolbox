"""Smoke test: run the standard_10x pipeline end-to-end on synthetic data.

Bypasses Streamlit by calling each step's `run()` directly. Pause steps need
their pre-computations (qc metrics, PCA, cluster previews) seeded manually,
since those normally run inside the corresponding `render()` function.
"""
import anndata as ad
import numpy as np
import scanpy as sc

from sc_toolbox.presets import standard_10x
from sc_toolbox.steps.cluster import PREVIEW_RES, _preview_key


def make_synthetic_adata(n_cells: int = 500, n_genes: int = 1500) -> ad.AnnData:
    rng = np.random.default_rng(42)
    X = rng.poisson(0.5, size=(n_cells, n_genes)).astype("float32")
    gene_names = [f"GENE{i}" for i in range(n_genes - 20)] + [
        f"MT-{i}" for i in range(20)
    ]
    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.var_names_make_unique()
    return adata


def test_smoke_standard_pipeline():
    adata = make_synthetic_adata()
    pipeline = standard_10x()

    # Pre-compute QC metrics (lives inside qc.render in the live app)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None
    )

    runtime_params = {
        "qc": {
            "min_genes": 0,
            "max_genes": 10**6,
            "min_counts": 0,
            "max_counts": 10**8,
            "max_pct_mt": 100.0,
        },
        "pca": {"n_pcs": 10},
        "cluster": {"resolution": 0.5, "algorithm": "leiden"},
    }

    for step in pipeline.active_steps():
        # Seed pre-computations that the live UI does inside render()
        if step.name == "pca" and "X_pca" not in adata.obsm:
            sc.tl.pca(adata, n_comps=min(50, min(adata.shape) - 1), use_highly_variable=True)
        if step.name == "cluster":
            if "X_umap" not in adata.obsm:
                sc.tl.umap(adata)
            for res in PREVIEW_RES:
                sc.tl.leiden(adata, resolution=res, key_added=_preview_key(res))

        adata = step.run(adata, runtime_params.get(step.name, {}))

    assert adata.n_obs > 0
    assert "X_umap" in adata.obsm
    assert "leiden" in adata.obs.columns
    # preview cols should have been cleaned up by cluster.run
    for res in PREVIEW_RES:
        assert _preview_key(res) not in adata.obs.columns
    # markers step populated rank_genes_groups in uns
    assert "rank_genes_groups" in adata.uns
    assert "names" in adata.uns["rank_genes_groups"]
