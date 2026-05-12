import scanpy as sc
import streamlit as st

from sc_toolbox import plots
from sc_toolbox.pipeline import Step


def _render(adata, defaults):
    if "X_pca" not in adata.obsm:
        n_comps = min(50, min(adata.shape) - 1)
        sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
    st.pyplot(plots.elbow_pca(adata))

    d = defaults.get("pca", {})
    max_pcs = int(adata.uns["pca"]["variance_ratio"].shape[0])
    n_pcs = st.slider(
        "n_pcs", min_value=5, max_value=max_pcs, value=int(d.get("n_pcs", 20))
    )
    return {"n_pcs": int(n_pcs)}


def _run(adata, params):
    adata.uns.setdefault("sc_toolbox", {})["n_pcs"] = int(params["n_pcs"])
    return adata


step = Step(name="pca", label="PCA", auto=False, run=_run, render=_render)
