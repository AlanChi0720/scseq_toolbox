import scanpy as sc
import streamlit as st

from sc_toolbox import plots
from sc_toolbox.pipeline import Step

PREVIEW_RES = (0.2, 0.5, 0.8, 1.2)


def _preview_key(res: float) -> str:
    return f"_preview_leiden_{res}"


def _render(adata, defaults):
    cache_key = "_cluster_previews_built"
    if not st.session_state.get(cache_key):
        if "X_umap" not in adata.obsm:
            with st.spinner("Computing UMAP for preview…"):
                sc.tl.umap(adata)
        for res in PREVIEW_RES:
            key = _preview_key(res)
            if key not in adata.obs.columns:
                sc.tl.leiden(adata, resolution=res, key_added=key)
        st.session_state[cache_key] = True

    cols = st.columns(len(PREVIEW_RES))
    for col, res in zip(cols, PREVIEW_RES):
        key = _preview_key(res)
        n_clusters = adata.obs[key].nunique()
        with col:
            st.pyplot(
                plots.umap_preview(
                    adata, color_key=key, title=f"res={res} → {n_clusters} clusters"
                )
            )

    d = defaults.get("cluster", {})
    resolution = st.slider(
        "resolution", 0.05, 2.0, float(d.get("resolution", 0.5)), 0.05
    )
    algorithm = st.radio("algorithm", ["leiden", "louvain"], horizontal=True, index=0)
    return {"resolution": float(resolution), "algorithm": algorithm}


def _run(adata, params):
    for res in PREVIEW_RES:
        key = _preview_key(res)
        if key in adata.obs.columns:
            del adata.obs[key]
    if params.get("algorithm", "leiden") == "louvain":
        sc.tl.louvain(adata, resolution=float(params["resolution"]), key_added="leiden")
    else:
        sc.tl.leiden(adata, resolution=float(params["resolution"]), key_added="leiden")
    return adata


step = Step(name="cluster", label="Cluster", auto=False, run=_run, render=_render)
