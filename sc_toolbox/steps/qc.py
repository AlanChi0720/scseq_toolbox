import scanpy as sc
import streamlit as st

from sc_toolbox import plots
from sc_toolbox.pipeline import Step


def _ensure_qc_metrics(adata, mito_prefix: str) -> None:
    adata.uns.setdefault("sc_toolbox", {})["mito_prefix"] = mito_prefix
    if "n_genes_by_counts" in adata.obs.columns:
        return
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None
    )


def _render(adata, defaults):
    mito_prefix = st.session_state.get("mito_prefix", "MT-")
    _ensure_qc_metrics(adata, mito_prefix)

    st.pyplot(plots.violin_qc(adata))

    d = defaults.get("qc", {})
    c1, c2 = st.columns(2)
    min_genes = c1.number_input("min_genes", value=int(d.get("min_genes", 200)), step=50, min_value=0)
    max_genes = c2.number_input("max_genes", value=int(d.get("max_genes", 6000)), step=100, min_value=1)
    min_counts = c1.number_input("min_counts", value=int(d.get("min_counts", 500)), step=100, min_value=0)
    max_counts = c2.number_input("max_counts", value=int(d.get("max_counts", 25000)), step=500, min_value=1)
    max_pct_mt = st.slider(
        "max_pct_mt (%)", 0.0, 100.0, float(d.get("max_pct_mt", 15.0)), 0.5
    )

    mask = (
        (adata.obs.n_genes_by_counts >= min_genes)
        & (adata.obs.n_genes_by_counts <= max_genes)
        & (adata.obs.total_counts >= min_counts)
        & (adata.obs.total_counts <= max_counts)
        & (adata.obs.pct_counts_mt <= max_pct_mt)
    )
    st.info(
        f"Cells passing: {int(mask.sum()):,} / {len(mask):,} "
        f"({100 * float(mask.mean()):.1f}%)"
    )

    return {
        "min_genes": int(min_genes),
        "max_genes": int(max_genes),
        "min_counts": int(min_counts),
        "max_counts": int(max_counts),
        "max_pct_mt": float(max_pct_mt),
    }


def _run(adata, params):
    mask = (
        (adata.obs.n_genes_by_counts >= params["min_genes"])
        & (adata.obs.n_genes_by_counts <= params["max_genes"])
        & (adata.obs.total_counts >= params["min_counts"])
        & (adata.obs.total_counts <= params["max_counts"])
        & (adata.obs.pct_counts_mt <= params["max_pct_mt"])
    )
    return adata[mask].copy()


step = Step(
    name="qc", label="Quality control", auto=False, run=_run, render=_render
)
