"""Cluster annotation step.

Paused step that shows top markers per cluster as a dot plot and lets the
user type a cell-type name for each cluster. Blank inputs fall back to
`cluster_<id>` so the resulting `adata.obs["cell_type"]` column is always
fully populated.
"""
import streamlit as st

from sc_toolbox import plots
from sc_toolbox.pipeline import Step


def _render(adata, defaults):
    if "rank_genes_groups" not in adata.uns:
        st.warning(
            "No marker gene results in adata.uns — make sure the 'Marker genes' "
            "step is enabled before this one. Continuing will assign default "
            "cluster_<id> labels."
        )
        labels = {c: "" for c in adata.obs["leiden"].cat.categories}
        return {"labels": labels}

    st.pyplot(plots.markers_dotplot(adata, n_genes=5))

    prior = defaults.get("annotate", {})
    labels: dict[str, str] = {}
    clusters = list(adata.obs["leiden"].cat.categories)

    st.write("Type a cell-type name for each cluster (blank → `cluster_<id>`).")
    cols = st.columns(2)
    for i, cluster in enumerate(clusters):
        with cols[i % 2]:
            labels[cluster] = st.text_input(
                f"Cluster {cluster}",
                value=str(prior.get(cluster, "")),
                key=f"annotate_{cluster}",
            ).strip()

    return {"labels": labels}


def _run(adata, params):
    labels = params.get("labels", {}) or {}
    mapping = {
        c: (labels.get(c) or f"cluster_{c}")
        for c in adata.obs["leiden"].cat.categories
    }
    adata.obs["cell_type"] = adata.obs["leiden"].map(mapping).astype("category")
    return adata


step = Step(
    name="annotate",
    label="Cluster annotation",
    auto=False,
    run=_run,
    render=_render,
)
