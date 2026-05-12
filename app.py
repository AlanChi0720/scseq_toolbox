"""sc_toolbox — Streamlit entry point.

Run with:  streamlit run app.py
"""
import json
import os
import tempfile
from datetime import datetime, timezone

import streamlit as st

from sc_toolbox import __version__
from sc_toolbox.io_utils import load_uploaded
from sc_toolbox.plots import markers_dotplot, umap_preview
from sc_toolbox.presets import PRESETS, standard_10x


# ----------------------------------------------------------------------------
# Session state
# ----------------------------------------------------------------------------

def init_state() -> None:
    s = st.session_state
    s.setdefault("screen", "config")
    s.setdefault("pipeline", None)
    s.setdefault("defaults", {})
    s.setdefault("sample_name", "")
    s.setdefault("mito_prefix", "MT-")
    s.setdefault("step_idx", 0)
    s.setdefault("step_params", {})
    s.setdefault("adata", None)


# ----------------------------------------------------------------------------
# Screen 0 — configure pipeline
# ----------------------------------------------------------------------------

def screen_config() -> None:
    st.title("sc_toolbox")
    st.caption(f"v{__version__} · Step 1 / 3 — Configure pipeline")

    preset_name = st.radio(
        "Preset",
        list(PRESETS.keys()) + ["custom"],
        horizontal=True,
        index=0,
    )
    base = PRESETS[preset_name]() if preset_name in PRESETS else standard_10x()

    st.subheader("Steps")
    enabled: dict[str, bool] = {}
    for step in base.steps:
        if step.name == "load":
            enabled["load"] = True
            continue
        enabled[step.name] = st.checkbox(
            step.label,
            value=base.enabled.get(step.name, True),
            key=f"en_{step.name}",
        )

    with st.expander("Advanced settings — initial defaults"):
        defaults: dict[str, dict] = {}

        st.markdown("**QC initial slider positions**")
        c1, c2 = st.columns(2)
        defaults["qc"] = {
            "min_genes": c1.number_input("min_genes", value=200, step=50, min_value=0),
            "max_genes": c2.number_input("max_genes", value=6000, step=100, min_value=1),
            "min_counts": c1.number_input("min_counts", value=500, step=100, min_value=0),
            "max_counts": c2.number_input("max_counts", value=25000, step=500, min_value=1),
            "max_pct_mt": st.number_input("max_pct_mt (%)", value=15.0, step=0.5, min_value=0.0, max_value=100.0),
        }

        st.markdown("**HVG**")
        defaults["hvg"] = {
            "n_top_genes": st.number_input("n_top_genes", value=2000, step=100, min_value=100),
            "flavor": st.selectbox("flavor", ["seurat", "cell_ranger", "seurat_v3"], index=0),
        }

        st.markdown("**PCA**")
        defaults["pca"] = {
            "n_pcs": st.number_input("n_pcs (initial)", value=20, step=1, min_value=5, max_value=50)
        }

        st.markdown("**Neighbors**")
        defaults["neighbors"] = {
            "n_neighbors": st.number_input("n_neighbors", value=15, step=1, min_value=2)
        }

        st.markdown("**Cluster**")
        defaults["cluster"] = {
            "resolution": st.number_input(
                "resolution (initial)", value=0.5, step=0.05, min_value=0.05, max_value=2.0
            )
        }

        st.markdown("**Markers**")
        defaults["markers"] = {
            "method": st.selectbox(
                "method", ["wilcoxon", "t-test", "logreg"], index=0
            ),
            "n_genes": st.number_input(
                "n_genes (per cluster)", value=50, step=10, min_value=5, max_value=500
            ),
        }

    if st.button("Next →", type="primary"):
        base.enabled = enabled
        st.session_state.pipeline = base
        st.session_state.defaults = defaults
        st.session_state.screen = "upload"
        st.rerun()


# ----------------------------------------------------------------------------
# Screen 1 — upload data
# ----------------------------------------------------------------------------

def screen_upload() -> None:
    st.title("sc_toolbox")
    st.caption(f"v{__version__} · Step 2 / 3 — Upload data")

    st.markdown(
        "Based on your pipeline, please upload **one** of:\n"
        "- 10x mtx folder packaged as `.tar.gz` / `.tgz` / `.zip`\n"
        "- `.h5ad` file\n"
        "- 10x `.h5` file"
    )

    uploaded = st.file_uploader(
        "Upload data",
        type=["h5ad", "h5", "tar.gz", "tgz", "zip"],
        accept_multiple_files=False,
    )

    sample_name = st.text_input(
        "Sample name",
        value=st.session_state.sample_name or "sample_1",
    )

    st.subheader("Species (for mito-gene prefix)")
    species = st.radio(
        "Species",
        ["Human (MT-)", "Mouse (mt-)", "Custom"],
        index=0,
        label_visibility="collapsed",
    )
    custom_prefix = ""
    if species == "Custom":
        custom_prefix = st.text_input("Custom mito prefix", value="MT-")

    c1, c2 = st.columns([1, 1])
    if c1.button("← Back"):
        st.session_state.screen = "config"
        st.rerun()
    if c2.button("Run pipeline →", type="primary", disabled=uploaded is None):
        with st.spinner("Loading data…"):
            adata = load_uploaded(uploaded)
        st.session_state.adata = adata
        st.session_state.sample_name = sample_name
        if species.startswith("Human"):
            st.session_state.mito_prefix = "MT-"
        elif species.startswith("Mouse"):
            st.session_state.mito_prefix = "mt-"
        else:
            st.session_state.mito_prefix = custom_prefix or "MT-"
        st.session_state.step_idx = 0
        st.session_state.step_params = {}
        st.session_state.screen = "run"
        st.rerun()


# ----------------------------------------------------------------------------
# Screen 2 — run pipeline
# ----------------------------------------------------------------------------

def _render_progress(steps, idx) -> None:
    crumbs = []
    for i, step in enumerate(steps):
        if i < idx:
            crumbs.append(f"{step.label} ✓")
        elif i == idx:
            crumbs.append(f"**{step.label}** ⏳" if step.auto else f"**{step.label}** ⏸")
        else:
            crumbs.append(step.label)
    st.caption(" → ".join(crumbs))


def screen_run() -> None:
    st.title("sc_toolbox")
    st.caption(f"v{__version__} · Step 3 / 3 — Run pipeline")

    s = st.session_state
    steps = s.pipeline.active_steps()

    if s.step_idx >= len(steps):
        s.screen = "done"
        st.rerun()
        return

    _render_progress(steps, s.step_idx)
    step = steps[s.step_idx]

    if step.auto:
        with st.spinner(f"Running {step.label}…"):
            defaults = s.defaults.get(step.name, {})
            s.adata = step.run(s.adata, defaults)
        s.step_params[step.name] = {}
        s.step_idx += 1
        st.rerun()
        return

    st.subheader(step.label)
    params = step.render(s.adata, s.defaults)
    if st.button("Continue →", type="primary", key=f"continue_{step.name}"):
        with st.spinner(f"Applying {step.label}…"):
            s.adata = step.run(s.adata, params)
        s.step_params[step.name] = params
        # Invalidate cluster preview cache if we just left a step that mutates obs
        s.pop("_cluster_previews_built", None)
        s.step_idx += 1
        st.rerun()


# ----------------------------------------------------------------------------
# Screen 3 — done
# ----------------------------------------------------------------------------

def _adata_to_bytes(adata) -> bytes:
    fd, path = tempfile.mkstemp(suffix=".h5ad")
    os.close(fd)
    try:
        adata.write_h5ad(path)
        with open(path, "rb") as f:
            return f.read()
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def screen_done() -> None:
    st.title("sc_toolbox — Done 🎉")
    s = st.session_state
    adata = s.adata

    st.success(f"Pipeline complete on sample `{s.sample_name}`")

    c1, c2, c3 = st.columns(3)
    c1.metric("Cells", adata.n_obs)
    c2.metric("Genes", adata.n_vars)
    if "leiden" in adata.obs.columns:
        c3.metric("Clusters", int(adata.obs.leiden.nunique()))

    if "X_umap" in adata.obsm and "leiden" in adata.obs.columns:
        st.pyplot(umap_preview(adata, color_key="leiden", title="Final UMAP"))

    if "rank_genes_groups" in adata.uns:
        st.subheader("Top markers per cluster")
        st.pyplot(markers_dotplot(adata, n_genes=5))

    params_payload = {
        "sc_toolbox_version": __version__,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "pipeline": s.pipeline.name,
        "enabled_steps": [step.name for step in s.pipeline.active_steps()],
        "step_params": s.step_params,
        "sample_name": s.sample_name,
        "mito_prefix": s.mito_prefix,
    }
    st.download_button(
        "Download params.json",
        json.dumps(params_payload, indent=2, default=str).encode("utf-8"),
        file_name="params.json",
        mime="application/json",
    )

    st.download_button(
        "Download result.h5ad",
        _adata_to_bytes(adata),
        file_name="result.h5ad",
        mime="application/octet-stream",
    )

    if st.button("Start a new run"):
        for k in list(s.keys()):
            del s[k]
        st.rerun()


# ----------------------------------------------------------------------------
# Router
# ----------------------------------------------------------------------------

SCREENS = {
    "config": screen_config,
    "upload": screen_upload,
    "run": screen_run,
    "done": screen_done,
}


def main() -> None:
    st.set_page_config(page_title="sc_toolbox", layout="wide")
    init_state()
    SCREENS[st.session_state.screen]()


if __name__ == "__main__":
    main()
