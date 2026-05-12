# sc_toolbox — Design Document (v1)

> Status: **Decisions captured. Ready to scaffold once user confirms.**
> Last updated: 2026-05-11

---

## 1. Goals

Build an interactive Python scRNA-seq pipeline runner with an "upload and wait" user experience. The user configures a pipeline, uploads data, then the program runs automatically — pausing only at human-decision points (QC thresholds, n_pcs, clustering resolution) where it shows the relevant plot and waits for the user's choice.

Conceptually: a Python equivalent of [scDrake](https://bioinfocz.github.io/scdrake/), but driven by a Streamlit web UI instead of a YAML+R config.

## 2. Non-goals (v1)

- Multi-sample integration / batch correction → v2
- Differential expression analysis → v2
- Cell-type annotation / marker discovery → v2
- Cloud deployment, multi-user → out of scope (local Streamlit only)
- Spatial transcriptomics → out of scope

## 3. User flow

```
┌──────────────────────────────────────────────────────────────┐
│ Screen 0 — Configure pipeline                                │
│   • Pick preset (standard_10x) or "custom"                   │
│   • Toggle optional steps (doublet, cell-cycle, …)           │
│   • [Next →]                                                 │
└──────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────┐
│ Screen 1 — Upload data                                       │
│   • Instructions dynamically reflect what the pipeline needs │
│       e.g. "Upload one of: 10x mtx (.tar.gz/.zip) or .h5ad"  │
│   • Sample name input                                        │
│   • [Run pipeline →]                                         │
└──────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────┐
│ Screen 2 — Run                                               │
│   Steps execute sequentially:                                │
│     load            → auto                                   │
│     qc              → ❚❚ pause: violin + sliders             │
│     filter          → auto                                   │
│     (doublet)       → auto (if enabled)                      │
│     normalize       → auto                                   │
│     hvg             → auto                                   │
│     pca             → ❚❚ pause: elbow + n_pcs slider         │
│     neighbors       → auto                                   │
│     cluster         → ❚❚ pause: UMAP previews + resolution   │
│     umap            → auto                                   │
└──────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────┐
│ Screen 3 — Done                                              │
│   • Final UMAP coloured by cluster                           │
│   • Summary stats (cells, genes, clusters)                   │
│   • Downloads: result.h5ad, params.json, report.html         │
└──────────────────────────────────────────────────────────────┘
```

Standard preset has **three** pause points. With Scrublet enabled, still three (Scrublet is auto).

## 4. UI mockups

### 4.1 Screen 0 — Pipeline config

```
┌─────────────────────────────────────────────┐
│ sc_toolbox                                  │
│ ─────────────────────────────────────────── │
│  Step 1 / 3 — Configure pipeline            │
│                                             │
│  Preset:  (•) standard_10x   ( ) custom    │
│                                             │
│  ☑ Quality control + filtering              │
│  ☐ Doublet detection (Scrublet)             │
│  ☐ Cell cycle scoring                       │
│  ☑ Normalization (log1p)                    │
│  ☑ Highly variable genes (HVG)              │
│  ☑ PCA                                      │
│  ☑ Neighbors + Leiden + UMAP                │
│                                             │
│  ▼ Advanced settings (initial defaults)     │
│    QC initial sliders:                      │
│      min_genes [200]   max_genes [6000]     │
│      min_counts[500]   max_counts[25000]    │
│      max_pct_mt [15]                        │
│    HVG: n_top_genes [2000]  flavor [seurat] │
│    PCA initial n_pcs [20]                   │
│    Neighbors: n_neighbors [15]              │
│    Cluster initial resolution [0.5]         │
│                                             │
│                          [ Next → ]         │
└─────────────────────────────────────────────┘
```

The advanced settings expander is collapsed by default. Users who don't open it get the standard defaults. Values set here populate (a) auto-step params, and (b) the initial position of pause-screen sliders — the user can still adjust at runtime.

### 4.2 Screen 1 — Upload

```
┌─────────────────────────────────────────────┐
│ Step 2 / 3 — Upload data                    │
│                                             │
│  Based on your pipeline, please upload:     │
│    • 10x mtx folder (.tar.gz / .zip), or    │
│    • .h5ad file, or                         │
│    • 10x .h5 file                           │
│                                             │
│  [ Drop file here or click to browse ]      │
│                                             │
│  Sample name: [pbmc_3k                  ]   │
│                                             │
│  Species (for mito gene prefix):            │
│    (•) Human   ( MT- )                      │
│    ( ) Mouse   ( mt- )                      │
│    ( ) Custom  [ ____  ]                    │
│                                             │
│         [ ← Back ]   [ Run pipeline → ]     │
└─────────────────────────────────────────────┘
```

### 4.3 Screen 2 — QC pause

```
┌─────────────────────────────────────────────┐
│ Running pipeline                            │
│ Load ✓  →  QC ❚❚ (waiting for you)          │
│ ─────────────────────────────────────────── │
│                                             │
│   [Violin: n_genes_by_counts]               │
│   [Violin: total_counts]                    │
│   [Violin: pct_counts_mt]                   │
│                                             │
│   n_genes_by_counts:                        │
│     [ 200 ] ─●──────────● [ 6000 ]          │
│   total_counts:                             │
│     [ 500 ] ─●──────────● [ 25000 ]         │
│   pct_counts_mt (max):                      │
│                          ●─── [ 15 % ]      │
│                                             │
│   Cells passing: 2,541 / 2,700 (94.1%)      │
│                                             │
│   [ Preview filter ]   [ Continue → ]       │
└─────────────────────────────────────────────┘
```

### 4.4 Screen 2 — PCA pause

```
┌─────────────────────────────────────────────┐
│ Running pipeline                            │
│ Load ✓ → QC ✓ → Filter ✓ → Norm ✓ → HVG ✓   │
│ → PCA ❚❚                                    │
│ ─────────────────────────────────────────── │
│                                             │
│   [Elbow plot: variance ratio for PCs 1-50] │
│                                             │
│   n_pcs: [ 20 ] ●────────                   │
│                                             │
│                       [ Continue → ]        │
└─────────────────────────────────────────────┘
```

### 4.5 Screen 2 — Cluster pause

```
┌─────────────────────────────────────────────┐
│ Cluster — choose resolution                 │
│                                             │
│   ┌────────┐ ┌────────┐ ┌────────┐ ┌──────┐ │
│   │UMAP@0.2│ │UMAP@0.5│ │UMAP@0.8│ │@1.2  │ │
│   │ 6 clus │ │ 9 clus │ │13 clus │ │17    │ │
│   └────────┘ └────────┘ └────────┘ └──────┘ │
│                                             │
│   Resolution: [ 0.5 ] ●─────                │
│   Algorithm:  (•) Leiden   ( ) Louvain      │
│                                             │
│                       [ Continue → ]        │
└─────────────────────────────────────────────┘
```

## 5. Architecture

### 5.1 Core abstractions

```python
# sc_toolbox/pipeline.py
from dataclasses import dataclass, field
from typing import Callable, Protocol
import anndata as ad

class RenderFn(Protocol):
    def __call__(self, adata: ad.AnnData, prior_params: dict) -> dict: ...
    # Renders Streamlit widgets/plots. Returns chosen params dict.

class RunFn(Protocol):
    def __call__(self, adata: ad.AnnData, params: dict) -> ad.AnnData: ...
    # Performs the actual scanpy work.

@dataclass
class Step:
    name: str               # stable id, e.g. "qc"
    label: str              # human label, e.g. "Quality Control"
    auto: bool              # True → no UI, run with defaults
    run: RunFn
    render: RenderFn | None = None   # required when auto=False
    needs: list[str] = field(default_factory=list)
        # declared input requirements; the upload screen aggregates these
        # e.g. ["counts_matrix"], or v2 ["counts_matrix", "sample_metadata"]

@dataclass
class Pipeline:
    name: str
    steps: list[Step]
    enabled: dict[str, bool] = field(default_factory=dict)

    def active_steps(self) -> list[Step]:
        return [s for s in self.steps if self.enabled.get(s.name, True)]

    def aggregated_needs(self) -> set[str]:
        out: set[str] = set()
        for s in self.active_steps():
            out.update(s.needs)
        return out
```

### 5.2 Streamlit session_state schema

```python
st.session_state = {
    "screen": "config" | "upload" | "run" | "done",
    "pipeline": Pipeline,                  # set on Screen 0 → Next
    "defaults": dict[str, dict],           # step_name → user-overridden defaults
                                           # from Screen 0's "Advanced settings"
    "sample_name": str,                    # set on Screen 1
    "mito_prefix": str,                    # "MT-" / "mt-" / user-supplied
    "step_idx": int,                       # index into pipeline.active_steps()
    "step_params": dict[str, dict],        # step_name → chosen params at runtime
    "adata": AnnData | None,
    "log": list[str],                      # human-readable progress log
}
```

### 5.3 Main loop (app.py, sketch)

```python
def main():
    s = st.session_state
    s.setdefault("screen", "config")

    {
        "config": screen_config,
        "upload": screen_upload,
        "run":    screen_run,
        "done":   screen_done,
    }[s.screen]()

def screen_run():
    steps = st.session_state.pipeline.active_steps()
    idx = st.session_state.step_idx

    if idx >= len(steps):
        st.session_state.screen = "done"
        st.rerun()
        return

    step = steps[idx]
    render_progress_bar(steps, idx)

    if step.auto:
        with st.spinner(f"Running {step.label}…"):
            st.session_state.adata = step.run(st.session_state.adata, {})
        st.session_state.step_idx += 1
        st.rerun()
    else:
        params = step.render(st.session_state.adata, st.session_state.step_params)
        if st.button("Continue →", key=f"continue_{step.name}"):
            with st.spinner(f"Applying {step.label}…"):
                st.session_state.adata = step.run(st.session_state.adata, params)
            st.session_state.step_params[step.name] = params
            st.session_state.step_idx += 1
            st.rerun()
```

## 6. Step specifications

Each step lists: AnnData input/output contract, scanpy functions used, UI controls (if any), and parameters returned.

### 6.1 `load` (auto)
- **Input:** uploaded archive/file from Screen 1
- **Output AnnData:** raw counts in `adata.X`; `adata.var_names` unique
- **Backend:** `sc.read_10x_mtx`, `sc.read_10x_h5`, `sc.read_h5ad`
- **Auto behaviour:** detect format by extension/magic bytes

### 6.2 `qc` (pause)
- **Pre-compute (in render):**
  - `adata.var["mt"] = adata.var_names.str.startswith(st.session_state.mito_prefix)`
    — prefix comes from the user's Screen 1 selection, NOT auto-detected.
  - `sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=False)`
- **UI:**
  - 3 violin plots (one each: n_genes_by_counts, total_counts, pct_counts_mt)
  - Range sliders for n_genes, total_counts; single slider (upper bound) for pct_counts_mt
  - Initial slider positions come from `st.session_state.defaults["qc"]` (set on Screen 0)
  - Live counter: "Cells passing: X / Y"
- **Params returned:** `{min_genes, max_genes, min_counts, max_counts, max_pct_mt}`
- **Run:** boolean-mask `adata` accordingly

### 6.3 `doublet` (optional, auto, off by default)
- **Backend:** Scrublet
- Adds `adata.obs["doublet_score"]`, `adata.obs["predicted_doublet"]`
- Removes predicted doublets

### 6.4 `normalize` (auto)
- `adata.layers["counts"] = adata.X.copy()` (preserve raw)
- `sc.pp.normalize_total(adata, target_sum=1e4)`
- `sc.pp.log1p(adata)`

### 6.5 `hvg` (auto)
- `sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")`
- **Does not subset** — only flags `adata.var["highly_variable"]` so downstream still has full gene info

### 6.6 `pca` (pause)
- **Pre-compute (in render):** `sc.tl.pca(adata, n_comps=50, use_highly_variable=True, mask_var="highly_variable")`
- **UI:** matplotlib elbow plot of `adata.uns["pca"]["variance_ratio"][:50]` + integer slider for `n_pcs` (default 20, range 5–50)
- **Params returned:** `{n_pcs}`
- **Run:** persist choice to `adata.uns["sc_toolbox"]["n_pcs"]`

### 6.7 `neighbors` (auto)
- `sc.pp.neighbors(adata, n_neighbors=15, n_pcs=adata.uns["sc_toolbox"]["n_pcs"])`

### 6.8 `cluster` (pause)
- **Pre-compute (in render):** for each `res ∈ [0.2, 0.5, 0.8, 1.2]`, run `sc.tl.leiden(adata, resolution=res, key_added=f"_preview_{res}")`. Run a temporary `sc.tl.umap(adata)` once and cache coords for previews.
- **UI:** 4 small UMAP scatter plots side by side, each labelled with cluster count; final slider for chosen resolution (continuous, default 0.5); radio for algorithm (Leiden / Louvain).
- **Params returned:** `{resolution, algorithm}`
- **Run:** clean up `_preview_*` columns, run the final `sc.tl.leiden(adata, resolution=params["resolution"], key_added="leiden")`

### 6.9 `umap` (auto)
- `sc.tl.umap(adata)` (final, with default min_dist/spread)

## 7. Done screen

- Final UMAP coloured by `leiden`, downloadable as PNG
- Summary table: n_cells (before/after QC), n_genes, n_HVG, n_pcs chosen, n_clusters
- Downloads:
  - **`result.h5ad`** — full AnnData
  - **`params.json`** — reproducibility manifest
  - **`report.html`** — static HTML with embedded plots (a simple Jinja template; PDF deferred to v1.1)

## 8. Reproducibility — `params.json`

```json
{
  "sc_toolbox_version": "0.1.0",
  "timestamp": "2026-05-11T14:30:00Z",
  "pipeline": "standard_10x",
  "enabled_steps": ["load", "qc", "filter", "normalize", "hvg",
                    "pca", "neighbors", "cluster", "umap"],
  "step_params": {
    "load":      {"format": "10x_mtx", "filename": "pbmc_3k.tar.gz"},
    "qc":        {"min_genes": 200, "max_genes": 6000,
                  "min_counts": 500, "max_counts": 25000,
                  "max_pct_mt": 15.0},
    "normalize": {"target_sum": 10000},
    "hvg":       {"n_top_genes": 2000, "flavor": "seurat"},
    "pca":       {"n_pcs": 20},
    "neighbors": {"n_neighbors": 15},
    "cluster":   {"resolution": 0.5, "algorithm": "leiden"}
  }
}
```

v1.1 will add a headless replay mode: `sc-toolbox replay params.json --input pbmc_3k.tar.gz`.

## 9. Final file layout

```
sc_toolbox/
├── pyproject.toml
├── .venv/                   ← created by `python -m venv .venv`
├── README.md
├── DESIGN.md                ← this file
├── app.py                   ← `streamlit run app.py`
├── sc_toolbox/
│   ├── __init__.py
│   ├── pipeline.py          ← Step, Pipeline, run-loop driver
│   ├── presets.py           ← standard_10x() factory
│   ├── plots.py             ← violin / elbow / UMAP helpers (return mpl Figures)
│   ├── io_utils.py          ← load archives, detect format
│   └── steps/
│       ├── __init__.py
│       ├── load.py
│       ├── qc.py
│       ├── doublet.py
│       ├── normalize.py
│       ├── hvg.py
│       ├── pca.py
│       ├── neighbors.py
│       ├── cluster.py
│       └── umap.py
└── tests/
    ├── test_pipeline.py
    └── test_steps_smoke.py  ← tiny synthetic AnnData fixture, runs full pipeline
```

## 10. Dependencies (pyproject.toml)

```toml
[project]
name = "sc-toolbox"
version = "0.1.0"
requires-python = ">=3.11,<3.13"
dependencies = [
    "streamlit>=1.30",
    "scanpy>=1.10",
    "anndata>=0.10",
    "leidenalg>=0.10",
    "louvain>=0.8",
    "python-igraph>=0.11",
    "scrublet>=0.2.3",
    "matplotlib>=3.8",
    "seaborn>=0.13",
    "pandas>=2.0",
    "numpy>=1.26",
    "scipy>=1.11",
    "h5py>=3.10",
    "jinja2>=3.1",
]

[project.optional-dependencies]
dev = ["pytest>=8.0", "ruff>=0.5"]
```

## 11. Decisions log

Resolved on 2026-05-11:

1. **Mitochondrial gene prefix → user selects on Screen 1.** Dropdown with Human (`MT-`) / Mouse (`mt-`) / Custom (text input). No auto-detection.
2. **No checkpoints in v1.** Reason: PBMC-3k-scale runs would produce ~9× the AnnData size if each step were saved, which bloats disk and adds little value. AnnData lives only in `st.session_state` (in-memory). If the browser refreshes or the Streamlit server restarts, the user reruns from the beginning. Re-evaluate if users frequently lose progress on long runs.
3. **Report format → HTML only in v1.** PDF export deferred to v1.1.
4. **Defaults → exposed via "Advanced settings" expander on Screen 0.** Standard values prefilled (`min_genes=200, max_pct_mt=15, n_top_genes=2000, n_pcs=20, resolution=0.5`); users can override in one place before the pipeline runs. These values feed both auto-step parameters and the initial position of pause-screen sliders.

## 12. Definition of done (v1)

- [ ] `python -m venv .venv && .venv\Scripts\Activate.ps1 && pip install -e ".[dev]"` works on fresh Python 3.11 env
- [ ] `streamlit run app.py` opens browser at `localhost:8501`
- [ ] PBMC 3k example data runs end-to-end through the standard preset
- [ ] All three pause screens (QC / PCA / cluster) render plots and accept input correctly
- [ ] Downloads on done screen produce valid `.h5ad` + `params.json` + `report.html`
- [ ] Smoke test in `tests/test_steps_smoke.py` runs full pipeline on synthetic data in <30s
