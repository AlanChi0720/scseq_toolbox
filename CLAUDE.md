# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`sc_toolbox` is an interactive Python scRNA-seq pipeline runner (conceptually a Python port of [scDrake](https://bioinfocz.github.io/scdrake/)). A Streamlit web UI drives a configurable scanpy pipeline that auto-runs, pausing only at human-decision steps (QC thresholds, n_pcs, clustering resolution) where it draws the relevant plot and waits for user input. v1 is single-sample only.

The high-level UX contract — every design change must preserve this — is "upload and wait": opinionated defaults, minimal prompts, only stop for genuinely subjective choices. See [DESIGN.md](DESIGN.md) for the full design and decision log.

## Common commands

The project uses traditional `venv` + `pip` (not uv).

```powershell
# First-time install (Python 3.11)
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"                       # installs runtime + dev deps (pytest, ruff)

# Launch the Streamlit app (with venv activated)
streamlit run app.py
# …or without activating:
.venv\Scripts\python.exe -m streamlit run app.py

# Run tests
pytest -v tests/
pytest tests/test_steps_smoke.py::test_smoke_standard_pipeline   # single test
```

Streamlit auto-reloads on file save; if it doesn't, click "Always rerun" in the upper-right banner.

## Windows / Python quirks

- On Windows the bare `python` command on PATH may resolve to a Microsoft Store App-execution-alias stub at `%LOCALAPPDATA%\Microsoft\WindowsApps\python.exe` that exits with "Python was not found" — do **not** treat that as evidence Python is absent. Probe with `Get-Command python | ft Source`, check the `py` launcher's registered installs via `%LOCALAPPDATA%\Programs\Python\Launcher\py.exe --list`, and look at install dirs under `%LOCALAPPDATA%\Programs\Python\` before declaring anything missing.
- Prefer calling `.venv\Scripts\python.exe` by full path in shell commands to bypass the alias.

## Architecture

### The `Step` / `Pipeline` abstraction

[`sc_toolbox/pipeline.py`](sc_toolbox/pipeline.py) defines two dataclasses that everything else hangs off:

- `Step` has a `name`, `label`, `auto` flag, `run(adata, params)` function, and optional `render(adata, defaults)` function.
- `Pipeline` is an ordered list of `Step`s plus an `enabled` dict toggling each by name. `active_steps()` filters by that toggle.

**Critical split: `render` vs `run`.** Paused steps (`auto=False`) must implement both. `render` is called repeatedly while the user is on the pause screen (every Streamlit rerun triggered by widget interaction); it does the heavy pre-computation needed to draw plots (computing QC metrics, running PCA, building cluster previews) AND draws the Streamlit widgets, returning the chosen params dict. `run` is called once when the user clicks Continue and finalizes the step (applies the filter, persists the chosen `n_pcs`, runs the final clustering). Because of this, calling `step.run()` directly (e.g. from tests) requires manually seeding the state that `render` would have set up — the smoke test does exactly this.

Auto steps have `render=None` and just consume params dicts from `st.session_state.defaults[step.name]`.

### Streamlit driver loop

[`app.py`](app.py) is a four-screen state machine routed by `st.session_state.screen` ∈ {`config`, `upload`, `run`, `done`}. The flow is **pipeline-config FIRST, upload SECOND** — the upload screen's accepted formats and required side-data are determined by the configured pipeline. Never reorder.

`screen_run` is the core loop: read `step_idx`, get the next active step, dispatch on `step.auto`. After every step transition it pops `_cluster_previews_built` from session_state so the cluster preview cache is invalidated whenever `adata` mutates.

`st.session_state` schema (set up by `init_state` in app.py):
- `screen`, `step_idx` — routing
- `pipeline: Pipeline` — chosen on Screen 0
- `defaults: dict[step_name, dict]` — initial values from Screen 0's "Advanced settings" expander (feeds auto-step params AND initial slider positions on pause screens)
- `step_params: dict[step_name, dict]` — params actually chosen at runtime, written to `params.json` on the done screen
- `mito_prefix: str` — picked on Screen 1, read by `qc.py` (NOT auto-detected — explicit user choice)
- `adata: AnnData | None` — the single mutating state object threaded through every step

There are **no checkpoints in v1** by design. AnnData lives only in `session_state`; browser refresh = restart from beginning. Don't add disk-backed step persistence without re-discussing — the decision was driven by AnnData bloat (9 step copies of a PBMC-3k-sized dataset ≈ 1 GB).

### Adding a new step

1. Create `sc_toolbox/steps/foo.py` with at minimum a `_run(adata, params) -> AnnData` and a module-level `step = Step(name="foo", label="…", auto=True/False, run=_run, render=…)`.
2. If `auto=False`, implement `_render(adata, defaults) -> dict` that draws Streamlit widgets and returns params. Do all expensive prep inside `_render` (cache via `st.session_state` if needed — see the `_cluster_previews_built` flag pattern in [cluster.py](sc_toolbox/steps/cluster.py)).
3. Import and register the step in [`sc_toolbox/presets.py`](sc_toolbox/presets.py) `standard_10x()` and add an `enabled` entry.
4. If the step depends on outputs from a prior pause step, read them from `adata.uns["sc_toolbox"]` (the convention; see how `neighbors.py` reads `n_pcs`).

### scanpy state conventions

- Raw counts are stashed to `adata.layers["counts"]` in `normalize.py` before log-normalization.
- HVG selection only **flags** `adata.var["highly_variable"]`; it does not subset. Downstream steps that need HVG-only data should pass `use_highly_variable=True`.
- The mito-prefix mask in `qc.py` writes `adata.var["mt"]`; `pct_counts_mt` etc. come from `sc.pp.calculate_qc_metrics`.
- The cluster step writes the final clustering to `adata.obs["leiden"]` regardless of whether the user picked Leiden or Louvain (downstream code only knows the key name).

## Known gotchas

- **`adata.obs.pop(key, default)` does not work.** pandas `DataFrame.pop` takes only the column name (no default); it raises `TypeError: DataFrame.pop() takes 2 positional arguments but 3 were given`. Use `if key in adata.obs.columns: del adata.obs[key]`. This bit cluster.py once already.
- **Always install with the `[dev]` extra** (`pip install -e ".[dev]"`) so `pytest` and `ruff` come along. Plain `pip install -e .` gives runtime deps only and `pytest` will be missing.
- **Leiden FutureWarning** about backend changing to igraph is benign (scanpy ≥1.10 deprecation notice). To silence, pass `flavor="igraph", n_iterations=2, directed=False` to every `sc.tl.leiden` call.
- **`numpy<2.0` pin** in pyproject is intentional — some scanpy-ecosystem deps still have numpy 2.x incompat issues.

## Test data

`test_data/pbmc3k_filtered_gene_bc_matrices.tar.gz` — canonical 10x PBMC 3k (2700 cells × 32738 genes, hg19, older format with `genes.tsv` not `features.tsv.gz`; scanpy handles both). Standard Seurat-tutorial QC thresholds: `min_genes=200, max_genes=2500, max_pct_mt=5` → ~2638 cells kept.

The synthetic AnnData fixture in `tests/test_steps_smoke.py::make_synthetic_adata` (500 cells × 1500 genes Poisson) runs the full pipeline in well under 30s and is what the smoke test exercises.
