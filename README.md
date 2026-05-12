# sc_toolbox

Interactive Python scRNA-seq pipeline runner — upload your data, configure a pipeline, and the tool runs automatically, pausing only when you need to make a decision (QC thresholds, n_pcs, clustering resolution). Conceptually a Python version of [scDrake](https://bioinfocz.github.io/scdrake/).

See [DESIGN.md](DESIGN.md) for the full architecture and decisions log.

## Quickstart

Requires Python 3.11. From the project root:

### Windows (PowerShell)
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
streamlit run app.py
```

### macOS / Linux
```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
streamlit run app.py
```

First install pulls scanpy, streamlit, leidenalg, scrublet, … — a few minutes on a fresh machine. After that, `streamlit run app.py` opens a browser tab at `http://localhost:8501`.

## Flow

1. **Configure pipeline** — pick a preset (`standard_10x`), toggle optional steps (doublet detection, cell cycle scoring), tweak default values under "Advanced settings".
2. **Upload data** — `.h5ad`, 10x `.h5`, or 10x mtx folder packaged as `.tar.gz` / `.zip`. Pick the species so the tool knows the mito-gene prefix.
3. **Run** — auto-runs through Load → QC → Normalize → HVG → PCA → Neighbors → Cluster → UMAP, pausing at three decision points:
   - **QC**: violin plots + sliders for `min_genes`, `max_genes`, `min_counts`, `max_counts`, `max_pct_mt`
   - **PCA**: variance-ratio elbow plot + slider for `n_pcs`
   - **Cluster**: four UMAP previews at different resolutions + slider + Leiden/Louvain radio
4. **Download** — `result.h5ad` (full AnnData) and `params.json` (reproducibility manifest).

## Layout

```
sc_toolbox/
├── app.py                  Streamlit entry point
├── pyproject.toml
├── DESIGN.md
├── sc_toolbox/
│   ├── pipeline.py         Step + Pipeline dataclasses
│   ├── presets.py          standard_10x() factory
│   ├── plots.py            mpl figure helpers
│   ├── io_utils.py         file readers (h5ad / h5 / mtx archives)
│   └── steps/              one module per pipeline step
└── tests/
    └── test_steps_smoke.py end-to-end smoke test on synthetic data
```

## Run tests

```bash
pytest -v tests/
```

## Status

v1 — single-sample only. Multi-sample integration, DE, and cell-type annotation are planned for v2.
