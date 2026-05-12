import scanpy as sc

from sc_toolbox.pipeline import Step


def _run(adata, params):
    sc.pp.scrublet(adata, verbose=False)
    if "predicted_doublet" in adata.obs.columns:
        adata = adata[~adata.obs.predicted_doublet.fillna(False)].copy()
    return adata


step = Step(
    name="doublet", label="Doublet detection (Scrublet)", auto=True, run=_run
)
