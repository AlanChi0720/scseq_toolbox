from sc_toolbox.pipeline import Step


def _run(adata, params):
    return adata


step = Step(name="load", label="Load data", auto=True, run=_run)
