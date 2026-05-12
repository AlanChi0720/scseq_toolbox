from sc_toolbox.pipeline import Pipeline
from sc_toolbox.steps import (
    cluster,
    doublet,
    hvg,
    load,
    neighbors,
    normalize,
    pca,
    qc,
    umap,
)


def standard_10x() -> Pipeline:
    return Pipeline(
        name="standard_10x",
        steps=[
            load.step,
            qc.step,
            doublet.step,
            normalize.step,
            hvg.step,
            pca.step,
            neighbors.step,
            cluster.step,
            umap.step,
        ],
        enabled={
            "load": True,
            "qc": True,
            "doublet": False,
            "normalize": True,
            "hvg": True,
            "pca": True,
            "neighbors": True,
            "cluster": True,
            "umap": True,
        },
    )


PRESETS = {"standard_10x": standard_10x}
