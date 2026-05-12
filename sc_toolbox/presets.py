from sc_toolbox.pipeline import Pipeline
from sc_toolbox.steps import (
    cell_cycle,
    cluster,
    doublet,
    hvg,
    load,
    markers,
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
            cell_cycle.step,
            hvg.step,
            pca.step,
            neighbors.step,
            cluster.step,
            markers.step,
            umap.step,
        ],
        enabled={
            "load": True,
            "qc": True,
            "doublet": False,
            "normalize": True,
            "cell_cycle": False,
            "hvg": True,
            "pca": True,
            "neighbors": True,
            "cluster": True,
            "markers": True,
            "umap": True,
        },
    )


PRESETS = {"standard_10x": standard_10x}
