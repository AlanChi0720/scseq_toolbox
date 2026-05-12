import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns


def violin_qc(adata: ad.AnnData):
    fig, axes = plt.subplots(1, 3, figsize=(10, 4))
    keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    for ax, key in zip(axes, keys):
        sns.violinplot(y=adata.obs[key], ax=ax, inner="quartile")
        ax.set_title(key)
        ax.set_xlabel("")
    fig.tight_layout()
    return fig


def elbow_pca(adata: ad.AnnData, n_show: int = 50):
    var_ratio = adata.uns["pca"]["variance_ratio"][:n_show]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, len(var_ratio) + 1), var_ratio, marker="o", markersize=3)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance ratio")
    ax.set_title("PCA elbow")
    fig.tight_layout()
    return fig


def umap_preview(adata: ad.AnnData, color_key: str, title: str = ""):
    fig, ax = plt.subplots(figsize=(4, 4))
    coords = adata.obsm["X_umap"]
    cats = adata.obs[color_key].astype("category")
    palette = sns.color_palette("tab20", n_colors=max(len(cats.cat.categories), 1))
    for i, cat in enumerate(cats.cat.categories):
        mask = (cats == cat).values
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, color=palette[i % len(palette)])
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    return fig
