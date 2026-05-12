"""Cell cycle scoring step (Tirosh et al. 2016 gene signatures).

Scores S-phase and G2/M-phase activity per cell and assigns a discrete phase
label (G1 / S / G2M) to `adata.obs`. Does not regress the signal out — that
remains opt-in and is deferred to a later version.
"""
import scanpy as sc

from sc_toolbox.pipeline import Step

# Tirosh et al. 2016 cell-cycle marker genes (human, uppercase symbols).
S_GENES_HUMAN = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
]
G2M_GENES_HUMAN = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


def _to_mouse(symbols: list[str]) -> list[str]:
    """Convert human gene symbols to mouse convention (Title-case)."""
    return [s.capitalize() for s in symbols]


def _gene_lists_for_prefix(mito_prefix: str) -> tuple[list[str], list[str]]:
    """Pick S / G2M gene lists from the species inferred via mito prefix."""
    if mito_prefix == "mt-":
        return _to_mouse(S_GENES_HUMAN), _to_mouse(G2M_GENES_HUMAN)
    # Human, or unknown custom prefix → default to human.
    return list(S_GENES_HUMAN), list(G2M_GENES_HUMAN)


def _run(adata, params):
    mito_prefix = adata.uns.get("sc_toolbox", {}).get("mito_prefix", "MT-")
    s_genes, g2m_genes = _gene_lists_for_prefix(mito_prefix)

    var_index = set(adata.var_names)
    s_present = [g for g in s_genes if g in var_index]
    g2m_present = [g for g in g2m_genes if g in var_index]

    if not s_present or not g2m_present:
        adata.uns.setdefault("sc_toolbox", {})["cell_cycle_skipped"] = (
            f"no gene overlap (S={len(s_present)}, G2M={len(g2m_present)})"
        )
        return adata

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_present, g2m_genes=g2m_present)
    return adata


step = Step(
    name="cell_cycle",
    label="Cell cycle scoring",
    auto=True,
    run=_run,
)
