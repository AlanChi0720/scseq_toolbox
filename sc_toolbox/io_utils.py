import tarfile
import tempfile
import zipfile
from pathlib import Path

import anndata as ad
import scanpy as sc


def load_uploaded(uploaded_file) -> ad.AnnData:
    """Read a Streamlit UploadedFile into AnnData.

    Supports: .h5ad, .h5 (10x), .tar.gz / .tgz, .zip (10x mtx packaged).
    """
    name = uploaded_file.name.lower()
    data = bytes(uploaded_file.getbuffer())
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        raw = tmp_path / uploaded_file.name
        raw.write_bytes(data)

        if name.endswith(".h5ad"):
            return sc.read_h5ad(raw)
        if name.endswith(".h5"):
            return sc.read_10x_h5(raw)
        if name.endswith(".tar.gz") or name.endswith(".tgz"):
            extract_to = tmp_path / "extracted"
            extract_to.mkdir()
            with tarfile.open(raw, "r:gz") as t:
                t.extractall(extract_to)
            return _read_10x_mtx_dir(extract_to)
        if name.endswith(".zip"):
            extract_to = tmp_path / "extracted"
            extract_to.mkdir()
            with zipfile.ZipFile(raw, "r") as z:
                z.extractall(extract_to)
            return _read_10x_mtx_dir(extract_to)
        raise ValueError(f"Unsupported file format: {uploaded_file.name}")


def _read_10x_mtx_dir(root: Path) -> ad.AnnData:
    candidates = [root, *(p for p in root.rglob("*") if p.is_dir())]
    for cand in candidates:
        contents = {p.name for p in cand.iterdir()}
        if any(n.startswith("matrix.mtx") for n in contents):
            return sc.read_10x_mtx(cand, var_names="gene_symbols", make_unique=True)
    raise ValueError("Could not locate matrix.mtx[.gz] inside the upload")
