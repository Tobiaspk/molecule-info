import os
import anndata as ad
from molecule_info import utils, MoleculeInfo
import pytest


def unique(x):
    return list(set(x))


@pytest.fixture
def molecule_info():
    path = os.path.join(
        os.path.dirname(utils.get_package_location()),
        "tests",
        "data",
        "molecule_info.h5",
    )
    return MoleculeInfo(path)


def test_operations(molecule_info):
    # Number UMIs
    n_umis = 1012
    assert len(molecule_info) == n_umis & len(molecule_info.umi) == n_umis, (
        f"Expected {n_umis} UMIs, got {len(molecule_info)} UMIs"
    )

    # Number Features
    n_features = 6
    unique_features = unique(molecule_info.feature_names)
    assert len(unique_features) == n_features, (
        f"Expected {n_features} features, got {len(unique_features)} features"
    )

    # Number Barcodes
    n_barcodes = 100
    unique_barcodes = unique(molecule_info.barcodes)
    assert len(unique_barcodes) == n_barcodes, (
        f"Expected {n_barcodes} barcodes, got {len(unique_barcodes)} barcodes"
    )

    # Test Subset Features
    keep = ["C0315", "C0316"]
    molecule_info.select_features(keep)
    unique_features = unique(molecule_info.feature_names)
    assert len(unique_features) == len(keep), (
        f"Expected {len(keep)} features, got {len(unique_features)} features"
    )

    # Test Subset Barcodes
    keep = molecule_info.barcodes[:10]
    molecule_info.select_barcodes(keep)
    unique_barcodes = unique(molecule_info.barcodes)
    assert len(unique_barcodes) == len(keep), (
        f"Expected {len(keep)} barcodes, got {len(unique_barcodes)} barcodes"
    )


def test_sampling(molecule_info):
    # Sample counts
    n_counts = 100
    molecule_info.sample_reads(n_counts, seed=42)
    assert molecule_info.is_sampled, "MoleculeInfo should be marked as sampled"
    assert len(molecule_info) <= 100, (
        f"Expected {n_counts} counts after sampling, got {len(molecule_info)} counts"
    )
    assert sum(molecule_info.count) == n_counts, (
        f"Expected sum of counts to be {n_counts}, got {sum(molecule_info.count)}"
    )

    adata = molecule_info.to_adata()
    assert adata.X.sum() <= n_counts, (
        f"Expected sum of counts in AnnData to be ~ {n_counts}, got {adata.X.sum()}"
    )
    assert adata.shape[0] <= n_counts, (
        f"Expected number of cells in AnnData to be <= {n_counts}, got {adata.shape[0]}"
    )

    # Reset
    molecule_info.reset()
    assert not molecule_info.is_sampled, "MoleculeInfo should not be marked as sampled"
    assert len(molecule_info) == 1012, (
        f"Expected {1012} UMIs after reset, got {len(molecule_info)} UMIs"
    )


def test_adata(molecule_info):
    # To anndata
    adata = molecule_info.to_adata()
    assert isinstance(adata, ad.AnnData), f"Expected AnnData object, got {type(adata)}"
    assert adata.shape[0] == len(molecule_info.barcodes), (
        f"Expected {len(molecule_info.barcodes)} cells, got {adata.shape[0]} cells"
    )
    assert adata.shape[1] == len(molecule_info.feature_names), (
        f"Expected {len(molecule_info.feature_names)} features, got {adata.shape[1]} features"
    )
