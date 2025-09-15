import os
import logging
from typing import Union, Sequence

import numpy as np
import pandas as pd
import h5py
import anndata
from anndata.utils import make_index_unique
from scipy.sparse import csr_matrix

# Logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # change to DEBUG/ERROR


# Sampling
def sample_multinomial(x: np.ndarray, n: int, seed: int = 0) -> np.ndarray:
    """Subsample counts using multinomial sampling. Will sample exactly `n` counts.

    Args:
        x (np.ndarray): 1D array of counts per feature/barcode.
        n (int): Number of counts to sample.
        seed (int, optional): Random seed for reproducibility. Defaults to 0.

    Returns:
        np.ndarray: Subsampled counts of the same shape as `x`.

    Raises:
        AssertionError: If `n` >= `x.sum()`.
    """
    total = int(x.sum())
    assert n < total, f"n ({n}) must be smaller than total counts ({total})"
    rng = np.random.default_rng(seed)
    probs = x / total
    return rng.multinomial(n, probs).astype(np.int32)


# MoleculeInfo class
class MoleculeInfo:
    """Represent, manipulate, filter, and summarise molecule info."""

    def __init__(self, path: str, v: int = 1) -> None:
        """
        Args:
            path (str): Path to HDF5 file containing molecule info.
            v (int, optional): Verbosity level. If >0, summary is logged. Defaults to 1.
        """
        self.path = path
        self.v = v
        self.is_sampled = False
        self.n: int = 0
        self._load()
        if v > 0:
            logger.info(str(self))

    def __len__(self) -> int:
        return int((self.count > 0).sum())

    def __str__(self) -> str:
        return f"MoleculeInfo: {len(self)} UMIs | {len(self.feature_names)} features"

    __repr__ = __str__

    # Internal helpers
    def _load(self) -> None:
        """Load molecule info lazily from an HDF5 file (keeps datasets on-demand)."""
        self._h5 = h5py.File(self.path, "r")

        # keep datasets as references (lazy loading)
        self.umi = self._h5["umi"]
        self.count = self._h5["count"][:]  # small enough to load eagerly
        self.count_raw = self.count.copy()
        self.barcode_idx = self._h5["barcode_idx"]
        self.feature_idx = self._h5["feature_idx"]
        self.barcodes = self._h5["barcodes"][:]

        feature_names = self._h5["features"]["name"][:].astype(str)
        feature_names = np.char.upper(feature_names)
        self.feature_names = make_index_unique(pd.Index(feature_names)).values

    def _save(self, path: str) -> None:
        """Save current molecule info to HDF5.

        Args:
            path (str): Output path.
        """
        if os.path.exists(path):
            logger.warning(f"Overwriting existing file at {path}")

        with h5py.File(path, "w") as h5:
            filter_ids = self.count > 0
            h5.create_dataset("umi", data=self.umi[filter_ids])
            h5.create_dataset("count", data=self.count[filter_ids])
            h5.create_dataset("barcode_idx", data=self.barcode_idx[filter_ids])
            h5.create_dataset("feature_idx", data=self.feature_idx[filter_ids])
            h5.create_dataset("barcodes", data=self.barcodes)

            features = h5.create_group("features")
            features.create_dataset("name", data=self.feature_names)

    def _get_counts(self) -> csr_matrix:
        """Return a cell x gene matrix with unique UMI counts.

        Returns:
            csr_matrix: Sparse counts matrix of shape (cells, features).
        """
        filter_ids = self.count > 0
        barcode_idx = self.barcode_idx[filter_ids]
        feature_idx = self.feature_idx[filter_ids]

        combined = np.column_stack(
            (barcode_idx.astype("uint64"), feature_idx.astype("uint64"))
        )
        coords, counts = np.unique(combined, axis=0, return_counts=True)
        return csr_matrix((counts, (coords[:, 0], coords[:, 1])))

    def _subset_data(self, keep_idx: np.ndarray) -> None:
        """Subset and reindex data.

        Args:
            keep_idx (np.ndarray): Boolean mask of UMIs to keep.
        """
        # subset arrays
        self.count_raw = self.count_raw[keep_idx]
        self.count = self.count[keep_idx]
        self.barcode_idx = self.barcode_idx[keep_idx]
        self.feature_idx = self.feature_idx[keep_idx]
        self.umi = self.umi[keep_idx]

        # get unique indices
        unique_barcodes = np.unique(self.barcode_idx)
        unique_features = np.unique(self.feature_idx)

        # subset names
        self.feature_names = self.feature_names[unique_features]
        self.barcodes = self.barcodes[unique_barcodes]

        # update indices efficiently
        self.feature_idx = np.searchsorted(unique_features, self.feature_idx)
        self.barcode_idx = np.searchsorted(unique_barcodes, self.barcode_idx)

        if self.is_sampled:
            logger.info("Data was subset, resetting counts to original.")
            self.is_sampled = False

        if self.v > 0:
            frac = round(sum(keep_idx) / len(keep_idx) * 100, 1)
            logger.info(
                f"Subsetting data to {round(sum(keep_idx) / 1e6)}M ({frac}%) UMIs."
            )

    def _to_array(
        self, x: Union[str, Sequence[str], np.ndarray], dtype=None
    ) -> np.ndarray:
        """Convert input to a numpy array of the given dtype."""
        if isinstance(x, str):
            x = [x]
        return np.array(x, dtype=dtype) if dtype else np.array(x)

    # Public methods
    def select_features(self, feature_names: Union[str, Sequence[str]]) -> None:
        """Select a subset of features by name.

        Args:
            feature_names (Union[str, Sequence[str]]): Names of features to keep.

        Raises:
            ValueError: If any feature names are missing.
        """
        feature_names = self._to_array(feature_names)
        missing = feature_names[~np.isin(feature_names, self.feature_names)]
        if len(missing) > 0:
            raise ValueError(f"Missing {len(missing)} features: {missing[:5]}")

        keep_features = np.isin(self.feature_names, feature_names)
        keep_idx = np.isin(self.feature_idx, np.where(keep_features)[0])
        self._subset_data(keep_idx)

    def select_barcodes(self, barcodes: Union[str, Sequence[str]]) -> None:
        """Select a subset of barcodes by name.

        Args:
            barcodes (Union[str, Sequence[str]]): Barcodes to keep.

        Raises:
            ValueError: If any barcodes are missing.
        """
        barcodes = self._to_array(barcodes).astype(self.barcodes.dtype)
        missing = barcodes[~np.isin(barcodes, self.barcodes)]
        if len(missing) > 0:
            raise ValueError(f"Missing {len(missing)} barcodes: {missing[:5]}")

        keep_barcodes = np.isin(self.barcodes, barcodes)
        keep_idx = np.isin(self.barcode_idx, np.where(keep_barcodes)[0])
        self._subset_data(keep_idx)

    def get_df(self, add_dash_1: bool = True) -> pd.DataFrame:
        """Summarise barcodes.

        Args:
            add_dash_1 (bool, optional): Whether to append "-1" to barcode names. Defaults to True.

        Returns:
            pd.DataFrame: Summary per barcode with columns `umis`, `counts`, `features`.
        """
        filter_ids = self.count > 0
        df = pd.DataFrame(
            {
                "counts": self.count[filter_ids],
                "barcode_idx": self.barcode_idx[filter_ids],
                "feature_idx": self.feature_idx[filter_ids],
            }
        )

        summary = df.groupby("barcode_idx").agg(
            umis=("counts", len),
            counts=("counts", "sum"),
            features=("feature_idx", "nunique"),
        )

        barcodes = self.barcodes[summary.index]
        if add_dash_1:
            barcodes = np.char.add(barcodes.astype(str), "-1")
        summary.index = barcodes
        return summary

    def sample_reads(self, n: int, seed: int = 0) -> None:
        """Subsample counts.

        Args:
            n (int): Number of reads to sample.
            seed (int, optional): Random seed. Defaults to 0.
        """
        self.count = sample_multinomial(self.count_raw, n, seed=seed)
        self.is_sampled = True
        self.n = n

    def reset(self) -> None:
        """Reset counts to original."""
        self.count = self.count_raw
        self.is_sampled = False

    def to_adata(self, add_dash_1: bool = True) -> anndata.AnnData:
        """Convert to AnnData object.

        Args:
            add_dash_1 (bool, optional): Whether to append "-1" to obs names. Defaults to True.

        Returns:
            anndata.AnnData: AnnData object with counts matrix, barcodes as obs, features as var.
        """
        counts = self._get_counts()
        barcodes = self.barcodes[: counts.shape[0]]
        feature_names = self.feature_names[: counts.shape[1]]

        adata = anndata.AnnData(
            X=counts,
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=feature_names),
        )
        if add_dash_1:
            adata.obs_names = adata.obs_names + "-1"
        return adata

    def save(self, path: str) -> None:
        """Save molecule info to HDF5.

        Args:
            path (str): Output file path.
        """
        self._save(path)
        msg = (
            f"Saving {len(self)} UMIs (sampled to {self.n}) to {os.path.basename(path)}."
            if self.is_sampled
            else f"Saving {len(self)} UMIs to {os.path.basename(path)}."
        )
        logger.info(msg)
