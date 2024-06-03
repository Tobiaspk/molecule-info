import os
import numpy as np
import pandas as pd
import h5py
import anndata
import anndata.utils
from scipy.sparse import csr_matrix

def sample_reads(x, n, seed=0):
    """
    Subsample counts.
    
    Example:
    # data looks like
    x = [1, 2, 3, 1]
    total_counts = 7
    
    # draw 4 counts
    sample_indices = [0, 1, 0, 1, 1, 0, 1]  # length 7, keep 1s, drop 0s
    
    # get indices
    count_indices = [0, 1, 1, 2, 2, 2, 3]  # length 7, used for grouping
    
    # group by indices and sum
    sample_counts = [0, 1, 2, 1]
    """
    assert n < x.sum(), f"n ({n}) must be smaller than total counts ({x.sum()})"
    np.random.seed(seed)
    total_count = x.sum()
    sample_indices = np.random.binomial(1, n / total_count, total_count)
    count_indices = np.repeat(np.arange(len(x)), x)
    sample_counts = np.bincount(count_indices, weights=sample_indices).astype("int32")
    return sample_counts

class MoleculeInfo:
    """
    Class to represent, manipulate, filter and summarise molecule info.
    """
    
    def __init__(self, path, v=1):
        self.path = path
        self.v = v
        self.is_sampled = False
        self.n = 0
        self._load()
        if v > 0:
            print(self.__str__())
    
    def __len__(self):
        return len(self.umi[self.count > 0])

    def __str__(self):
        return f"MoleculeInfo: {len(self)} UMIs | {len(self.feature_names)} features"
    
    def __repr__(self):
        return f"MoleculeInfo: {len(self)} UMIs | {len(self.feature_names)} features"
    
    def _load(self):
        with h5py.File(self.path, "r") as h5:
            self.umi = h5["umi"][()]
            self.count = h5["count"][()]
            self.count_raw = self.count.copy()  # count can be overwritten by sample()
            self.barcode_idx = h5["barcode_idx"][()]
            self.feature_idx = h5["feature_idx"][()]
            self.barcodes = h5["barcodes"][()]
            self.feature_names = h5["features"]["name"][()]
            self.feature_names = np.char.upper(self.feature_names.astype(str))
            self.feature_names = anndata.utils.make_index_unique(pd.Index(self.feature_names)).values

    def _save(self, path):
        with h5py.File(path, "w") as h5:
            filter_ids = self.count > 0
            h5.create_dataset("umi", data=self.umi[filter_ids])
            h5.create_dataset("count", data=self.count[filter_ids])
            h5.create_dataset("barcode_idx", data=self.barcode_idx[filter_ids])
            h5.create_dataset("feature_idx", data=self.feature_idx[filter_ids])
            h5.create_dataset("barcodes", data=self.barcodes)
            h5.create_group("features")
            h5["features"].create_dataset("name", data=self.feature_names)

    def _get_counts(self):
        """Return a cell x gene matrix with unique UMI counts for each cell and features"""
        # filter data
        filter_ids = self.count > 0
        barcode_idx = self.barcode_idx[filter_ids]
        feature_idx = self.feature_idx[filter_ids]
        
        # get unique counts
        combined = np.column_stack((barcode_idx.astype("uint64"), feature_idx.astype("uint64")))
        coordinates, counts = np.unique(combined, axis=0, return_counts=True)
        csr_counts = csr_matrix((counts, (coordinates[:, 0], coordinates[:, 1])))
        return csr_counts
    
    def _subset_data(self, keep_idx):
        """Subset and reindex data"""
        
        # subset
        self.count_raw = self.count_raw[keep_idx]
        self.barcode_idx = self.barcode_idx[keep_idx]
        self.feature_idx = self.feature_idx[keep_idx]
        self.umi = self.umi[keep_idx]
        self.count = self.count[keep_idx]
        
        # get unique indices
        unique_barcode_idx = np.unique(self.barcode_idx)
        unique_feature_idx = np.unique(self.feature_idx)
        
        # subset names
        self.feature_names = self.feature_names[unique_feature_idx]
        self.barcodes = self.barcodes[unique_barcode_idx]
        
        # update indices
        map_feature_idx = {old: new for new, old in enumerate(unique_feature_idx)}
        map_barcode_idx = {old: new for new, old in enumerate(unique_barcode_idx)}
        self.feature_idx = np.array([map_feature_idx[i] for i in self.feature_idx])
        self.barcode_idx = np.array([map_barcode_idx[i] for i in self.barcode_idx])

        # outputs
        if self.is_sampled:
            print(f"Data was subset, resetting counts to original.")
            self.is_sampled = False
        
        if self.v > 0:
            print(f"Subsetting data to {round(sum(keep_idx) / 1000000)}M ({round(sum(keep_idx) / len(keep_idx) * 100, 1)}%) UMIs.")
    
    def _to_int(self, x, dtype=None):
        """Convert to array of type"""
        
        # to array
        if isinstance(x, list):
            x = np.array(x, dtype=dtype)
        if isinstance(x, str):
            x = np.array([x], dtype=dtype)
        
        # cast
        if dtype is None:
            return x
        return x.astype(dtype)
                
    def select_features(self, feature_names):
        """Select features by name"""
                
        # init and assert
        feature_names = self._to_int(feature_names)
        missing = feature_names[~np.isin(feature_names, self.feature_names)]
        if len(missing) > 0:
            raise ValueError(f"Missing {len(missing)} features: {missing[:5]}")
        
        # get indices and filter
        keep_features = np.isin(self.feature_names, feature_names)
        keep_idx = np.isin(self.feature_idx, np.where(keep_features)[0])
        self._subset_data(keep_idx)
    
    def select_barcodes(self, barcodes):
        """Select barcodes by name"""
        
        # init and assert
        barcodes = self._to_int(barcodes)
        missing = barcodes[~np.isin(barcodes, self.barcodes)]
        if len(missing) > 0:
            raise ValueError(f"Missing {len(missing)} barcodes: {missing[:5]}")
        
        # try casting to same type
        barcodes = barcodes.astype(self.barcodes.dtype)
        
        # get indices and filter
        keep_barcodes = np.isin(self.barcodes, barcodes)
        keep_idx = np.isin(self.barcode_idx, np.where(keep_barcodes)[0])
        self._subset_data(keep_idx)
    
    def get_df(self, by_cell=True):
        df = pd.DataFrame({
            "counts": self.count,
            "barcode_idx": self.barcode_idx,
            "feature_idx": self.feature_idx,
        })
        if by_cell:
            df = df.groupby("barcode_idx").agg(
                umis=("counts", len),
                counts=("counts", "sum"),
                features=("feature_idx", "nunique")
            )

    def sample_reads(self, n, seed=0):
        """Subsample counts"""
        self.count = sample_reads(self.count_raw, n, seed=seed)
        self.is_sampled = True
        self.n = n
    
    def reset(self):
        """Reset counts to original"""
        self.count = self.count_raw
        self.is_sampled = False
    
    def get_df(self, by_cell=True, add_dash_1=True):
        """Summarise barcodes"""
        
        # TODO: filter features?
        
        if not by_cell:
            raise NotImplementedError("by_cell=False not implemented yet")
        
        # filter
        filter_ids = self.count > 0
        count = self.count[filter_ids]
        barcode_idx = self.barcode_idx[filter_ids]
        feature_idx = self.feature_idx[filter_ids]
        barcodes = self.barcodes[barcode_idx]
        if self.v > 0:
            print(f"Aggregating {len(barcodes)} barcodes...")

        # create df
        df = pd.DataFrame({
            "counts": count,
            "barcode_idx": barcode_idx,
            "feature_idx": feature_idx,
        })
        
        # aggregate by cell
        df = df.groupby("barcode_idx").agg(
            umis=("counts", len),
            counts=("counts", "sum"),
            features=("feature_idx", "nunique")
        )
        
        # set index to barcode names
        index_new = barcodes[df.index]
        if add_dash_1:
            index_new = np.core.defchararray.add(barcodes, "-1")
        df.index = index_new
        return df
    
    def to_adata(self, add_dash_1=True):
        """Convert to adata. Calculates unique counts."""
        
        # get counts
        counts = self._get_counts()
        
        # filter barcodes (index cannot be higher than counts dimension)
        barcodes = self.barcodes[:counts.shape[0]]
        feature_names = self.feature_names[:counts.shape[1]]
        
        # convert to adata
        adata = anndata.AnnData(
            X=counts,
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=feature_names)
        )
        if add_dash_1:
            adata.obs_names = adata.obs_names + "-1"
        return adata
            
    # TODO: Saving function
    def save(self, path):
        self._save(path)
        
        if not self.is_sampled:
            msg = f"Saving {len(self)} UMIs to {os.path.basename(path)}."
        else:
            msg = f"Saving {len(self)} UMIs (sampled to {self.n}) to {os.path.basename(path)}."
        print(msg)