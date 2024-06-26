# Molecule Info H5 Helpers

This repository contains a set of helper functions for reading, subsampling and in future versions visualizing the data in the molecule info h5 files as generated by cellranger.

## Installation

```bash
pip install molecule_info
```

## Usage

### Subsampling

```python
import molecule_info as mi

# read
m = mi.MoleculeInfo(path)

# subsample
m.sample_reads(1000, seed=42)

# convert to anndata
adata = m.to_adata()
```

### Select subset of features

```python
import molecule_info as mi

# read
m = mi.MoleculeInfo(path)

# select features (inplace!)
features = ["C0310", "C0311", "C0312"]
m.select_features(features)

# whitelist
barcodes = ["AAACCTGAGGAGTCTG-1", "AAACCTGAGGAGTCTC-1"]
m.select_barcodes(barcodes)

# convert to anndata
adata = m.to_adata()
```