# sc2Xenium: Single-Cell to Xenium Reference Converter

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![anndata](https://img.shields.io/badge/Built%20with-AnnData-orange)](https://anndata.readthedocs.io/)

A Python toolkit for converting single-cell data from CSV or 10X Genomics format to Xenium-compatible AnnData objects and MTX references.

## Overview

sc2Xenium provides streamlined conversion of single-cell RNA-seq data into standardized formats compatible with 10X Genomics Xenium platforms. The toolkit offers:

- Conversion from raw CSVs or 10X Genomics outputs to AnnData objects
- Automated raw count matrix handling
- Generation of Xenium Custom Panel-ready MTX references
- Data validation and integrity checks

## Key Features

🔧 **Three Core Components:**
1. `AnndataConverter`: Creates AnnData objects from CSV inputs
2. `TenXConverter`: Processes 10X Genomics data to AnnData/Xenium formats
3. `DataTester`: Validates data integrity and format compliance

📊 **Input Support:**
- Raw count matrices (`raw_counts.csv.gz`)
- Cell metadata (`cells.csv`)
- Gene metadata (`genes.csv`)
- 10X Genomics directory structure

🎯 **Outputs:**
- Validated AnnData objects
- Xenium-ready MTX files (features, matrix, barcodes)
- QC reports and conversion metadata

