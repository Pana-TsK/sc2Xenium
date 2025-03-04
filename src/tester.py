import pandas as pd
import os
import tempfile
import gzip
import scipy
import scipy.sparse as sp
import zipfile
import scanpy as sc
import logging
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class TenXConverter:
    def __init__(self, input_path):
        self.adata = sc.read(input_path)
        self.raw_ad = self._get_raw_counts()

    def _get_raw_counts(self):
        if self.adata.raw is not None:
            return self.adata.raw.to_adata()
        else:
            logging.warning("No `.raw` attribute found. Using `.X` matrix.")
            return self.adata.copy()

    def _force_integer_matrix(self, X):
        """Convert matrix to integers and validate."""
        if sp.issparse(X):
            X = X.copy()
            X.data = np.round(X.data).astype(np.int32)
            # Check for non-integer values
            if not np.all(X.data == X.data.astype(int)):
                logging.error("Non-integer values detected in sparse matrix!")
        else:
            X = np.round(X).astype(np.int32)
            if not np.all(X == X.astype(int)):
                logging.error("Non-integer values detected in dense matrix!")
        return X

    def _validate_files(self, tmp_dir):
        """Validate the generated files before zipping."""
        # Check matrix.mtx.gz
        matrix_path = os.path.join(tmp_dir, "matrix.mtx.gz")
        with gzip.open(matrix_path, "rb") as f:
            header = f.readline().decode().strip()
            if "integer" not in header:
                logging.error("MTX header does not specify integer format!")

        # Check features.tsv.gz and barcodes.tsv.gz
        features = pd.read_csv(os.path.join(tmp_dir, "features.tsv.gz"), sep="\t", header=None)
        barcodes = pd.read_csv(os.path.join(tmp_dir, "barcodes.tsv.gz"), sep="\t", header=None)
        if len(features) != self.raw_ad.n_vars:
            logging.error("Features.tsv does not match the number of genes!")
        if len(barcodes) != self.raw_ad.n_obs:
            logging.error("Barcodes.tsv does not match the number of cells!")

    def h5ad_to_10x(self, output_path="matrix.zip"):
        """Convert AnnData to 10x format with strict integer checks."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Gene metadata
            genes = self.raw_ad.var[["gene_id", "gene"]].copy()
            genes.columns = ["gene_id", "gene_name"]
            genes["feature_type"] = self.raw_ad.var.get("feature_types", "Gene Expression")

            # Barcodes
            barcodes = pd.DataFrame(self.raw_ad.obs.index, columns=["barcode"])

            # Cell type annotations
            if "cell_type" in self.raw_ad.obs.columns:
                celltypes = self.raw_ad.obs[["cell_type"]].reset_index()
                celltypes.columns = ["barcode", "annotation"]
            else:
                celltypes = pd.DataFrame({"barcode": self.raw_ad.obs.index, "annotation": "Unknown"})

            # Write celltypes.csv
            celltypes_path = os.path.join(tmp_dir, "celltypes.csv")
            celltypes.to_csv(celltypes_path, index=False)

            # Matrix
            X = self._force_integer_matrix(self.raw_ad.X.T)
            matrix_path = os.path.join(tmp_dir, "matrix.mtx.gz")
            with gzip.open(matrix_path, "wb") as f:
                scipy.io.mmwrite(f, sp.csc_matrix(X), field='integer')

            # Write features and barcodes
            genes.to_csv(os.path.join(tmp_dir, "features.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")
            barcodes.to_csv(os.path.join(tmp_dir, "barcodes.tsv.gz"), sep="\t", index=False, header=False, compression="gzip")

            # Validate files before zipping
            self._validate_files(tmp_dir)

            # Zip all files
            with zipfile.ZipFile(output_path, "w") as zipf:
                for fname in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "celltypes.csv"]:
                    zipf.write(os.path.join(tmp_dir, fname), arcname=fname)

        logging.info(f"10x files saved to {output_path}")

if __name__ == "__main__":
    input_path = r"C:\Users\panag\OneDrive\Documents\coding\Projects\sc2Xenium\data\adata_file\adata.h5ad"
    converter = TenXConverter(input_path)
    converter.h5ad_to_10x(output_path="matrix.zip")