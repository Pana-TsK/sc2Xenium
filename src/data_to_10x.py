import pandas as pd
import os
import tempfile
import gzip
import scipy
import scipy.sparse as sp
import zipfile
import scanpy as sc
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
class TenXConverter:
    def __init__(self, input_path,):
        self.adata = self.load_anndata(input_path)
        self.raw_ad = self.retrieve_rawcounts()

    def load_anndata(self, path):
        """
        Load in the Anndata object from the specified path.
        """
        try:
            ad = sc.read(path) # Specify the path to the AnnData object
        except Exception as e:
            raise ValueError(f"Could not read AnnData object: {e}")
        
        return ad
    
    def retrieve_rawcounts(self):
        """
        Attempt to retrieve the raw counts from the Anndata object.
        If they are not retrieved, the function continues with the full dataset.
        """
        if self.adata.raw is not None:
            raw_ad = self.adata.raw.to_adata()
        else:
            print("No raw counts found. Using the X matrix instead, with the assumption this contains non-normalized data.")
            raw_ad = self.adata.copy()
        
        return raw_ad

    def h5ad_to_10x(self, 
                    gene_id_key="gene_id", 
                    gene_name_key="gene", 
                    cell_type_key="cell_type", 
                    output_path="matrix.zip", 
                    barcode_key=None, 
                    subsample_rate=None):

                    
        """
        Convert an AnnData object to 10x Genomics format and package the results into a zip file.
        
        Parameters:
        gene_id_key (str): Column name in ad.var for gene IDs. Standard is "gene_id".
        gene_name_key (str): Column name in ad.var for gene names. Standard is "gene".
        cell_type_key (str): Column name in ad.obs for cell type annotations. Standard is "cell_type".
        output_path (str): Path for the output zip file. Will be overwritten if it already exists.
        barcode_key (str): Optional key in ad.obs for barcodes; if None, uses ad.obs.index.
        subsample_rate (float): Optional fraction (0-1) for subsampling cells. 
        """

        # Inform that this function may take a while
        logging.info("Converting AnnData object to 10x format. This may take a while for large datasets.")

        if subsample_rate:
            sc.pp.subsample(self.raw_ad, fraction=subsample_rate)
        
        logging.info("Extracting gene information")
        
        # Extract gene information
        genes = self.raw_ad.var.reset_index()[[gene_id_key, gene_name_key]].copy()
        genes.columns = ["gene_id", "gene_name"]  # Rename columns explicitly

        logging.info("Adding feature types")
        
        # Add feature type
        if "feature_types" in self.raw_ad.var.columns:
            genes["feature_type"] = self.raw_ad.var["feature_types"].values
        else:
            genes["feature_type"] = "Gene Expression"

        logging.info("Extracting barcode and cell type annotations")
        
        # Extract barcode information
        if barcode_key and barcode_key in self.raw_ad.obs.columns:
            barcodes = self.raw_ad.obs[[barcode_key]].copy()
        else:
            barcodes = pd.DataFrame(self.raw_ad.obs.index, columns=["barcode"])

        logging.info("Creating zip files")
        
        # Extract cell type annotations
        if cell_type_key in self.raw_ad.obs.columns:
            celltypes = self.raw_ad.obs[[cell_type_key]].reset_index()
            celltypes.columns = ["barcode", "annotation"]
        else:
            celltypes = pd.DataFrame({"barcode": self.raw_ad.obs.index, "annotation": "Unknown"})
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Write matrix in compressed MTX format (transposed to genes x cells)
            matrix_file = os.path.join(tmp_dir, "matrix.mtx.gz")
            with gzip.open(matrix_file, "wb") as handle:
                scipy.io.mmwrite(handle, sp.csc_matrix(self.raw_ad.X.T))
            
            # Write features file (genes)
            features_file = os.path.join(tmp_dir, "features.tsv.gz")
            genes.to_csv(features_file, sep="\t", index=False, header=False, compression="gzip")
            
            # Write barcodes file
            barcodes_file = os.path.join(tmp_dir, "barcodes.tsv.gz")
            barcodes.to_csv(barcodes_file, sep="\t", index=False, header=False, compression="gzip")
            
            # Write cell types file
            celltypes_file = os.path.join(tmp_dir, "celltypes.csv")
            celltypes.to_csv(celltypes_file, index=False)
            
            # Validate file contents
            with gzip.open(features_file, "rt") as f:
                lines = [line.strip() for line in f]
                logging.info(f"Features.tsv first line: {lines[0]}")
                logging.info(f"Features.tsv has {len(lines)} rows (should match matrix.mtx rows).")
            
            with gzip.open(matrix_file, "rt") as f:
                header = next(line for line in f if not line.startswith("%"))
                logging.info(f"Matrix header: {header.strip()}")
            
            # Package all files into a zip archive
            with zipfile.ZipFile(output_path, "w") as zip_handle:
                for file_name in ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "celltypes.csv"]:
                    zip_handle.write(os.path.join(tmp_dir, file_name), arcname=file_name)
        
        logging.info(f"10x formatted files zipped successfully to {output_path}.")

if __name__ == "__main__":
    # Example usage
    import scanpy as sc

    input_path = r"C:\Users\panag\OneDrive\Documents\coding\Projects\Liliana\data\output_anndata.h5ad"
    converter = TenXConverter(input_path)

    
    # Convert to 10x format
    converter.h5ad_to_10x( 
        gene_id_key="gene_id", 
        gene_name_key="gene",  # Ensure this is distinct from gene_id_key
        cell_type_key="cell_type", 
        output_path=r"C:\Users\panag\OneDrive\Documents\coding\Projects\Liliana\matrix.zip"
    )