import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import h5py
import logging
import scipy.io
import gzip
import tempfile
import zipfile

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class AnndataConverter:
    def __init__(self, cells_input, genes_input, matrix_input):
        self.cells_input = cells_input
        self.genes_input = genes_input
        self.matrix_input = matrix_input
        self.h5_path = None
        self.mtx_path = None
        self.validate_input_files()
    
    def validate_input_files(self):
        """Ensure input files exist and are not empty."""
        for path in [self.cells_input, self.genes_input, self.matrix_input]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"File not found: {path}")
            if os.path.getsize(path) == 0:
                raise ValueError(f"File is empty: {path}")
        logging.info("Input files validated successfully.")
    
    def read_csv(self, path):
        """Read a CSV file into a DataFrame."""
        try:
            df = pd.read_csv(path, index_col=None)
            logging.info(f"File read from {path}. Shape: {df.shape}")
            return df
        except Exception as e:
            logging.error(f"Error reading file {path}: {e}")
            raise
    
    def convert_to_mtx(self, output_dir=None):
        """
        Convert the CSV matrix file to Matrix Market (MTX) format.
        Generates three files: matrix.mtx, barcodes.tsv, and features.tsv.
        """
        if output_dir is None:
            output_dir = os.path.join(os.path.dirname(self.matrix_input), "mtx_files")
        os.makedirs(output_dir, exist_ok=True)
        
        matrix_path = os.path.join(output_dir, "matrix.mtx")
        barcodes_path = os.path.join(output_dir, "barcodes.tsv")
        features_path = os.path.join(output_dir, "features.tsv")
        
        if os.path.exists(matrix_path) and os.path.exists(barcodes_path) and os.path.exists(features_path):
            logging.info(f"MTX files already exist in {output_dir}.")
            self.mtx_path = output_dir
            return output_dir
        
        logging.info("Converting CSV matrix to MTX format...")
        
        # Read gene metadata from genes_input
        genes_df = self.read_csv(self.genes_input)
        # Ensure the first column is named "gene_id"
        if genes_df.columns[0] != "gene_id":
            logging.info(f"Renaming first column of genes file to 'gene_id'.")
            new_cols = ["gene_id"] + list(genes_df.columns[1:])
            genes_df.columns = new_cols
        # Save gene information (features.tsv) without header/index
        genes_df.to_csv(features_path, sep='\t', index=False, header=False)
        
        # Process the matrix CSV file (first column: cell barcode; remaining columns: expression values)
        cell_barcodes = []
        data = []
        row_indices = []  # Gene indices (1-indexed for MTX)
        col_indices = []  # Cell indices (1-indexed for MTX)
        
        chunk_size = 1000
        row_count = 0
        
        for chunk in pd.read_csv(self.matrix_input, chunksize=chunk_size):
            # The first column contains cell barcodes.
            cell_barcodes.extend(chunk.iloc[:, 0].tolist())
            for _, row in chunk.iterrows():
                # Expression values for genes
                expression_values = row.iloc[1:].values
                non_zero_mask = expression_values > 0
                if np.any(non_zero_mask):
                    non_zero_values = expression_values[non_zero_mask]
                    data.extend(non_zero_values)
                    gene_indices = np.where(non_zero_mask)[0] + 1  # Convert to 1-indexed
                    row_indices.extend(gene_indices)
                    col_indices.extend([row_count + 1] * np.sum(non_zero_mask))
                row_count += 1
                if row_count % 1000 == 0:
                    logging.info(f"Processed {row_count} cells.")
        
        # Write matrix.mtx in Matrix Market format
        with open(matrix_path, 'w') as f:
            f.write("%%MatrixMarket matrix coordinate real general\n")
            # Header: rows (genes) x cols (cells) x nonzeros
            f.write(f"{genes_df.shape[0]} {row_count} {len(data)}\n")
            for i in range(len(data)):
                f.write(f"{row_indices[i]} {col_indices[i]} {data[i]}\n")
        
        # Write barcodes.tsv
        with open(barcodes_path, 'w') as f:
            for barcode in cell_barcodes:
                f.write(f"{barcode}\n")
        
        logging.info(f"MTX conversion complete: {row_count} cells, {genes_df.shape[0]} genes, {len(data)} nonzero entries.")
        self.mtx_path = output_dir
        return output_dir
    
    def create_anndata_from_mtx(self):
        """
        Read the generated MTX files (matrix, barcodes, features) and create an AnnData object.
        This method adds cell and gene metadata.
        """
        if self.mtx_path is None:
            self.convert_to_mtx()
        
        logging.info(f"Creating AnnData object from MTX files in {self.mtx_path}...")
        
        # Read matrix.mtx and transpose
        adata = sc.read_mtx(os.path.join(self.mtx_path, "matrix.mtx")).T
        
        # Read barcodes (cell identifiers)
        with open(os.path.join(self.mtx_path, "barcodes.tsv"), 'r') as f:
            barcodes = [line.strip() for line in f]
        adata.obs_names = pd.Index(barcodes)
        
        # Read features (gene metadata)
        # Get column names from genes_input
        genes_df = self.read_csv(self.genes_input)
        if genes_df.columns[0] != "gene_id":
            new_cols = ["gene_id"] + list(genes_df.columns[1:])
            genes_df.columns = new_cols
        feature_columns = genes_df.columns.tolist()
        
        features = pd.read_csv(
            os.path.join(self.mtx_path, "features.tsv"), 
            sep='\t', 
            header=None, 
            names=feature_columns
        )
        
        # Set var_names to gene_id and assign metadata to var
        adata.var_names = features['gene_id'].astype(str)
        if adata.var_names.duplicated().any():
            logging.warning("Duplicate gene names found in features.tsv. Making them unique.")
            adata.var_names_make_unique()
        
        # Assign remaining columns to var
        adata.var = features.copy() # for the love of God please make this assign all columns 
        
        # Add cell metadata from cells_input
        cells_df = self.read_csv(self.cells_input)
        if cells_df.columns[0] != "barcode":
            logging.info(f"Renaming first column of cells file to 'barcode'.")
            new_cols = ["barcode"] + list(cells_df.columns[1:])
            cells_df.columns = new_cols
        cells_df.set_index("barcode", inplace=True)
        
        # Ensure cells in adata.obs exist in cells_df
        common_cells = adata.obs_names.intersection(cells_df.index)
        if len(common_cells) < len(adata.obs_names):
            logging.warning("Some cells in matrix are not present in cells metadata. Filtering...")
            adata = adata[common_cells, :].copy()
        
        adata.obs = cells_df.loc[adata.obs_names]
        
        logging.info(f"Created AnnData object with shape {adata.shape}.")
        return adata


def main():
    # Specify your file paths here
    cells_input = r"C:\Users\panag\OneDrive\Documents\coding\Projects\Liliana\data\SC_GEO_cells.csv"
    genes_input = r"C:\Users\panag\OneDrive\Documents\coding\Projects\Liliana\data\SC_GEO_genes.csv"
    matrix_input = r"C:\Users\panag\OneDrive\Documents\coding\Projects\Liliana\data\GSE192456_SC_GEO_raw_counts.csv.gz"
    
    # Create an instance of the converter
    converter = AnndataConverter(cells_input, genes_input, matrix_input)
    
    # Option 1: Convert to MTX format and create an AnnData object
    converter.convert_to_mtx()
    adata = converter.create_anndata_from_mtx()
    
    # Save the AnnData object to disk
    output_h5ad = os.path.join(os.path.dirname(matrix_input), "output_anndata.h5ad")
    adata.write_h5ad(output_h5ad)
    logging.info(f"AnnData object saved to {output_h5ad}.")
    
    return adata

if __name__ == "__main__":
    adata = main()

    print(adata)