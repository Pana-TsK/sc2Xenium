import zipfile
import scanpy as sc
import tempfile
import pandas as pd

class DataTester:
    """
    A class for testing data integrity.
    """
    def __init__(self, zip_path):
        self.zip_path = zip_path

    def list_zip_contents(self):
        """List the contents of a zip file."""
        try:
            with zipfile.ZipFile(self.zip_path, 'r') as zipf:
                print(f"Contents of {self.zip_path}:")
                for file in zipf.namelist():
                    print(f"  - {file}")
        except zipfile.BadZipFile:
            print(f"Error: {self.zip_path} is not a valid zip file.")
        except FileNotFoundError:
            print(f"Error: File not found - {self.zip_path}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
    

    def matrix_reader(self):
        # Read the matrix file
        with tempfile.TemporaryDirectory() as tmp_dir:
            with zipfile.ZipFile(self.zip_path, "r") as zip_handle:
                zip_handle.extractall(tmp_dir)
            mtx = sc.read_10x_mtx(tmp_dir)
        
        return mtx.to_df().head()

