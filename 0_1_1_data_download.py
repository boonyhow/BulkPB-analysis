import os
import tarfile
import urllib.request
from pathlib import Path
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,  # You can set to DEBUG for more messages
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

logger = logging.getLogger(__name__)

def download_and_untar(ftp_url, filename, errors):
    try:
        # Define paths
        script_dir = Path(__file__).parent.resolve()
        data_dir = script_dir / 'data' / filename
        raw_tar_dir = data_dir / 'raw_tars'
        extracted_dir = data_dir / 'extracted'

        # Create directories if they don't exist
        raw_tar_dir.mkdir(parents=True, exist_ok=True)
        extracted_dir.mkdir(parents=True, exist_ok=True)

        # Parse filename from URL
        parsed_url = urlparse(ftp_url)
        base_filename = os.path.basename(parsed_url.path)

        # Fallback if URL has no valid basename
        if not base_filename or base_filename.strip() == '':
            base_filename = f"{filename}.tar.gz"
            logger.warning(f"No valid filename in URL for {filename}, fallback to {base_filename}")

        file_path = raw_tar_dir / base_filename

        # Check if file already downloaded
        if file_path.exists():
            logger.info(f"File already downloaded: {file_path}")
        else:
            logger.info(f"Downloading {ftp_url} -> {file_path}")
            urllib.request.urlretrieve(ftp_url, file_path)
            logger.info(f"Download complete: {file_path}")

        # Check if already extracted
        if any(extracted_dir.iterdir()):
            logger.info(f"Already extracted: {extracted_dir}")
        else:
            # Extract if tar
            if tarfile.is_tarfile(file_path):
                logger.info(f"Extracting {file_path} -> {extracted_dir}")
                with tarfile.open(file_path, 'r:*') as tar:
                    tar.extractall(path=extracted_dir)
                logger.info(f"Extraction complete: {extracted_dir}")
            else:
                logger.warning(f"Downloaded file {file_path} is not a tar archive. Skipping extraction.")

    except Exception as e:
        logger.error(f"Error processing {filename}: {e}")
        errors.append((filename, str(e)))

def batch_download_and_untar(url_filename_list, max_workers=4):
    errors = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for ftp_url, filename in url_filename_list:
            futures.append(executor.submit(download_and_untar, ftp_url, filename, errors))

        for _ in tqdm(as_completed(futures), total=len(futures), desc="Downloading and extracting"):
            pass

    # Summary
    if errors:
        logger.error("\nSome downloads/extractions failed:")
        for filename, error in errors:
            logger.error(f"- {filename}: {error}")
    else:
        logger.info("\nAll files downloaded and extracted successfully.")

if __name__ == "__main__":
    # Example list of (url, filename)
    url_filename_list = [
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt.gz", "GSE176078"),
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz", "GSE176078"),
        ("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217517&format=file", 'GSE217517'),
        ("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226163&format=file", 'GSE226163')
    ]

    batch_download_and_untar(url_filename_list, max_workers=4)
