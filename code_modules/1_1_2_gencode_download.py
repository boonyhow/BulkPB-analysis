import os
import gzip
import requests
import pandas as pd
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,  # You can set to DEBUG for more messages
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

logger = logging.getLogger(__name__)

def download_gtf(gtf_url, output_path):
    print(f"Downloading GTF from {gtf_url}...")
    response = requests.get(gtf_url, stream=True)
    response.raise_for_status()
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    logger.info("Download complete.")

def parse_gtf(gtf_file):
    logger.info("Parsing GTF for gene_id and gene_name...")
    gene_data = []

    opener = gzip.open if str(gtf_file).endswith(".gz") else open
    with opener(gtf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            attributes = fields[8]
            gene_id = extract_attribute(attributes, "gene_id")
            gene_name = extract_attribute(attributes, "gene_name")
            if gene_id and gene_name:
                gene_data.append((gene_id, gene_name))

    df = pd.DataFrame(gene_data, columns=["gene_id", "gene_name"]).drop_duplicates()
    df['gene_id'] = df['gene_id'].apply(lambda x: x.split(".")[0])
    logger.info(f"Parsed {len(df)} gene entries.")
    return df


def extract_attribute(attr_str, key):
    # e.g., key = "gene_id" or "gene_name"
    for item in attr_str.split(";"):
        item = item.strip()
        if item.startswith(key):
            return item.split('"')[1]
    return None

def main():
    
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir / 'data' / 'gencode_data'
    
    data_dir.mkdir(parents=True, exist_ok=True)
    release = "47"  # or latest GENCODE release
    species = "human"
    gtf_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{release}/gencode.v{release}.basic.annotation.gtf.gz"
    gtf_path = data_dir / f"gencode.v{release}.basic.annotation.gtf.gz"
    csv_output = data_dir / f"gene_id_mapping.csv"

    if not os.path.exists(gtf_path):
        download_gtf(gtf_url, gtf_path)

    gene_df = parse_gtf(gtf_path)
    gene_df.to_csv(csv_output, index=False)
    logger.info(f"Saved gene ID-name mapping to: {csv_output}")

if __name__ == "__main__":
    main()
