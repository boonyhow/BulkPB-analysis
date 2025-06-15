import pandas as pd
from pathlib import Path
from tqdm import tqdm
import pyarrow.csv as pv
from scipy import sparse
from scipy.io import mmwrite
import gzip
import subprocess
import shutil

def gzip_file(filepath):
    gz_path = filepath.with_suffix(filepath.suffix + ".gz")
    with open(filepath, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    filepath.unlink()  # delete uncompressed file
    return gz_path


def gzip_except_metadata(folder_path, use_pigz=True):
    print("Compressing files...")

    compressor = "pigz" if use_pigz else "gzip"
    cmd = [
        "find", str(folder_path),
        "-type", "f", "!", "-name", "metadata.csv",
        "-exec", compressor, "-f", "{}", ";"
    ]

    result = subprocess.run(cmd, check=True)
    print(f"✅ Compressed files (excluding metadata.csv) using {compressor} with overwrite enabled.")

    # Just to confirm cleanup
    print("Verifying uncompressed files were removed...")
    for path in Path(folder_path).rglob("*"):
        if path.is_file() and not str(path).endswith(".gz") and path.name != "metadata.csv":
            print(f"⚠️ Uncompressed file still exists: {path}")

def organize_extracted_data_GSE226163(main_folder):
    main_dir = Path(main_folder)
    input_dir = main_dir / 'extracted'
    output_dir = main_dir / 'data_for_study'
    output_dir.mkdir(parents=True, exist_ok=True)

    bulk_outdir = output_dir / "bulk_rna_seq"
    sc_outdir = output_dir / "single_cell"
    bulk_outdir.mkdir(parents=True, exist_ok=True)
    sc_outdir.mkdir(parents=True, exist_ok=True)

    all_files = list(input_dir.glob("*.gz"))
    bulk_files = [f for f in all_files if "mR" in f.name]
    sc_files = [f for f in all_files if "SC" in f.name]

    # === Bulk RNA-seq ===
    bulk_data = {}
    for f in bulk_files:
        name_part = f.name.split("_")[1]
        sample_name = name_part.replace("mR.pergene.counts.txt.gz", "")
        df = pv.read_csv(f, parse_options=pv.ParseOptions(delimiter="\t")).to_pandas()
        df.columns = ["Gene", sample_name]
        bulk_data[sample_name] = df

    merged_bulk = None
    for sample_name, df in bulk_data.items():
        merged_bulk = df if merged_bulk is None else pd.merge(merged_bulk, df, on="Gene", how="inner")

    merged_bulk.to_csv(bulk_outdir / "combined_bulk.csv", index=False)
    print("✅ Saved: bulk_rna_seq/combined_bulk.csv")

    # === scRNA-seq to 10X-style output ===
    dfs = []
    metadata_rows = []

    for file in tqdm(sorted(sc_files), desc="Processing scRNA-seq files"):
        df = pv.read_csv(file).to_pandas()
        df = df.set_index(df.columns[0])
        name_part = file.name.split("_")[1]
        clean_sample = name_part.replace("SC.count.csv.gz", "")
        df.columns = [f"{clean_sample}_{col.split('.')[0]}" for col in df.columns]
        for col in df.columns:
            metadata_rows.append({"barcode": col, "orig.ident": clean_sample})
        dfs.append(df)

    merged_counts = pd.concat(dfs, axis=1)

    # Save metadata uncompressed
    pd.DataFrame(metadata_rows).to_csv(sc_outdir / "metadata.csv", index=False)
    print("✅ Saved: single_cell/metadata.csv")

    # Create sparse matrix and save
    matrix_sparse = sparse.csc_matrix(merged_counts.values)
    mm_path = sc_outdir / "matrix.mtx"
    mmwrite(mm_path, matrix_sparse)
    gzip_file(mm_path)

    # Save barcodes
    barcodes_path = sc_outdir / "barcodes.tsv"
    pd.DataFrame(merged_counts.columns).to_csv(barcodes_path, header=False, index=False)
    gzip_file(barcodes_path)

    # Save features (3-column)
    features_path = sc_outdir / "features.tsv"
    genes = merged_counts.index
    features = pd.DataFrame({
        0: [f"GENE{i+1:06d}" for i in range(len(genes))],
        1: genes,
        2: "Gene Expression"
    })
    features.to_csv(features_path, sep="\t", header=False, index=False)
    gzip_file(features_path)
    # gzip_except_metadata(sc_outdir)
    print("✅ Saved: matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz in 10X format")
    
    # === Subtype metadata generation ===
    all_sample_names = list(merged_bulk.columns[1:])
    subtype_records = []

    for sample in all_sample_names:
        if sample.startswith("iVC"):
            subtype = "Smooth Muscle cells"
        elif sample.startswith("iEC"):
            subtype = "Endothelial cells"
        else:
            subtype = "Unknown"
        subtype_records.append({"Sample": sample, "Subtype": subtype})

    subtype_df = pd.DataFrame(subtype_records).drop_duplicates("Sample")
    subtype_df.to_csv(output_dir / "sample_subtype.csv", index=False)
    print("✅ Saved: data_for_study/sample_subtype.csv")


if __name__ == "__main__":
    organize_extracted_data_GSE226163("data/GSE226163")
