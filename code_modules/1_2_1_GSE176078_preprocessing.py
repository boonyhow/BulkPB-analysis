import pandas as pd
from pathlib import Path
import shutil
import subprocess

def gzip_except_metadata(folder_path, use_pigz=True):
    cmd = [
        "find", str(folder_path),
        "-type", "f", "!", "-name", "metadata.csv",
        "-exec", "pigz" if use_pigz else "gzip", "-f", "{}", ";"
    ]
    subprocess.run(cmd, check=True)
    print(f"✅ Compressed files (excluding metadata.csv) using {'pigz' if use_pigz else 'gzip'} with overwrite enabled.")

def split_and_organize_by_subtype(base_folder):
    print("starting...")
    base_folder = Path(base_folder)
    extracted_dir = base_folder / "extracted"
    sc_folder = extracted_dir / "Wu_etal_2021_BRCA_scRNASeq"
    metadata_path = sc_folder / "metadata.csv"

    if not metadata_path.exists():
        print("⚠️ metadata.csv not found.")
        return

    metadata_df = pd.read_csv(metadata_path)
    if not {"orig.ident", "subtype"}.issubset(metadata_df.columns):
        print("⚠️ metadata.csv missing required columns.")
        return

    bulk_file = next(extracted_dir.glob("*bulkRNAseq_raw_counts.txt"), None)
    if bulk_file is None:
        raise FileNotFoundError("Bulk RNA-seq file not found.")

    bulk_df = pd.read_csv(bulk_file, sep="\t", index_col=0)

    subtypes = metadata_df["subtype"].unique()

    for subtype in subtypes:
        print(f'Running {subtype}...')
        subtype_ids = metadata_df.query("subtype == @subtype")["orig.ident"].unique().tolist()
        subtype_str = subtype.replace(" ", "_")
        subtype_path = base_folder.parent / f"GSE176078_{subtype_str}"
        output_dir = subtype_path / "data_for_study"
        bulk_dir = output_dir / "bulk_rna_seq"
        sc_dir = output_dir / "single_cell"
        intermediate_dir = output_dir / "intermediate_data"

        for d in [bulk_dir, sc_dir, intermediate_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Filter and save bulk data
        filtered_bulk = bulk_df.loc[:, bulk_df.columns.isin(subtype_ids)]
        filtered_bulk.to_csv(bulk_dir / "combined_bulk.csv")

        # Filter metadata
        subtype_metadata = metadata_df[metadata_df["orig.ident"].isin(subtype_ids)]
        subtype_metadata.to_csv(sc_dir / "metadata.csv", index=False)

        # Copy and format scRNA files
        sc_file_map = {
            "count_matrix_genes.tsv": "features.tsv",
            "count_matrix_barcodes.tsv": "barcodes.tsv",
            "count_matrix_sparse.mtx": "matrix.mtx"
        }

        for src_name, dst_name in sc_file_map.items():
            src_file = sc_folder / src_name
            dst_file = sc_dir / dst_name
            if not src_file.exists():
                raise FileNotFoundError(f"{src_file} not found.")
            if dst_name == "features.tsv":
                with open(src_file, "r") as f:
                    lines = [line.strip() for line in f]
                with open(dst_file, "w") as f_out:
                    for i, gene_name in enumerate(lines):
                        gene_id = f"GENE{i+1:06d}"
                        f_out.write(f"{gene_id}\t{gene_name}\tGene Expression\n")
                print(f"✅ Rewritten: {src_name} → {dst_name} (2-column format)")
            else:
                shutil.copy(src_file, dst_file)
                print(f"✅ Copied: {src_name} → {dst_name}")

        # Gzip except metadata
        gzip_except_metadata(sc_dir, use_pigz=False)

        # Save subtype mapping
        subtype_df = (
            subtype_metadata[["orig.ident", "subtype"]]
            .drop_duplicates()
            .rename(columns={"orig.ident": "Sample", "subtype": "Subtype"})
        )
        subtype_df.to_csv(output_dir / "sample_subtype.csv", index=False)
        print(f"🎉 Finished organizing for subtype: {subtype}")
        
    print("📦 Saving full dataset (no subtype split)...")
    all_path = base_folder.parent / "GSE176078_all"
    all_output_dir = all_path / "data_for_study"
    all_bulk_dir = all_output_dir / "bulk_rna_seq"
    all_sc_dir = all_output_dir / "single_cell"
    all_intermediate_dir = all_output_dir / "intermediate_data"

    for d in [all_bulk_dir, all_sc_dir, all_intermediate_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Save full bulk data
    bulk_df.to_csv(all_bulk_dir / "combined_bulk.csv")

    # Save full metadata
    metadata_df.to_csv(all_sc_dir / "metadata.csv", index=False)

    # Copy and format scRNA files
    for src_name, dst_name in sc_file_map.items():
        src_file = sc_folder / src_name
        dst_file = all_sc_dir / dst_name
        if not src_file.exists():
            raise FileNotFoundError(f"{src_file} not found.")
        if dst_name == "features.tsv":
            with open(src_file, "r") as f:
                lines = [line.strip() for line in f]
            with open(dst_file, "w") as f_out:
                for i, gene_name in enumerate(lines):
                    gene_id = f"GENE{i+1:06d}"
                    f_out.write(f"{gene_id}\t{gene_name}\tGene Expression\n")
            print(f"✅ Rewritten: {src_name} → {dst_name} (2-column format)")
        else:
            shutil.copy(src_file, dst_file)
            print(f"✅ Copied: {src_name} → {dst_name}")

    gzip_except_metadata(all_sc_dir, use_pigz=False)

    # Save subtype mapping
    subtype_map_df = (
        metadata_df[["orig.ident", "subtype"]]
        .drop_duplicates()
        .rename(columns={"orig.ident": "Sample", "subtype": "Subtype"})
    )
    subtype_map_df.to_csv(all_output_dir / "sample_subtype.csv", index=False)
    print("🎉 Finished saving full dataset in GSE176078_all")

if __name__ == '__main__':
    split_and_organize_by_subtype("data/GSE176078")
