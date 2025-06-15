import pandas as pd
import numpy as np
import pyarrow.csv as pv
import scipy.io
import scipy.sparse
import json
import re
from pathlib import Path
from tqdm import tqdm
import gzip
import shutil as sh

def organize_extracted_data_GSE217517(main_folder):
    main_dir = Path(main_folder)
    extracted_folder = main_dir / 'extracted'
    output_base = main_dir / 'data_for_study'
    output_base.mkdir(parents=True, exist_ok=True)

    # Create folders
    bulk_dir = output_base / 'bulk_rna_seq'
    sc_dir = output_base / 'single_cell'
    bulk_dir.mkdir(parents=True, exist_ok=True)
    sc_dir.mkdir(parents=True, exist_ok=True)

    # Regular expressions
    id_pattern = re.compile(r'_(\d{4})')
    gsm_pattern = re.compile(r'(GSM\d{6,})')

    bulk_files = []
    single_cell_files = []

    # Collect files
    for file in extracted_folder.glob('*'):
        if file.is_file():
            filename = file.name
            if 'bulk_dissociated_polyA' in filename:
                bulk_files.append(file)
            elif 'single_cell' in filename and 'pooled_single_cell' not in filename:
                single_cell_files.append(file)

    # Build mapping: GSM -> 4-digit sample ID
    gsm_to_sample_id = {}
    for file in bulk_files + single_cell_files:
        filename = file.name
        gsm_match = gsm_pattern.search(filename)
        id_match = id_pattern.search(filename)
        if gsm_match and id_match:
            gsm_id = gsm_match.group(1)
            sample_id = id_match.group(1)
            gsm_to_sample_id[gsm_id] = sample_id

    # Sort GSMs, assign P1, P2...
    sorted_gsm_sample_pairs = sorted(gsm_to_sample_id.items(), key=lambda x: int(x[0].replace('GSM', '')))
    id_mapping = {}
    current_p_number = 1
    for gsm_id, sample_id in sorted_gsm_sample_pairs:
        if sample_id not in id_mapping:
            id_mapping[sample_id] = f"P{current_p_number}"
            current_p_number += 1

    # === Bulk RNA-seq Processing ===
    unmapped_ratios = {}
    bulk_data = {}

    for file in tqdm(bulk_files, desc="Processing bulk RNA-seq"):
        gsm_id = gsm_pattern.search(file.name).group(1)
        sample_id = gsm_to_sample_id[gsm_id]
        p_number = id_mapping[sample_id]
        
        # Read bulk file (no header)
        parse_options = pv.ParseOptions(delimiter="\t")
        df = pv.read_csv(file, parse_options=parse_options, read_options=pv.ReadOptions(autogenerate_column_names=True)).to_pandas()
        df.columns = ['Gene'] + [f"Rep_{i}" for i in range(1, df.shape[1])]

        # Calculate unmapped ratio
        unmapped_row = df[df['Gene'] == 'N_noFeature']
        if unmapped_row.empty:
            raise ValueError(f"N_noFeature row not found in {file.name}")
        unmapped_counts = unmapped_row.iloc[0, 1:].astype(float)
        unmapped_ratios[p_number] = unmapped_counts.to_dict()

        # Keep only ENSG rows
        df = df[df['Gene'].str.startswith('ENSG')].reset_index(drop=True)

        bulk_data[p_number] = df


    # Identify and drop bad replicates and optionally entire samples
    samples_to_keep = {}
    justification_text = ""
    sample_level_drops = ""

    # First, replicate-level filtering
    replicate_std_devs = {}

    for sample_id, counts in unmapped_ratios.items():
        count_values = np.array(list(counts.values()))
        median_count = np.median(count_values)
        stdev_count = np.std(count_values)
        replicate_std_devs[sample_id] = stdev_count

        bad_replicates = []
        for rep_name, count in counts.items():
            if count > 3 * median_count or count < median_count / 3:
                bad_replicates.append(rep_name)

        if bad_replicates:
            justification_text += f"Sample {sample_id}: Dropped replicates {bad_replicates} due to N_noFeature count deviation.\n"

        good_replicates = [col for col in bulk_data[sample_id].columns if col != 'Gene' and col.replace("Rep_", "") not in [r.split("_")[1] for r in bad_replicates]]
        good_df = bulk_data[sample_id][['Gene'] + good_replicates]
        samples_to_keep[sample_id] = good_df

    # Then, sample-level filtering (based on stdev of unmapped ratios)
    median_std = np.median(list(replicate_std_devs.values()))
    threshold_std = 2 * median_std

    final_samples = {}
    for sample_id, df in samples_to_keep.items():
        if replicate_std_devs[sample_id] > threshold_std:
            sample_level_drops += f"Sample {sample_id} dropped entirely due to high variance in N_noFeature counts (stddev = {replicate_std_devs[sample_id]:.2f}, median = {median_std:.2f}).\n"
            continue
        final_samples[sample_id] = df

    samples_to_keep = final_samples


    # Merge and average replicates
    merged_bulk = None
    for sample_id, df in samples_to_keep.items():
        averaged = df.set_index('Gene').mean(axis=1)
        averaged.name = sample_id
        if merged_bulk is None:
            merged_bulk = averaged.to_frame()
        else:
            merged_bulk = merged_bulk.join(averaged, how='inner')

    # Round and convert to integer
    merged_bulk = np.round(merged_bulk).astype(int)

    # Sort columns as P1, P2, P3, ...
    merged_bulk = merged_bulk[sorted(merged_bulk.columns)]

    merged_bulk.to_csv(bulk_dir / "combined_bulk.csv")
    print("Saved combined bulk RNA-seq.")


    # === Single-cell RNA-seq Processing ===
    feature_list = []
    barcode_list = []
    count_matrices = []

    for file in tqdm(single_cell_files, desc="Processing single-cell RNA-seq"):
        gsm_id = gsm_pattern.search(file.name).group(1)
        sample_id = gsm_to_sample_id[gsm_id]
        p_number = id_mapping[sample_id]
        
        if 'features' in file.name:
            df = pd.read_csv(file, header=None, sep="\t")
            df.columns = ['feature_id', 'feature_name', 'feature_type']
            feature_list.append(df)
        elif 'barcodes' in file.name:
            df = pd.read_csv(file, header=None)
            df[0] = df[0].apply(lambda x: f"{p_number}_{x}")
            barcode_list.append(df)
        elif 'matrix' in file.name:
            mat = scipy.io.mmread(file)
            count_matrices.append(mat)

    merged_counts = scipy.sparse.hstack(count_matrices)
    merged_barcodes = pd.concat(barcode_list, ignore_index=True)
    final_features = feature_list[0]

    # Save single-cell files
    scipy.io.mmwrite(str(sc_dir / "matrix.mtx"), merged_counts)
    final_features.to_csv(sc_dir / "features.tsv", sep="\t", header=False, index=False)
    merged_barcodes.to_csv(sc_dir / "barcodes.tsv", header=False, index=False)

    # # Gzip them
    # for filename in ["matrix.mtx", "features.tsv", "barcodes.tsv"]:
    #     with open(sc_dir / filename, 'rb') as f_in:
    #         with gzip.open(sc_dir / f"{filename}.gz", 'wb') as f_out:
    #             sh.copyfileobj(f_in, f_out)
    #     (sc_dir / filename).unlink()

    # Metadata
    metadata = pd.DataFrame({
        "barcode": merged_barcodes[0],
        "orig.ident": merged_barcodes[0].apply(lambda x: x.split('_')[0])
    })
    metadata.to_csv(sc_dir / "metadata.csv", index=False)
    print("Saved merged single-cell RNA-seq.")

    # Save mapping
    with open(output_base / 'sample_id_mapping.json', 'w') as f:
        json.dump(id_mapping, f, indent=4)

    with open(output_base / 'README.txt', 'w') as f:
        f.write("Sample IDs have been renamed to P1, P2, etc., sorted by GSM ID.\n")
        f.write("Unmapped ratio deviations >3-fold resulted in replicate removal.\n")
        f.write("Samples with high overall variance (stddev > 2× median across samples) were removed entirely.\n\n")

        f.write("==== Unmapped Ratios by Sample ====\n")
        for sample_id, ratios in unmapped_ratios.items():
            f.write(f"{sample_id}:\n")
            for rep, ratio in ratios.items():
                f.write(f"  {rep}: {ratio:.5f}\n")

        f.write("\n==== Dropped Replicates ====\n")
        f.write(justification_text)

        f.write("\n==== Dropped Samples Due to High Variance ====\n")
        f.write(sample_level_drops if sample_level_drops else "None\n")
    print("Saved mapping and README.")
    
    # === Save Sample Subtype Table ===
    subtype_df = pd.DataFrame({
        "Sample": sorted(samples_to_keep.keys(), key=lambda x: int(x[1:])),  # Sort as P1, P2, ...
        "Subtype": ["HGSOC"] * len(samples_to_keep)
    })
    subtype_df.to_csv(output_base / "sample_subtype.csv", index=False)
    print("Saved sample_subtype.csv.")


if __name__ == "__main__":
    organize_extracted_data_GSE217517(
        "data/GSE217517"
    )
