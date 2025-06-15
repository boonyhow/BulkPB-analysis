suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(edgeR)
  library(glue)
  library(dplyr)
  library(futile.logger)
  library(yaml)
  library(ggplot2)
  library(patchwork)
})

# Set up logger format (similar to Python's logging)
flog.layout(layout.format("[~l] ~t | ~m"))
flog.threshold(INFO)  # Adjust to DEBUG for more verbosity

read_data <- function(exp_id) {
  b <- read.table(glue("data/{exp_id}/data_for_study/bulk_rna_seq/combined_bulk.csv"),
                  row.names = 1, sep = ',', header = TRUE)
  
  sc_f <- Read10X(glue("data/{exp_id}/data_for_study/single_cell/"))
  sc <- CreateSeuratObject(sc_f, min.cells = 10, min.features = 200) ## Standard data filtering here
  
  metadata <- read.csv(glue("data/{exp_id}/data_for_study/single_cell/metadata.csv"), row.names = 1)
  sc <- AddMetaData(sc, metadata = metadata)
  
  pb <- AggregateExpression(sc, assays = "RNA", group.by = "orig.ident")$RNA
  pb <- pb[, colnames(b)]
  
  return(list(Bulk = list(raw = b), PseudoBulk = list(raw = pb)))
}

preprocessing_data <- function(data, gene_map_path = "data/gencode_data/gene_id_mapping.csv", non_zeros_required = 1) {
  ### non_zeros_required (int) ===> Ratio of samples needed for genes to be non-zero
  ###                               to be considered relevant in data 
  ###                               Maxmimum value is 1, minimum value is 0
  ###                               Default value is 1. 
  data <- data[rowSums(data > 0) >= (non_zeros_required * ncol(data)), , drop = FALSE]
  
  if (any(grepl("^ENSG", rownames(data)))) {
    flog.warn("Some ENSG gene IDs detected — performing gene name mapping where possible.")
  
    gene_map <- read.csv(gene_map_path, stringsAsFactors = FALSE)
    gene_map <- gene_map[gene_map$gene_id %in% rownames(data) & !is.na(gene_map$gene_name), ]
    
    data$gene_id <- rownames(data)
    merged_data <- merge(data, gene_map, by.x = "gene_id", by.y = "gene_id")
    
    processed_data <- merged_data %>%
      select(-gene_id) %>%
      group_by(gene_name) %>%
      summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
      as.data.frame()

    
    rownames(processed_data) <- processed_data$gene_name
    processed_data$gene_name <- NULL
    
    return(processed_data)
  } else {
    return(data)
  }
}


preprocess_all <- function(data, gene_map_path = "data/gencode_data/gene_id_mapping.csv", non_zeros_required = 1, mode = 'intersect') {
  flog.info("Preprocessing Bulk and PseudoBulk data")
  bulk_pre <- preprocessing_data(data$Bulk$raw, gene_map_path, non_zeros_required)
  pseudo_pre <- preprocessing_data(data$PseudoBulk$raw, gene_map_path, non_zeros_required)
  if (mode == 'union') {
    all_genes <- union(rownames(bulk_pre), rownames(pseudo_pre))
    flog.info(glue("Found {length(all_genes)} genes expressed in either Bulk or PseudoBulk"))

    # Fill in missing genes with zeros so matrices can be matched later
    bulk_pre <- bulk_pre[intersect(rownames(bulk_pre), all_genes), , drop = FALSE]
    pseudo_pre <- pseudo_pre[intersect(rownames(pseudo_pre), all_genes), , drop = FALSE]

    # Add missing genes with 0s (pad to full union)
    missing_bulk_genes <- setdiff(all_genes, rownames(bulk_pre))
    missing_pseudo_genes <- setdiff(all_genes, rownames(pseudo_pre))

    if (length(missing_bulk_genes) > 0) {
      bulk_pre <- rbind(bulk_pre, matrix(0, nrow = length(missing_bulk_genes), ncol = ncol(bulk_pre),
                                        dimnames = list(missing_bulk_genes, colnames(bulk_pre))))
    }

    if (length(missing_pseudo_genes) > 0) {
      pseudo_pre <- rbind(pseudo_pre, matrix(0, nrow = length(missing_pseudo_genes), ncol = ncol(pseudo_pre),
                                            dimnames = list(missing_pseudo_genes, colnames(pseudo_pre))))
    }

    # Reorder rows to ensure same gene order
    bulk_pre <- bulk_pre[all_genes, , drop = FALSE]
    pseudo_pre <- pseudo_pre[all_genes, , drop = FALSE]
  
  } else if (mode == 'intersect') {
    common_genes <- intersect(rownames(bulk_pre), rownames(pseudo_pre))
    flog.info(glue("Found {length(common_genes)} genes shared between Bulk and PseudoBulk"))
    if (length(common_genes) < 0.5 * dim(bulk_pre)[1]) {
      flog.warn("Preprocessing mode 'intersect' caused the number of genes retained to be less than 50% of the original number of genes in bulk, consider mode 'union' if applicable, or lower the non_zeros_required parameter in preprocessing.")
    }
    if (length(common_genes) < 0.5 * dim(pseudo_pre)[1]) {
      flog.warn("Preprocessing mode 'intersect' caused the number of genes retained to be less than 50% of the original number of genes in pseudo-bulk, consider mode 'union' if applicable, or lower the non_zeros_required parameter in preprocessing.")
    }
    bulk_pre <- bulk_pre[common_genes, , drop = FALSE]
    pseudo_pre <- pseudo_pre[common_genes, , drop = FALSE]
      
  } else {
    flog.error('Mode in data preprocessing step is not well defined, stopping the code.')
    stop()
  }
  # keep_genes <- rowSums(edgeR::cpm(bulk_pre) > 1) | rowSums(edgeR::cpm(pseudo_pre) > 1)

  data$Bulk$preprocessed <- round(bulk_pre)
  data$PseudoBulk$preprocessed <- round(pseudo_pre)
  return(data)
}

normalise_data <- function(data) {
  
  norm_bulk <- log1p(edgeR::cpm(calcNormFactors(DGEList(counts = data$Bulk$preprocessed), method = "TMM")))
  norm_pseudo <- log1p(edgeR::cpm(calcNormFactors(DGEList(counts = data$PseudoBulk$preprocessed), method = "TMM")))
  
  data$Bulk$normalised <- norm_bulk
  data$PseudoBulk$normalised <- norm_pseudo
  return(data)
}


main <- function(config_path = "configs/1_config.yaml") {
  if (!file.exists(config_path)) {
    flog.error(glue("Config file '{config_path}' not found."))
    stop()
  }
  
  config <- yaml::read_yaml(config_path)
  exp_id <- config$exp_id
  mode <- config$preprocessing$mode %||% "intersect"
  non_zeros_required <- config$preprocessing$non_zeros_required %||% 1

  file_path <- config$save_file_path %||% '1_preprocessed.RData'

  data_dir <- glue("data/{exp_id}/data_for_study")
  if (!dir.exists(data_dir)) {
    flog.error(glue("Experiment ID '{exp_id}' is not valid. Directory {data_dir} not found."))
    stop()
  }
  
  flog.info(glue("Starting preprocessing for experiment ID: {exp_id}"))
  
  flog.info(glue("Reading data for experiment: {exp_id}"))
  data_obj <- read_data(exp_id)


  flog.info(glue("Preprocessing data with mode='{mode}', non_zeros_required={non_zeros_required}"))
  data_obj <- preprocess_all(data_obj, non_zeros_required = non_zeros_required, mode = mode)
  flog.info("Normalising data using log-CPM (TMM)")
  data_obj <- normalise_data(data_obj)

  output_dir <- glue("{data_dir}/intermediate_data")
  if (!dir.exists(output_dir)) {
    flog.info(glue("Creating output directory: {output_dir}"))
    dir.create(output_dir, recursive = TRUE)
  }
  
  save_path <- glue("{output_dir}/{file_path}")
  saveRDS(data_obj, file = save_path)
  flog.info(glue("Saved preprocessed data to: {save_path}"))
}

# === Entry Point for CLI ===
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "1_config.yaml"
main(config_path)