suppressPackageStartupMessages({
    library(Seurat)
    library(oposSOM)
    library(pROC)
    library(sparsepca)
    library(limma)
    library(philentropy)    # For Jensen-Shannon Divergence
    library(matrixStats)    # For row/col stats
    library(transport)      # For Optimal Transport
    library(stats)          # For Mahalanobis
    library(glue)
    library(dplyr)
    library(futile.logger)
    library(DESeq2)
    library(edgeR)
    library(ggplot2)
    library(pheatmap)
    library(ComplexHeatmap)
    library(ggrepel)
    library(circlize)
    library(pbmcapply)
    library(patchwork)
    library(parallel)
    library(topGO)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(tidyverse)
    library(future)
    library(future.apply)
    library(progressr)
    library(tidyr)
})

### UTITLITY 
save_named_plot_list <- function(
  plot_list,
  folder_path,
  width = 8,
  height = 6,
  dpi = 300,
  format = "png",
  custom_sizes = list()  # named list: list(plot_name = c(width, height))
) {
  # Ensure folder exists
  if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

  # Save each plot with custom or default size
  purrr::iwalk(plot_list, function(p, name) {
    dims <- custom_sizes[[name]]
    this_width <- if (!is.null(dims)) dims[1] else width
    this_height <- if (!is.null(dims)) dims[2] else height

    file_path <- file.path(folder_path, paste0(name, ".", format))
    ggsave(filename = file_path, plot = p, width = this_width, height = this_height, dpi = dpi, units = "in")
  })
}


#  ============================
#  =         CODES FOR        =
#  =        PCA ANALYSIS      =
#  ============================
test_spca_separation <- function(spca_layer_obj, method = c("t.test", "wilcox"), p.adjust.method = "BH") {
  method <- match.arg(method)

  # Extract data
  scores <- spca_layer_obj$pca_obj$scores  # matrix of sample x PC scores
  rownames(scores) <- rownames(spca_layer_obj$pca_obj$scores)
  pca_data <- spca_layer_obj$pca_data       # contains Type ("Bulk" or "Pseudo")

  # Ensure order matches
  sample_types <- pca_data$Type[match(rownames(scores), pca_data$Data_Name)]

  # Run tests on each PC
  pvals <- apply(scores, 2, function(pc_scores) {
    df <- data.frame(Score = pc_scores, Type = sample_types)
    if (length(unique(df$Type)) < 2) return(NA)  # skip if only one group
    test_result <- switch(method,
                          t.test = t.test(Score ~ Type, data = df),
                          wilcox = wilcox.test(Score ~ Type, data = df))
    test_result$p.value
  })

  # Adjust for multiple testing
  adj_pvals <- p.adjust(pvals, method = p.adjust.method)

  # Combine results
  result_df <- data.frame(
    PC = paste0("PC", seq_along(pvals)),
    pval = pvals,
    adj_pval = adj_pvals
  )

  return(result_df)
}


run_all_sparse_pca <- function(data_obj, subtype_df = NULL, k = 10, alpha = 1e-3, beta = 1e-3, max_iter = 500) {
  flog.info("Running Sparse PCA (SparsePCA package) for all layers dynamically...")

  get_layer_spca <- function(bulk_layer, pseudo_layer, layer_name, subtype_df = NULL) {
    stopifnot(all(colnames(bulk_layer) == colnames(pseudo_layer)))
    common_genes <- intersect(rownames(bulk_layer), rownames(pseudo_layer))
    bulk_layer <- bulk_layer[common_genes, , drop = FALSE]
    pseudo_layer <- pseudo_layer[common_genes, , drop = FALSE]

    colnames(bulk_layer) <- paste0(colnames(bulk_layer), "_Bulk")
    colnames(pseudo_layer) <- paste0(colnames(pseudo_layer), "_PseudoBulk")

    expr <- cbind(bulk_layer, pseudo_layer)

    # Remove zero-variance genes
    gene_vars <- apply(expr, 1, var)
    if (any(gene_vars == 0)) {
      flog.warn("Layer '%s' has %d zero-variance genes. Removing them before Sparse PCA.", 
                layer_name, sum(gene_vars == 0))
      expr <- expr[gene_vars > 0, , drop = FALSE]
    }

    if (nrow(expr) == 0) {
      flog.error("All genes have zero variance in layer '%s'. Skipping Sparse PCA.", layer_name)
      return(NULL)
    }

    # Run Sparse PCA
    spca_result <- spca(
      X = t(expr),      # samples x genes
      k = k,
      alpha = alpha,
      beta = beta,
      center = TRUE,
      scale = TRUE,
      max_iter = max_iter
    )

    all_samples <- colnames(expr)
    spca_df <- data.frame(
      Data_Name = all_samples,
      PC1 = spca_result$scores[, 1],
      PC2 = spca_result$scores[, 2],
      Type = ifelse(grepl("_PseudoBulk$", all_samples), "Pseudo", "Bulk"),
      Sample = sub("_.*$", "", all_samples),
      data_type = layer_name
    )

    if (!is.null(subtype_df)) {
      spca_df <- merge(spca_df, subtype_df, by = "Sample", all.x = TRUE)
    }

    return(list(pca_obj = spca_result, pca_data = spca_df))
  }

  collect_spca <- function(bulk_list, pseudo_list, layer_name, subtype_df = NULL) {
    if (is.matrix(bulk_list) || is.data.frame(bulk_list)) {
      get_layer_spca(bulk_list, pseudo_list, layer_name, subtype_df)
    } else if (is.list(bulk_list)) {
      common_subkeys <- intersect(names(bulk_list), names(pseudo_list))
      result <- lapply(common_subkeys, function(subkey) {
        get_layer_spca(bulk_list[[subkey]], pseudo_list[[subkey]], subkey, subtype_df)
      })
      names(result) <- common_subkeys
      return(result)
    }
  }

  bulk <- data_obj$Bulk
  pseudo <- data_obj$PseudoBulk
  shared_layers <- setdiff(intersect(names(bulk), names(pseudo)), "batch_corrected")

  spca_result <- lapply(shared_layers, function(layer) {
    collect_spca(bulk[[layer]], pseudo[[layer]], layer, subtype_df)
  })
  names(spca_result) <- shared_layers

  return(spca_result)
}


sparse_pca_to_de <- function(spca_obj, pc = 1, threshold = 0, column_name = "logFC") {
  loadings <- spca_obj$loadings[, pc]
  gene_names <- names(spca_obj$center)

  df <- data.frame(
    gene = gene_names,
    sig = ifelse(abs(loadings) > threshold, "Significant",  "Not Significant"),
    stringsAsFactors = FALSE
  )
  df[[column_name]] <- loadings
  rownames(df) <- df$gene

  return(df)
}

run_all_pca <- function(data_obj, subtype_df = NULL) {
  flog.info("Running PCA for all layers dynamically...")

  get_layer_pca <- function(bulk_layer, pseudo_layer, layer_name, subtype_df = NULL) {
    stopifnot(all(colnames(bulk_layer) == colnames(pseudo_layer)))
    common_genes <- intersect(rownames(bulk_layer), rownames(pseudo_layer))
    bulk_layer <- bulk_layer[common_genes, , drop = FALSE]
    pseudo_layer <- pseudo_layer[common_genes, , drop = FALSE]

    colnames(bulk_layer) <- paste0(colnames(bulk_layer), "_Bulk")
    colnames(pseudo_layer) <- paste0(colnames(pseudo_layer), "_PseudoBulk")
    
    expr <- cbind(bulk_layer, pseudo_layer)

    # Check for zero-variance genes
    gene_vars <- apply(expr, 1, var)
    zero_var_genes <- sum(gene_vars == 0)
    if (zero_var_genes > 0) {
      flog.warn("Layer '%s' has %d zero-variance genes. Removing them before PCA.", layer_name, zero_var_genes)
      expr <- expr[gene_vars > 0, , drop = FALSE]
    }

    if (nrow(expr) == 0) {
      flog.error("All genes have zero variance in layer '%s'. Skipping PCA.", layer_name)
      return(NULL)
    }

    # Run PCA
    pca <- prcomp(t(expr), scale. = TRUE)

    all_samples <- colnames(expr)
    pca_df <- data.frame(
      Data_Name = all_samples,
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      Type = ifelse(grepl("_PseudoBulk$", all_samples), "Pseudo", "Bulk"),
      Sample = sub("_.*$", "", all_samples),
      data_type = layer_name
    )

    # Optional: merge in subtype if provided
    if (!is.null(subtype_df)) {
      pca_df <- merge(pca_df, subtype_df, by = "Sample", all.x = TRUE)
    }


    return(list(pca_obj = pca, pca_data = pca_df))
  }


  collect_pca <- function(bulk_list, pseudo_list, layer_name, subtype_df = NULL) {
    result <- list()

    # Case: direct matrix/data.frame (e.g., "normalised")
    if (is.matrix(bulk_list) || is.data.frame(bulk_list)) {

      result <- get_layer_pca(bulk_list, pseudo_list, layer_name, subtype_df)
    
    # Case: nested list (e.g., "batch_corrected")
    } else if (is.list(bulk_list)) {
      common_subkeys <- intersect(names(bulk_list), names(pseudo_list))
      result <- lapply(common_subkeys, function(subkey) {
        get_layer_pca(bulk_list[[subkey]], pseudo_list[[subkey]], subkey, subtype_df)
      })
      names(result) <- common_subkeys
    }

    return(result)
  }

  bulk <- data_obj$Bulk
  pseudo <- data_obj$PseudoBulk
  shared_layers <- setdiff(intersect(names(bulk), names(pseudo)), "raw")

  pca_result <- lapply(shared_layers, function(layer) {
    collect_pca(bulk[[layer]], pseudo[[layer]], layer, subtype_df)
  })
  names(pca_result) <- shared_layers

  return(pca_result)
}



#  ============================
#  =         CODES FOR        =
#  =        KDE ANALYSIS      =
#  ============================

run_all_kde <- function(data_obj, use_all_samples = TRUE) {
  flog.info("Running KDE for all layers dynamically...")

  extract_kde_data <- function(bulk_layer, pseudo_layer, layer_name, use_all_samples) {
    common_genes <- intersect(rownames(bulk_layer), rownames(pseudo_layer))
    bulk_layer <- bulk_layer[common_genes, , drop = FALSE]
    pseudo_layer <- pseudo_layer[common_genes, , drop = FALSE]

    if (!use_all_samples) {
      bulk_data <- rowMeans(bulk_layer)
      pseudo_data <- rowMeans(pseudo_layer)
      df <- data.frame(
        Gene = rep(common_genes, 2),
        Expression = c(bulk_data, pseudo_data),
        Modality = rep(c("Bulk", "PseudoBulk"), each = length(bulk_data)),
        Layer = layer_name
      )
    } else {
      colnames(bulk_layer) <- paste0(colnames(bulk_layer), "_Bulk")
      colnames(pseudo_layer) <- paste0(colnames(pseudo_layer), "_PseudoBulk")
      combined <- cbind(bulk_layer, pseudo_layer)

      df <- reshape2::melt(combined, variable.name = "Sample", value.name = "Expression")
      df$Gene <- rep(rownames(combined), times = ncol(combined))
      df$Modality <- ifelse(grepl("PseudoBulk$", df$Sample), "PseudoBulk", "Bulk")
      df$Layer <- layer_name
    }

    return(df)
  }

  collect_kde <- function(bulk_list, pseudo_list, layer_name, use_all_samples) {
    if (is.matrix(bulk_list) || is.data.frame(bulk_list)) {
      return(extract_kde_data(bulk_list, pseudo_list, layer_name, use_all_samples))
    } else if (is.list(bulk_list)) {
      common_subkeys <- intersect(names(bulk_list), names(pseudo_list))
      dfs <- lapply(common_subkeys, function(subkey) {
        extract_kde_data(bulk_list[[subkey]], pseudo_list[[subkey]], subkey, use_all_samples)
      })
      return(do.call(rbind, dfs))
    }
  }

  bulk <- data_obj$Bulk
  pseudo <- data_obj$PseudoBulk
  shared_layers <- setdiff(intersect(names(bulk), names(pseudo)), "raw")

  kde_df_list <- lapply(shared_layers, function(layer) {
    collect_kde(bulk[[layer]], pseudo[[layer]], layer, use_all_samples)
  })
  names(kde_df_list) <- shared_layers

  return(do.call(rbind, kde_df_list))
}



#  ============================
#  =         CODES FOR        =
#  =        COR ANALYSIS      =
#  ============================

run_all_correlation <- function(data_obj, samplewise_correlation = TRUE) {
  flog.info("Running Pearson Correlations for all layers dynamically...")

  extract_correlation <- function(bulk_layer, pseudo_layer, layer_name, samplewise_correlation) {
    common_genes <- intersect(rownames(bulk_layer), rownames(pseudo_layer))
    bulk_layer <- bulk_layer[common_genes, , drop = FALSE]
    pseudo_layer <- pseudo_layer[common_genes, , drop = FALSE]

    if (!samplewise_correlation) {
      flatten_corr_list <- function(corr_list, layer_name) {
        do.call(rbind, lapply(names(corr_list), function(type_label) {
          data.frame(
            Layer = layer_name,
            Type = type_label,
            Correlation = corr_list[[type_label]]
          )
        }))
      }

      # BULK vs BULK
      bulk_corr <- cor(bulk_layer)

      # PB vs PB
      pb_corr <- cor(as.data.frame(pseudo_layer))

      # BULK vs PB
      colnames(bulk_layer) <- paste0(colnames(bulk_layer), "_Bulk")
      colnames(pseudo_layer) <- paste0(colnames(pseudo_layer), "_PB")
      combined <- cbind(bulk_layer, pseudo_layer)
      all_corr <- cor(combined)

      # Extract only Bulk vs PB correlations
      bulk_cols <- grep("_Bulk$", colnames(combined), value = TRUE)
      pb_cols <- grep("_PB$", colnames(combined), value = TRUE)
      cross_corr <- all_corr[bulk_cols, pb_cols]

      tmp <- list(
        Bulk_vs_Bulk = bulk_corr[upper.tri(bulk_corr)],
        PB_vs_PB = pb_corr[upper.tri(pb_corr)],
        Bulk_vs_PB = as.vector(cross_corr)
      )

      return(flatten_corr_list(tmp, layer_name))
    } else {
      colnames(bulk_layer) <- paste0(colnames(bulk_layer), "_Bulk")
      colnames(pseudo_layer) <- paste0(colnames(pseudo_layer), "_PseudoBulk")
      combined <- cbind(bulk_layer, pseudo_layer)

      df <- cor(combined)
    }

    return(df)
  }

  collect_correlation <- function(bulk_list, pseudo_list, layer_name, samplewise_correlation) {
    if (is.matrix(bulk_list) || is.data.frame(bulk_list)) {
      return(extract_correlation(bulk_list, pseudo_list, layer_name, samplewise_correlation))
    } else if (is.list(bulk_list)) {
      common_subkeys <- intersect(names(bulk_list), names(pseudo_list))
      dfs <- lapply(common_subkeys, function(subkey) {
        extract_correlation(bulk_list[[subkey]], pseudo_list[[subkey]], subkey, samplewise_correlation)
      })
      names(dfs) <- common_subkeys
      return(dfs)
    }
  }

  bulk <- data_obj$Bulk
  pseudo <- data_obj$PseudoBulk
  shared_layers <- setdiff(intersect(names(bulk), names(pseudo)), "raw")

  corr_df_list <- lapply(shared_layers, function(layer) {
    collect_correlation(bulk[[layer]], pseudo[[layer]], layer, samplewise_correlation)
  })
  names(corr_df_list) <- shared_layers

  return(corr_df_list)
}



#  ============================
#  =         CODES FOR        =
#  =        DE ANALYSIS       =
#  ============================

make_metadata <- function(counts, subtypes , subtype_df) {
  samples <- colnames(counts)
  
  # Extract Individual and Modality from sample names
  indiv <- sub("(_Bulk.*|_PseudoBulk.*)$", "", samples)

  modality <- ifelse(grepl("Original$", samples), "Original", "ComBatSeq")
  data_type <- ifelse(grepl("_PseudoBulk", samples), "PseudoBulk", "Bulk")

  md <- data.frame(
    Sample = samples,
    Individual = indiv,
    Modality = modality,
    DataType = data_type,
    stringsAsFactors = FALSE
  )


  # Join with subtype metadata
  md <- left_join(md, subtype_df, by = c("Individual" = "Sample"))

  if (!is.null(subtypes)) {
    # Filter to only include specified subtypes
    md <- md[md$Subtype %in% subtypes, ]
  }
  

  # md$Subtype <- factor(make.names(md$Subtype), levels = make.names(subtypes))
  rownames(md) <- md$Sample
  return(md)
}

# Generate full expression matrix
get_combined_counts <- function(data_list) {

  rename_cols <- function(mat, suffix) {
    colnames(mat) <- paste0(colnames(mat), suffix)
    mat
  }
    # Generalised function to get renamed count matrices
  get_counts <- function(data_list, modality) {
   
    path <- data_list[[modality]]
  
    path <- path$preprocessed
    suffix <- paste0("_", modality, "_Original")
  
    rename_cols(path, suffix)
  }

  bulk <- get_counts(data_list,  "Bulk")
  pseudo <- get_counts(data_list,  "PseudoBulk")
  round(cbind(bulk, pseudo))
}


run_deseq2_de_pipeline <- function(counts, subtype_df,subtype = NULL,  metadata_column = "DataType") {
  # Build metadata
  meta <- make_metadata(counts, subtype, subtype_df)
  # group <- factor(meta[[metadata_column]])
  meta[[metadata_column]] <- factor(meta[[metadata_column]])
  # design <- model.matrix(group)
  # colnames(design) <- levels(group)
  # Construct DESeq2 object with dynamic design
  design_formula <- as.formula(paste("~", metadata_column))
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = meta[colnames(counts), ],
                                 design = design_formula)

  dds <- DESeq(dds)

  # Get default contrast: second vs first level of selected metadata column
  levels_vec <- levels(meta[[metadata_column]])
  if (length(levels_vec) < 2) stop("Need at least 2 groups in metadata for DESeq2.")

  res <- results(dds, contrast = c(metadata_column, levels_vec[2], levels_vec[1]))
  res <- as.data.frame(res)

  # Annotate
  res <- res %>%
    mutate(
      logFC = log2FoldChange,
      negLogP = -log10(padj),
      sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant")
    )

  return(res)
}


#  ============================
#  =        CODES FOR         =
#  =       SPCA AND DE        =
#  ============================


summarise_spca_vs_de_overlap <- function(results_list, spca_list, pc = 1, column_name = "logFC", de_method = "deseq2_res") {
  all_results <- list()
  updated_de_tables <- list()

  for (exp_id in names(results_list)) {
    if (!de_method %in% names(results_list[[exp_id]]$de)) {
      warning(glue::glue("⚠️ Skipping {exp_id}: DE method '{de_method}' not found."))
      next
    }

    de_res <- results_list[[exp_id]]$de[[de_method]]
    spca_obj <- spca_list[[exp_id]]$normalised$pca_obj

    # Convert SPCA to DE-like dataframe
    spca_df <- sparse_pca_to_de(spca_obj, pc = pc, column_name = column_name)

    # Match rownames for SPCA-DE comparison
    common_genes <- intersect(rownames(de_res), rownames(spca_df))

    # Add SPCA-DE overlap column
    spca_sig_genes <- rownames(subset(spca_df, sig == "Significant"))
    de_sig_genes <- rownames(subset(de_res, sig == "Significant"))

    both_sig <- intersect(spca_sig_genes, de_sig_genes)

    # Create new column in DE result
    de_res$spca_de_sig <- ifelse(rownames(de_res) %in% both_sig, "Significant", "Not Significant")

    # Save updated DE result
    updated_de_tables[[exp_id]] <- de_res

    # Now generate the SPCA-DE overlap metrics per direction
    for (direction in c("Bulk", "PseudoBulk")) {
      spca_genes <- rownames(spca_df[
        spca_df$sig == "Significant" &
          ((direction == "Bulk" & spca_df[[column_name]] < 0) |
           (direction == "PseudoBulk" & spca_df[[column_name]] > 0)),
      ])

      de_genes <- get_intersect_significant_genes(de_res, direction = direction)

      eval_df <- evaluate_gene_overlap(de_genes, spca_genes)
      eval_df$Direction <- paste(direction, "Upregulated")
      eval_df$Experiment <- exp_id
      eval_df$DEMethod <- de_method

      all_results[[paste(exp_id, direction, sep = "_")]] <- eval_df
    }
  }

  result_df <- do.call(rbind, all_results)
  rownames(result_df) <- NULL

  return(list(
    summary_df = result_df,
    updated_de = updated_de_tables
  ))
}



evaluate_gene_overlap <- function(de_genes, spca_genes) {
  tp <- length(intersect(de_genes, spca_genes))
  fp <- length(setdiff(spca_genes, de_genes))
  fn <- length(setdiff(de_genes, spca_genes))

  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0

  data.frame(
    DE_Genes = length(de_genes),
    SPCA_Genes = length(spca_genes),
    Overlap = tp,
    Precision = round(precision, 3),
    Recall = round(recall, 3),
    F1_Score = round(f1, 3)
  )
}



get_intersect_significant_genes <- function(de_result, sig_label = "Significant", direction = NULL, logfc_col = "logFC") {
  if (is.null(de_result)) return(character(0))
  filtered <- subset(de_result, sig == sig_label)
  
  if (!is.null(direction)) {
    if (!logfc_col %in% colnames(filtered)) return(character(0))
    
    if (direction == "PseudoBulk") {
      filtered <- filtered[filtered[[logfc_col]] > 0, ]
    } else if (direction == "Bulk") {
      filtered <- filtered[filtered[[logfc_col]] < 0, ]
    }
  }
  
  rownames(filtered)
}


run_de_pipelines <- function(data_obj, data_subtype_df, run_go = T) {
  # Step 1: Get combined counts
  counts <- get_combined_counts(data_obj)

    # Step 2: Run DE pipelines (sequential)
  de_results <- list(
    deseq2_res = run_deseq2_de_pipeline(counts, data_subtype_df)
  )

  # Step 2b: Volcano plots
  volcano_plots <- list(
    deseq2_res = plot_volcano(de_results$deseq2_res, logfc_col = 'log2FoldChange', padj_col = 'padj')
  )

  if (run_go) {
    # Step 3: Run GO analyses in parallel with message logs
    message("⚙️ Starting GO analyses in parallel...")
    go_results <- future_lapply(names(de_results), function(method) {
      message(glue("🧬 Running GO: {method}"))
      run_go_by_direction(de_results[[method]])
    })
  } else {
    go_results <- lapply(names(de_results), function(method) {
     NULL
    })
  }
  # Preserve names
  names(go_results) <- names(de_results)

  return(list(
    de = de_results,
    go = go_results,
    volcano_plots = volcano_plots
  ))
}


#  ============================
#  =         CODES FOR        =
#  =        GO ANALYSIS       =
#  ============================

# Wrapper to run GO analysis for up and down genes
run_go_by_direction <- function(de_result,
                                logFC_threshold = 1,
                                sig_col = "sig",
                                logfc_col = "logFC",
                                sig_label = "Significant") {
  safe_go <- function(gene_set, bg_genes, direction_label) {
    tryCatch({
      run_go_analysis_pipeline(gene_set, bg_genes)
    }, error = function(e) {
      message(glue("⚠️ GO analysis failed for {direction_label}: {e$message}"))
      return(NULL)
    })
  }

  # Filter based on custom column names
  up_genes <- rownames(subset(de_result,
                              !!de_result[[sig_col]] == sig_label &
                              !!de_result[[logfc_col]] > logFC_threshold))

  down_genes <- rownames(subset(de_result,
                                !!de_result[[sig_col]] == sig_label &
                                !!de_result[[logfc_col]] < -logFC_threshold))

  all_genes <- rownames(de_result)

  list(
    pseudo_up = safe_go(up_genes, all_genes, "pseudo_up"),
    bulk_up   = safe_go(down_genes, all_genes, "bulk_up")
  )
}




# Function to convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_list) {
  entrez_ids <- mapIds(
    org.Hs.eg.db, 
    keys = gene_list, 
    column = "ENTREZID", 
    keytype = "SYMBOL", 
    multiVals = "first"
  )

  entrez_ids <- na.omit(entrez_ids)  # Remove missing IDs
  message("Number of mapped genes: ", length(entrez_ids))

  return(entrez_ids)
}

# Function to prepare data for topGO
prepare_topgo_data <- function(entrez_gene_list, all_entrez_genes, ontology_type = "BP") {
  # Initialize all genes as 0 (background)
  gene_list <- rep(0, length(all_entrez_genes))
  names(gene_list) <- all_entrez_genes

  # Mark significant genes as 1
  gene_list[entrez_gene_list] <- 1

  # Convert to factor (required by topGO)
  gene_list <- as.factor(gene_list)

  # Check if both levels (0 and 1) exist
  if (length(unique(gene_list)) < 2) {
    stop("Error: allGenes must contain both 1 (significant) and 0 (background). Check input.")
  }

  # Create the `topGOdata` object
  go_data <- new(
    "topGOdata",
    ontology = ontology_type,
    allGenes = gene_list,
    nodeSize = 10,  # Minimum number of genes in a GO term
    annot = annFUN.org,
    mapping = "org.Hs.eg.db",
    ID = "entrez"
  )

  return(go_data)
}

perform_go_analysis <- function(go_data) {
  message("Running GO enrichment using weight01 algorithm...")
  
  # Run GO enrichment test
  result_weight01 <- runTest(go_data, algorithm = "weight01", statistic = "fisher")

  # Get significant GO terms
  go_results <- GenTable(go_data, weight01 = result_weight01, topNodes = 10, numChar = 1000)
  
  # Extract genes mapped to each GO term
  go_results <- go_results %>%
    mutate(MappedGenes = map(GO.ID, ~ {
      entrez_ids <- genesInTerm(go_data, .x)
      
      if (length(entrez_ids) > 0 && !is.null(entrez_ids[[1]])) {
        # Convert Entrez IDs to Gene Symbols (get all mappings)
        gene_symbols <- mapIds(
          org.Hs.eg.db, 
          keys = entrez_ids[[1]], 
          keytype = "ENTREZID", 
          column = "SYMBOL", 
          multiVals = "list"  # Return all possible gene symbols
        )
        
        # Flatten list, remove NAs, and return as comma-separated string
        gene_symbols <- na.omit(unlist(gene_symbols))
        if (length(gene_symbols) > 0) {
          return(paste(gene_symbols, collapse = ", "))
        }
      }
      
      return(NA_character_)  # Return NA if no genes found
    })) 

  return(go_results)
}

# Pipeline function to run the entire GO analysis
run_go_analysis_pipeline <- function(significant_genes, background_genes) {

  if (length(significant_genes) == 0) {
    stop("No significant genes found for GO analysis.")
  }

  # Convert to Entrez IDs
  entrez_gene_list <- convert_to_entrez(significant_genes)

  # Convert background gene list to Entrez IDs (for topGO)
  all_entrez_genes <- convert_to_entrez(background_genes)

  # Prepare `topGOdata` object
  go_data <- prepare_topgo_data(entrez_gene_list, all_entrez_genes, ontology_type = "BP")

  # Perform GO enrichment
  go_results <- perform_go_analysis(go_data)

  return(go_results)
}


#  ============================
#  =         CODES FOR        =
#  =      oposSOM ANALYSIS    =
#  ============================

run_opposom_layer <- function(data_obj, layer_name = 'normalised', full_opossom = FALSE, subtype_df = NULL) {
  flog.info(glue("Running oposSOM on layer: {layer_name}"))

  run_opposom_on_matrix <- function(bulk_layer, pseudo_layer, name, full_opossom = FALSE, subtype_df = NULL) {
    if (!identical(rownames(bulk_layer), rownames(pseudo_layer))) {
      flog.warn(glue("Skipping {name} due to mismatched gene names."))
      return(NULL)
    }

    common_genes <- intersect(rownames(bulk_layer), rownames(pseudo_layer))
    if (length(common_genes) == 0) {
      flog.warn(glue("No common genes found in {name}"))
      return(NULL)
    }

    bulk_layer <- bulk_layer[common_genes, , drop = FALSE]
    pseudo_layer <- pseudo_layer[common_genes, , drop = FALSE]
    # Explicitly suffix with _Bulk or _PseudoBulk to avoid ambiguity
    bulk_ids <- paste0(colnames(bulk_layer), "_Bulk")
    pseudo_ids <- paste0(colnames(pseudo_layer), "_PseudoBulk")

    colnames(bulk_layer) <- bulk_ids
    colnames(pseudo_layer) <- pseudo_ids

    # Compute subtype-wise averages for Bulk and PseudoBulk
    avg_bulk_df <- NULL
    avg_pseudo_df <- NULL
    avg_labels <- c()

    if (!is.null(subtype_df)) {
      # Extract clean sample names from column names
      base_ids <- sub("_(Bulk|PseudoBulk)$", "", colnames(bulk_layer))
      
      subtype_map <- subtype_df$Subtype
      names(subtype_map) <- subtype_df$Sample
      
      # Bulk averages
      for (subtype in unique(subtype_df$Subtype)) {
        samples <- subtype_df$Sample[subtype_df$Subtype == subtype]
        matched_cols <- intersect(paste0(samples, "_Bulk"), colnames(bulk_layer))
        if (length(matched_cols) > 0) {
          avg_expr <- rowMeans(bulk_layer[, matched_cols, drop = FALSE])
          avg_bulk_df <- cbind(avg_bulk_df, avg_expr)
          colnames(avg_bulk_df)[ncol(avg_bulk_df)] <- paste0(subtype, "_Bulk_Avg")
          avg_labels <- c(avg_labels, paste0(subtype, " (Bulk)"))
        }
      }

      # PseudoBulk averages
      for (subtype in unique(subtype_df$Subtype)) {
        samples <- subtype_df$Sample[subtype_df$Subtype == subtype]
        matched_cols <- intersect(paste0(samples, "_PseudoBulk"), colnames(pseudo_layer))
        if (length(matched_cols) > 0) {
          avg_expr <- rowMeans(pseudo_layer[, matched_cols, drop = FALSE])
          avg_pseudo_df <- cbind(avg_pseudo_df, avg_expr)
          colnames(avg_pseudo_df)[ncol(avg_pseudo_df)] <- paste0(subtype, "_PseudoBulk_Avg")
          avg_labels <- c(avg_labels, paste0(subtype, " (PseudoBulk)"))
        }
      }
    }

    # Merge all
    combined <- cbind(bulk_layer, pseudo_layer, avg_bulk_df, avg_pseudo_df)
    # Build metadata for group labels: Subtype (Modality)
    sample_ids <- sub("_(Bulk|PseudoBulk)$", "", colnames(combined))
    modality <- ifelse(grepl("PseudoBulk", colnames(combined)), "PseudoBulk", "Bulk")

    # Match subtype if available
    # Assign subtype from subtype_df
    subtype <- rep("Unknown", length(sample_ids))
    if (!is.null(subtype_df)) {
      match_idx <- match(sample_ids, subtype_df$Sample)
      subtype[!is.na(match_idx)] <- subtype_df$Subtype[match_idx[!is.na(match_idx)]]
    }

    # Override subtype for average columns by extracting directly
    is_avg <- grepl("_Avg$", colnames(combined))
    subtype[is_avg] <- gsub("_(Bulk|PseudoBulk)_Avg$", "", colnames(combined)[is_avg])
    # Final group label: Subtype (Modality)
    group_labels <- paste0(subtype, " (", modality, ")")
    names(group_labels) <- colnames(combined)

    # Print mapping table
    cat("Sample → Group Label Mapping Preview:\n")
    print(data.frame(Sample = names(group_labels), Group = group_labels, check.names = FALSE))
    
    env <- opossom.new(list(dataset.name = paste0("opossom_", name),
                            dim.1stLvlSom = 40))
    if (!full_opossom) {
        env$preferences$activated.modules$geneset.analysis <- FALSE  
        env$preferences$activated.modules$sample.similarity.analysis <- FALSE  
        env$preferences$activated.modules$psf.analysis <- FALSE  
        env$preferences$activated.modules$group.analysis <- FALSE  
        env$preferences$activated.modules$difference.analysis <- FALSE  
    }
    env$preferences$activated.modules$reporting <- FALSE
    env$indata <- as.matrix(combined)
    env$group.labels <- group_labels
    env$preferences$database.dataset <- "hsapiens_gene_ensembl"
    env$preferences$database.id.type <- "hgnc_symbol"

    opossom.run(env)

    output_file <- paste0("opossom_", name, ".RData")
    if (file.exists(output_file)) file.remove(output_file)

    return(env)
  }

  get_layer_by_name <- function(layer_name, bulk, pseudo) {
    parts <- unlist(strsplit(layer_name, " - "))
    for (p in parts) {
      if (!is.null(bulk[[p]]) && !is.null(pseudo[[p]])) {
        bulk <- bulk[[p]]
        pseudo <- pseudo[[p]]
      } else {
        flog.error(glue("Layer '{layer_name}' not found in data_obj."))
        return(NULL)
      }
    }
    return(list(bulk = bulk, pseudo = pseudo))
  }

  if (tolower(layer_name) %in% c("raw", "raw_counts")) {
    flog.info(glue("Skipping {layer_name} because it's raw counts."))
    return(NULL)
  }

  layer_data <- get_layer_by_name(layer_name, data_obj$Bulk, data_obj$PseudoBulk)
  if (is.null(layer_data)) return(NULL)

  bulk_layer <- layer_data$bulk
  pseudo_layer <- layer_data$pseudo

  if (is.matrix(bulk_layer) || is.data.frame(bulk_layer)) {
    return(run_opposom_on_matrix(bulk_layer, pseudo_layer, layer_name, full_opossom, subtype_df = subtype_df))
  } else {
    flog.warn(glue("Layer '{layer_name}' is not a matrix/data.frame."))
    return(NULL)
  }
}



#  ============================
#  =        CODES FOR         =
#  =        PLOTTING          =
#  ============================
plot_multi_experiment_go_split <- function(go_results_named_list, top_n = 5, de_method_param = c("deseq2_res", "limma_res", "edger_res", NULL), force_top_n = FALSE, wrap_bulk = 45, wrap_pb = 45) {
  de_method <- if (is.null(de_method_param)) NULL else match.arg(de_method_param)

  process_go_df <- function(df, direction_label, experiment_name, wrap_width) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      mutate(
        Direction = direction_label,
        GeneRatio = Significant / Annotated,
        Significance = suppressWarnings(as.numeric(gsub("< ", "", weight01))),
        Significance = ifelse(is.na(Significance), 1e-30, Significance),
        NegLogP = -log10(Significance),
        Term = stringr::str_wrap(Term, wrap_width),
        Experiment = experiment_name
      )
  }

  all_data <- purrr::imap_dfr(go_results_named_list, function(go_obj, exp_name) {
    go_layer <- if (!is.null(de_method)) go_obj$go[[de_method]] else go_obj
    bind_rows(
      process_go_df(go_layer$bulk_up, "Prominent in bulk", exp_name, wrap_bulk),
      process_go_df(go_layer$pseudo_up, "Prominent in PB", exp_name, wrap_pb)
    )
  })

  use_overlap_mode <- !force_top_n
  if (use_overlap_mode) {
    overlapping_terms <- all_data %>%
      group_by(Direction, Term) %>%
      summarise(Freq = n_distinct(Experiment), .groups = "drop") %>%
      filter(Freq >= 2) %>%
      group_by(Direction) %>%
      slice_max(order_by = Freq, n = top_n, with_ties = FALSE)

    if (nrow(overlapping_terms) < 2) {
      message("🔍 Not enough overlapping GO terms. Falling back to top ", top_n, " per experiment.")
      use_overlap_mode <- FALSE
    }
  }

  if (!use_overlap_mode) {
    plot_df <- all_data %>%
      group_by(Direction, Experiment) %>%
      slice_max(order_by = NegLogP, n = top_n, with_ties = FALSE)
  } else {
    plot_df <- all_data %>%
      semi_join(overlapping_terms, by = c("Direction", "Term"))
  }

  plot_df <- plot_df %>%
    mutate(
      Term = factor(Term),
      Experiment = factor(
        stringr::str_replace_all(Experiment, "_", " "),
        levels = stringr::str_replace_all(names(go_results_named_list), "_", " ")
      )
    )

  # Get global ranges to share scales
  range_color <- range(plot_df$NegLogP, na.rm = TRUE)
  range_size <- range(plot_df$GeneRatio, na.rm = TRUE)

  shared_colour_scale <- scale_color_gradient(low = "blue", high = "red", limits = range_color)
  shared_size_scale <- scale_size(range = c(3, 10), limits = range_size)

  make_go_plot <- function(direction_label) {
    ggplot(plot_df %>% filter(Direction == direction_label),
           aes(x = Experiment, y = Term, size = GeneRatio, color = NegLogP)) +
      geom_point(alpha = 0.9) +
      shared_colour_scale +
      shared_size_scale +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18),
        legend.position = "right"
      ) +
      labs(
        title = direction_label,
        x = "Experiment",
        y = "GO Term",
        color = "-log10(p-value)",
        size = "Gene Ratio"
      ) +
      guides(
        color = guide_colorbar(title = "-log10(p-value)"),
        size = guide_legend(title = "Gene Ratio")
      )
  }

  list(
    bulk_up = make_go_plot("Prominent in bulk"),
    pseudo_up = make_go_plot("Prominent in PB")
  )
}

plot_single_experiment_go_split <- function(go_obj, 
                                            top_n = 10, 
                                            de_method_param = c("deseq2_res", "limma_res", "edger_res", NULL), 
                                            plt_title = NULL,
                                            string_wrap = 40) {
  de_method <- if (is.null(de_method_param)) NULL else match.arg(de_method_param)
  plt_title <- if (is.null(plt_title)) "" else plt_title
  go_layer <- if (!is.null(de_method)) go_obj$go[[de_method]] else go_obj

  # Helper to process each GO term set
  process_go_df <- function(df, modality_label) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      mutate(
        Modality = modality_label,
        GeneRatio = Significant / Annotated,
        Significance = suppressWarnings(as.numeric(gsub("< ", "", weight01))),
        Significance = ifelse(is.na(Significance), 1e-30, Significance),
        NegLogP = -log10(Significance),
        Term = stringr::str_wrap(Term, string_wrap)
      )
  }

  # Combine both modalities
  plot_df <- bind_rows(
    process_go_df(go_layer$bulk_up, "Bulk"),
    process_go_df(go_layer$pseudo_up, "PB")
  ) %>%
    group_by(Modality) %>%
    slice_max(order_by = NegLogP, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Term = factor(Term, levels = rev(unique(Term))))

  # Global scale ranges
  range_color <- range(plot_df$NegLogP, na.rm = TRUE)
  range_size <- range(plot_df$GeneRatio, na.rm = TRUE)

  # Plot
  ggplot(plot_df, aes(x = Modality, y = Term, size = GeneRatio, color = NegLogP)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "blue", high = "red", limits = range_color) +
    scale_size(range = c(3, 10), limits = range_size) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    ) +
    labs(
      title = plt_title,
      x = "Modality",
      y = "GO Term",
      color = "-log10(p-value)",
      size = "Gene Ratio"
    )
}



plot_pca_results <- function(pca_collection, label, title = NULL, colour_by = NULL, label_points = FALSE, sample_metadata = NULL) {
    
  find_pca_result_by_label <- function(pca_obj, label) {
    if (is.list(pca_obj)) {
      # Base case: current list has a `pca_data$data_type` matching label
      if (!is.null(pca_obj$pca_data) && any(pca_obj$pca_data$data_type == label)) {
        return(pca_obj)
      }

      # Recursive case: check sub-elements
      for (sub in pca_obj) {
        result <- find_pca_result_by_label(sub, label)
        if (!is.null(result)) return(result)
      }
    }

    return(NULL)
  }

  pca_result <- find_pca_result_by_label(pca_collection, label)

  if (is.null(pca_result)) {
    stop(glue("PCA result for label '{label}' not found in the object."))
  }

  pca_df <- pca_result$pca_data
  pca_obj <- pca_result$pca_obj

  if (!is.null(sample_metadata)) {
    pca_df <- left_join(pca_df, sample_metadata, by = c("Sample" = "Sample"))
  }

  var_explained <- (pca_obj$sdev)^2 / sum(pca_obj$sdev^2)
  pc1_var <- round(var_explained[1] * 100, 2)
  pc2_var <- round(var_explained[2] * 100, 2)

  plot_title <- if (is.null(title)) "" else title

  # Aesthetic mapping
  base_aes <- aes(x = PC1, y = PC2, shape = Type)
  if (!is.null(colour_by) && colour_by %in% colnames(pca_df)) {
    base_aes <- modifyList(base_aes, aes(colour = .data[[colour_by]]))
  }

   p <- ggplot(pca_df, base_aes) +
    geom_point(size = 6, alpha = 0.85, na.rm = TRUE) +
    labs(
      title = plot_title,
      x = paste0("PC1 (", pc1_var, "% variance)"),
      y = paste0("PC2 (", pc2_var, "% variance)"),
      colour = colour_by
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.box = "vertical",         # Stack legends vertically
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    ) +
    guides(
      colour = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)
    )

  if (label_points) {
    p <- p + ggrepel::geom_text_repel(aes(label = Sample), size = 7, max.overlaps = 20)
  }

  return(p)
}


plot_volcano <- function(de_df,
                         logfc_col = "log2FoldChange",
                         padj_col = "padj",
                         logfc_thresh = 1,
                         padj_thresh = 0.05,
                         highlight_col = NULL,
                         title =  NULL) {
  if (is.null(de_df)) return(NULL)
  de_df <- as.data.frame(de_df)

  if (!all(c(logfc_col, padj_col) %in% colnames(de_df))) {
    stop("Missing required columns in DE result dataframe.")
  }

  title <- if (is.null(title)) "" else title

  de_df <- de_df %>%
    mutate(
      Significance = case_when(
        !!sym(padj_col) < padj_thresh & abs(!!sym(logfc_col)) > logfc_thresh ~ "Significant",
        TRUE ~ "Not significant"
      ),
      Highlight = case_when(
        !is.null(highlight_col) & .data[[highlight_col]] == "Significant" & Significance == "Significant" ~ "Informative DEGs",
        Significance == "Significant" ~ "Initial DEGs",
        TRUE ~ "Not DEGs"
      ),
      log10_padj = -log10(!!sym(padj_col))
    )

  # Count stats
  sig_count <- sum(de_df$Significance == "Significant", na.rm = TRUE)
  spca_count <- sum(de_df$Highlight == "Informative DEGs", na.rm = TRUE)
  all_gene_count <- nrow(de_df)
  subtitle_text <- paste("Total Gene Count:", all_gene_count ,"\nInitial DEGs:", sig_count, "\nInformative DEGs:", spca_count)

  ggplot(de_df, aes(x = !!sym(logfc_col), y = log10_padj)) +
    geom_point(aes(color = Highlight), size = 3, alpha = 0.85, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Informative DEGs" = "royalblue",
        "Initial DEGs" = "coral",
        "Not DEGs" = "grey80"
      )
    ) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = "Log2 Fold Change",
      y = expression(-log[10]("Adjusted P-value")),
      color = "Gene Classication"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    )
}

plot_som_portrait_by_pos <- function(env,
                                     save_path = "pos_SOM_Portraits.png",
                                     width = 16, height = 10, dpi = 300,
                                     sample_pos = c(1, 2),
                                     font_size = 18,
                                     sub_plt_title = FALSE,
                                     force_rename = FALSE,
                                     force_rename_suffix = NULL) {
  message("Generating SOM portraits by sample_pos indexing for illustration...")

  all_dfs <- list()
  sample_names <- names(env$group.labels)
  group_labels <- env$group.labels

  parsed <- stringr::str_match(group_labels, "^(.*) \\((.*)\\)$")
  subtype <- parsed[, 2]
  modality <- parsed[, 3]

  metadata_matrix <- env$metadata
  som_dim <- env$preferences$dim.1stLvlSom
  palette <- env$color.palette.portraits(1000)

  for (i in sample_pos) {
    sample <- sample_names[i]
    expr_vector <- metadata_matrix[, sample]
    expression_matrix <- matrix(expr_vector, som_dim, som_dim)

    df <- expand.grid(x = 1:som_dim, y = 1:som_dim)
    df$expression <- as.vector(expression_matrix)
    df$Subtype <- subtype[i]
    df$Modality <- modality[i]

    clean_sample <- gsub("(_Bulk|_PseudoBulk|^Average_)", "", sample)
    is_avg <- grepl("Avg|Average", sample, ignore.case = TRUE)

    df$SampleCore <- clean_sample
    df$is_avg <- is_avg
    all_dfs[[sample]] <- df
  }

  full_df <- bind_rows(all_dfs)
  full_df$Modality <- dplyr::recode(full_df$Modality, PseudoBulk = "PB")
  full_df$Modality <- factor(full_df$Modality, levels = c("Bulk", "PB"))
  full_df$Subtype <- factor(full_df$Subtype, levels = sort(unique(full_df$Subtype)))

  # Create plots by subtype
  plot_list <- list()

  for (sub in levels(full_df$Subtype)) {
    df_sub <- full_df %>% filter(Subtype == sub)

    # Assign labels
    averaged_label <- paste0(sub, "\nAveraged")
    patient_samples <- df_sub %>% filter(!is_avg) %>% distinct(SampleCore)

    if (force_rename) {
      rename_suffix <- if (is.null(force_rename_suffix)) "Sample " else force_rename_suffix
      patient_samples <- patient_samples %>% mutate(Rename = paste0(rename_suffix, row_number()))
    } else {
      patient_samples <- patient_samples %>% mutate(Rename = SampleCore)
    }

    df_sub <- df_sub %>%
      left_join(patient_samples, by = "SampleCore") %>%
      mutate(Sample = ifelse(is_avg, averaged_label, Rename)) %>%
      group_by(Modality) %>%
      arrange(Sample, .by_group = TRUE) %>%
      ungroup() %>%
      mutate(Sample = factor(Sample, levels = c(averaged_label, patient_samples$Rename)))

    p <- ggplot(df_sub, aes(x = x, y = y, fill = expression)) +
      geom_tile() +
      facet_grid(rows = vars(Modality), cols = vars(Sample),
                 scales = "free", space = "free", switch = "y") +
      scale_fill_gradientn(colors = palette) +
      coord_cartesian() +
      theme_minimal(base_size = font_size) +
      theme(
        strip.text.x = element_text(size = font_size, face = 'bold'),
        strip.text.y.left = element_text(size = font_size, face = 'bold'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none",
        panel.spacing = unit(0.2, "lines"),
        strip.placement = "outside"
      ) +
      labs(x = NULL, y = NULL)

    plot_list[[sub]] <- p
  }

  final_plot <- wrap_plots(plot_list, ncol = 1, guides = "collect") &
    theme(legend.position = "none")

  ggsave(save_path, plot = final_plot, width = width, height = height, dpi = dpi)
  message("Saved to ", save_path)
  final_plot <- final_plot + plot_layout(widths = 1) & theme(plot.margin = margin())
  return(final_plot)
}


plot_dynamic_subtype_som_portraits <- function(env,
                                               save_path = "Subtype_SOM_Portraits.png",
                                               width = 16, height = 10, dpi = 300,
                                               font_size = 18, sub_plt_title = FALSE,
                                               force_rename = FALSE,
                                               force_rename_suffix = NULL,
                                               make_title_som_sample_name = F,
                                               som_sample_name_wrap = 30,
                                               force_som_title = NULL) {
  message("Generating SOM portraits with patchwork layout by subtype...")

  all_dfs <- list()
  sample_names <- names(env$group.labels)
  group_labels <- env$group.labels

  parsed <- stringr::str_match(group_labels, "^(.*) \\((.*)\\)$")
  subtype <- parsed[, 2]
  modality <- parsed[, 3]

  metadata_matrix <- env$metadata
  som_dim <- env$preferences$dim.1stLvlSom
  palette <- env$color.palette.portraits(1000)

  for (i in seq_along(sample_names)) {
    sample <- sample_names[i]
    expr_vector <- metadata_matrix[, sample]
    expression_matrix <- matrix(expr_vector, som_dim, som_dim)

    df <- expand.grid(x = 1:som_dim, y = 1:som_dim)
    df$expression <- as.vector(expression_matrix)
    df$Subtype <- subtype[i]
    df$Modality <- modality[i]

    clean_sample <- gsub("(_Bulk|_PseudoBulk|^Average_)", "", sample)
    is_avg <- grepl("Avg|Average", sample, ignore.case = TRUE)

    df$SampleCore <- clean_sample
    df$is_avg <- is_avg
    all_dfs[[sample]] <- df
  }

  full_df <- bind_rows(all_dfs)
  full_df$Modality <- dplyr::recode(full_df$Modality, PseudoBulk = "PB")
  full_df$Modality <- factor(full_df$Modality, levels = c("Bulk", "PB"))
  full_df$Subtype <- factor(full_df$Subtype, levels = sort(unique(full_df$Subtype)))

  plot_list <- list()

  for (sub in levels(full_df$Subtype)) {
    df_sub <- full_df %>% filter(Subtype == sub)

    title_txt <- if (!is.null(force_som_title)) force_som_title else sub
    # Identify and label averaged and patient samples
    if (make_title_som_sample_name) {
      averaged_label <- paste0(stringr::str_wrap(stringr::str_to_title(title_txt), som_sample_name_wrap), '\nMean')
    } else {
      averaged_label <- paste0(stringr::str_wrap(title_txt, som_sample_name_wrap), '\nMean')
    }


    patient_samples <- df_sub %>% filter(!is_avg) %>% distinct(SampleCore)

    if (force_rename) {
      rename_suffix <- if (is.null(force_rename_suffix)) "Sample " else force_rename_suffix
      patient_samples <- patient_samples %>% mutate(Rename = paste0(rename_suffix, row_number()))
    } else {
      patient_samples <- patient_samples %>% mutate(Rename = SampleCore)
    }

    df_sub <- df_sub %>% 
      left_join(patient_samples, by = "SampleCore") %>%
      mutate(Sample = ifelse(is_avg, averaged_label, Rename)) %>%
      group_by(Modality) %>%
      arrange(Sample, .by_group = TRUE) %>%
      ungroup() %>%
      mutate(Sample = factor(Sample, levels = c(averaged_label, patient_samples$Rename)))


    p <- ggplot(df_sub, aes(x = x, y = y, fill = expression)) +
      geom_tile() +
      facet_grid(rows = vars(Modality), cols = vars(Sample),
                 scales = "free", space = "free", switch = "y") +
      scale_fill_gradientn(colors = palette) +
      coord_cartesian() +
      theme_minimal(base_size = font_size) +
      theme(
        strip.text.x = element_text(size = font_size, face = 'bold'),
        strip.text.y.left = element_text(size = font_size, face = 'bold'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none",
        panel.spacing = unit(0.2, "lines"),
        strip.placement = "outside"
      ) +
      labs(x = NULL, y = NULL)

    if (!sub_plt_title) sub <- ""

    title_plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = sub,
               size = font_size / 3, fontface = "bold", hjust = 0.5, vjust = 0.5)

    patch <- title_plot / p + plot_layout(heights = c(0.12, 1))
    plot_list[[sub]] <- patch
  }

  final_plot <- wrap_plots(plot_list, ncol = 1, guides = "collect") &
    theme(legend.position = "none")

  ggsave(save_path, plot = final_plot, width = width, height = height, dpi = dpi)
  message("Saved to ", save_path)
  final_plot <- final_plot + plot_layout(widths = 1) & theme(plot.margin = margin())
  return(final_plot)
}

plot_kde_results <- function(kde_df, label = NULL, title = NULL, custom_genes = NULL) {
  plot_data <- if (!is.null(label)) {
    kde_df %>% filter(Layer == label)
  } else {
    kde_df
  }
  if (!is.null(custom_genes) && "Gene" %in% colnames(plot_data)) {
    plot_data <- plot_data %>% filter(Gene %in% custom_genes)
  }

  plot_title <- if (!is.null(title)) title else ""

  ggplot(plot_data, aes(x = Expression, fill = Modality, colour = Modality)) +
    geom_density(alpha = 0.4, adjust = 1.2) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    ) +
    labs(
      title = plot_title,
      x = "Expression Value",
      y = "Density",
      fill = "Modality",
      colour = "Modality"
    )
}

plot_pca_variance_explained <- function(pca_result, top_n = NULL, title =NULL ) {
  # Check for sdev or eigenvalues
  if (!"sdev" %in% names(pca_result)) {
    stop("The input must be a PCA object with a 'sdev' component (e.g., from prcomp() or similar).")
  }
  title <- if (is.null(title)) "" else title
  # Compute variance explained
  var_explained <- pca_result$sdev^2
  var_explained_ratio <- var_explained / sum(var_explained)

  # Default to all PCs if top_n is not specified
  if (is.null(top_n)) top_n <- length(var_explained_ratio)

  # Build data frame
  df <- data.frame(
    PC = factor(paste0("PC", seq_along(var_explained_ratio)), levels = paste0("PC", seq_along(var_explained_ratio)))
    # Variance Explained = var_explained_ratio,
    # Cumulative Explained = cumsum(var_explained_ratio)
  )
  df['Variance Explained'] <- var_explained_ratio
  df['Cumulative Explained'] <- cumsum(var_explained_ratio)
  df_top <- df[1:top_n, ]


  df_long <- df_top %>%
    pivot_longer(cols = c("Variance Explained", "Cumulative Explained"),
                 names_to = "Metric", values_to = "Value")

  
  ggplot(df_long, aes(x = PC, y = Value, group = Metric, fill = Metric, color = Metric)) +
    geom_col(data = subset(df_long, Metric == "Variance Explained"), aes(fill = Metric), width = 0.7) +
    geom_line(data = subset(df_long, Metric == "Cumulative Explained"), aes(color = Metric), size = 1) +
    geom_point(data = subset(df_long, Metric == "Cumulative Explained"), aes(color = Metric), size = 2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("Variance Explained" = "steelblue")) +
    scale_color_manual(values = c("Cumulative Explained" = "red")) +
    labs(
      title = title,
      x = "Principal Components",
      y = "Variance Explained (%)",
      fill = "Metric",
      color = "Metric"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold",size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 18)
    )
}


plot_correlation_boxplot <- function(correlation_data, plt_title = '', save_path = '', width = 14, height = 6, dpi = 300) {
  # Define plot title
  plot_title <- if (nchar(plt_title) > 0) plt_title else ''

  # Rename types for improved readability
  correlation_data$Type <- recode(
    correlation_data$Type,
    "Bulk_vs_Bulk" = "Bulk–Bulk",
    "PB_vs_PB" = "PB–PB",
    "Bulk_vs_PB" = "Bulk–PB"
  )

  # Set factor levels to control plotting order
  correlation_data$Type <- factor(
    correlation_data$Type,
    levels = c("Bulk–Bulk", "PB–PB", "Bulk–PB")
  )

  # Define color mapping
  color_map <- c(
    "Bulk–Bulk" = "#F8766D",
    "PB–PB" = "#00BFC4",
    "Bulk–PB" = "#66c2a5"
  )

  quantiles_df <- correlation_data %>%
    group_by(Type) %>%
    summarise(
      Min = min(Correlation),
      Q1 = quantile(Correlation, 0.25),
      Median = median(Correlation),
      Q3 = quantile(Correlation, 0.75),
      Max = max(Correlation)
    )


  print(quantiles_df)
  # Build plot
  correlation_plot <- ggplot(correlation_data, aes(x = Type, y = Correlation, fill = Type, colour = Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    scale_fill_manual(values = color_map, name = "Correlation Type") +
    scale_colour_manual(values = color_map, name = "Correlation Type") +
    theme_classic(base_size = 18) +
    labs(title = plot_title, x = NULL, y = "Pearson Correlation") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.ticks.x = element_blank(), # Remove x-axis ticks
      legend.position = "top",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    )

  # Save if requested
  if (nchar(save_path) > 0) {
    ggsave(save_path, plot = correlation_plot, width = width, height = height, dpi = dpi)
  }

  return(correlation_plot)
}
