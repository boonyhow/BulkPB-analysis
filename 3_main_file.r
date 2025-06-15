source('code_modules/3_functions.r')
# === Check if precomputed objects exist ===
need_to_save_objects <- TRUE
if (file.exists("3_full_analysis_objects.rds")) {
  message("Found existing 3_full_analysis_objects.rds. Loading into environment...")
  list2env(readRDS("3_full_analysis_objects.rds"), .GlobalEnv)
  need_to_save_objects <- FALSE
}

#  ============================
#  =       Loading all        =
#  =          data            =
#  ============================
## GSE176078
d_176078 <- readRDS('data/GSE176078/data_for_study/intermediate_data/2_datapreprocessing_intersect_0-2.rds')
subtype_d_176078 <- read.csv('data/GSE176078/data_for_study/sample_subtype.csv')

## GSE217517
d_217517 <- readRDS('data/GSE217517/data_for_study/intermediate_data/2_datapreprocessing_intersect_0-2.rds')
subtype_d_217517 <- read.csv('data/GSE217517/data_for_study/sample_subtype.csv')

## GSE226163
d_226163 <- readRDS('data/GSE226163/data_for_study/intermediate_data/2_datapreprocessing_intersect_0-2.rds')
subtype_d_226163 <- read.csv('data/GSE226163/data_for_study/sample_subtype.csv')


#  ============================
#  =       Running all        =
#  =       calculations       =
#  ============================

## KDE
d_176078_kde <- run_all_kde(d_176078, use_all_samples = FALSE)
d_217517_kde <- run_all_kde(d_217517, use_all_samples = FALSE)
d_226163_kde <- run_all_kde(d_226163, use_all_samples = FALSE)


d_176078_sparse_pca <- run_all_sparse_pca(d_176078, subtype_d_176078)
d_217517_sparse_pca <- run_all_sparse_pca(d_217517, subtype_d_217517) 
d_226163_sparse_pca <- run_all_sparse_pca(d_226163, subtype_d_226163)


d_176078_opossom <- run_opposom_layer(d_176078, subtype_df = subtype_d_176078)
d_217517_opossom <- run_opposom_layer(d_217517, subtype_df = subtype_d_217517)
d_226163_opossom <- run_opposom_layer(d_226163, subtype_df = subtype_d_226163)


d_176078_de_go <- run_de_pipelines(d_176078, subtype_d_176078, run_go = T)
d_217517_de_go <- run_de_pipelines(d_217517, subtype_d_217517, run_go = T)
d_226163_de_go <- run_de_pipelines(d_226163, subtype_d_226163, run_go = T)


d_176078_corr <- run_all_correlation(d_176078, samplewise_correlation = F)
d_217517_corr <- run_all_correlation(d_217517, samplewise_correlation = F)
d_226163_corr <- run_all_correlation(d_226163, samplewise_correlation = F)

de_go_results_list <- list(
  GSE217517 = d_217517_de_go,
  GSE176078 = d_176078_de_go,
  GSE226163 = d_226163_de_go
)

sparse_pca_results_list <- list(
  GSE217517 = d_217517_sparse_pca,
  GSE176078 = d_176078_sparse_pca,
  GSE226163 = d_226163_sparse_pca
)

deseq_de_comparison <- summarise_spca_vs_de_overlap(de_go_results_list, sparse_pca_results_list, de_method = 'deseq2_res')

d_176078_sparse_pca_go <-  run_go_by_direction(deseq_de_comparison$updated_de$GSE176078, sig_col = "spca_de_sig", logFC_threshold = 1)
d_217517_sparse_pca_go <-  run_go_by_direction(deseq_de_comparison$updated_de$GSE217517, sig_col = "spca_de_sig", logFC_threshold = 1)
d_226163_sparse_pca_go <-  run_go_by_direction(deseq_de_comparison$updated_de$GSE226163, sig_col = "spca_de_sig", logFC_threshold = 1)

sparse_pca_go_results_list <- list(
  OV = d_217517_sparse_pca_go,
  BC = d_176078_sparse_pca_go
  # iPSC_derived_cells = d_226163_sparse_pca_go
)

#  ============================
#  =       Plotting all       =
#  =          graphs          =
#  ============================

## CREATING PLOT OBJECTS 

p_corr_217517 <- plot_correlation_boxplot(d_217517_corr$normalised, plt_title = 'OV')
p_corr_176078 <- plot_correlation_boxplot(d_176078_corr$normalised, plt_title = 'BC')
p_corr_226163 <- plot_correlation_boxplot(d_226163_corr$normalised)

p_kde_217517 <-   plot_kde_results(d_217517_kde, label = "normalised", title = 'OV') 
p_kde_176078 <-   plot_kde_results(d_176078_kde, label = "normalised", title = 'BC') 
p_kde_226163 <-   plot_kde_results(d_226163_kde, label = "normalised") 

p_spca_217517 <-   plot_pca_results(d_217517_sparse_pca, colour_by = "Sample", label_points = F, 'normalised', title = 'OV') 
p_spca_176078 <-   plot_pca_results(d_176078_sparse_pca, colour_by = "Subtype", label_points = F, 'normalised', title = 'BC') 
p_spca_226163 <-   plot_pca_results(d_226163_sparse_pca, colour_by = "Subtype", label_points = F, 'normalised') 
      
p_217517_opossom <- plot_dynamic_subtype_som_portraits(d_217517_opossom, width = 30, height = 10, force_som_title = 'OV')
p_176078_opossom <- plot_dynamic_subtype_som_portraits(d_176078_opossom, save_path = 'supp_plot_opossom_176078.png', width = 30 / 2, height = 30 / 2, force_rename =  T, force_rename_suffix = "P") # 3 Subtypes, hence height is 30
p_226163_opossom <- plot_dynamic_subtype_som_portraits(d_226163_opossom, save_path = 'supp_plot_opossom_226163.png', 
                                                          width = 30 / 2, 
                                                          height = 20/ 2, 
                                                          force_rename =  T, 
                                                          force_rename_suffix = "Replicate ", 
                                                          make_title_som_sample_name = T,
                                                          som_sample_name_wrap = 12
                                                        ) # 2 Subtypes, hence height is 20

p_176078_sparse_pca_go <- wrap_elements(full = plot_single_experiment_go_split(d_176078_sparse_pca_go, top_n= 5, de_method_param =  NULL)) 
p_217517_sparse_pca_go <- wrap_elements(full = plot_single_experiment_go_split(d_217517_sparse_pca_go, top_n= 5, de_method_param =  NULL))
p_226163_sparse_pca_go <- wrap_elements(full = plot_single_experiment_go_split(d_226163_sparse_pca_go, top_n= 5, de_method_param =  NULL, string_wrap = 30))

p_217517_volcano <- plot_volcano(deseq_de_comparison$updated_de$GSE217517, highlight_col = "spca_de_sig", title = 'OV') 
p_226163_volcano <- plot_volcano(deseq_de_comparison$updated_de$GSE226163, highlight_col = "spca_de_sig") 
p_176078_volcano <- plot_volcano(deseq_de_comparison$updated_de$GSE176078, highlight_col = "spca_de_sig", title = 'BC') 

spare_pca_data_go <- plot_multi_experiment_go_split(sparse_pca_go_results_list, top_n = 5,  de_method_param = NULL,
  force_top_n = TRUE, wrap_bulk = 40, wrap_pb = 30)

p_spca_217517 <- p_spca_217517 + guides(color = "none", shape = "none")
p_spca_176078 <- p_spca_176078 + guides(color = "none")
# p_spca_226163 <- p_spca_226163 + guides(color = "none", shape = "none")

## GENERATING LIST OF PLOT OBJECTS 

overall_plots_no_ipsc <- list(
  sparse_pca = (
    (p_spca_217517 + p_spca_176078 ) + 
    plot_layout(guides = "collect") +
    plot_annotation() & 
    theme(legend.position = "top")
  )
  ,
  correlation_boxplot = (
    (p_corr_217517 + p_corr_176078 )+ 
    plot_layout(guides = "collect") +
    plot_annotation() & 
    theme(legend.position = "top")
  ),
  kde_density = (
      (p_kde_217517 + p_kde_176078 )+ 
      plot_layout(guides = "collect") +
      plot_annotation() & 
      theme(legend.position = "top")
    ),
  volcano_plot = (
      (p_217517_volcano + p_176078_volcano )+ 
      plot_layout(guides = "collect") +
      plot_annotation() & 
      theme(legend.position = "top")
    ),
  go_plot =((spare_pca_data_go$bulk | spare_pca_data_go$pseudo)+ 
      plot_layout(guides = "collect") +
      plot_annotation() & 
      theme(legend.position = "right")) , 
  som_portraits = p_217517_opossom
)


scale_factor <- 2.5
save_named_plot_list(
  overall_plots_no_ipsc,
  folder_path = "plots/no_ipsc_overall_large",
  width = 35 /scale_factor, height = 15 / scale_factor, dpi = 300,
  format = 'pdf',
  custom_sizes = list(
    som_portraits = c(35 / scale_factor, 10 / scale_factor),  # 3× default width
    go_plot = c(35 / scale_factor, 20 / scale_factor)        # taller PCA
  )
)




#  ============================
#  =       Supplementary      =
#  =         materials        =
#  ============================

#### EXPLAINED VARIANCE PCA PLOT

p_176078_sparse_pca_variance_explained <- plot_pca_variance_explained(d_176078_sparse_pca$normalised$pca_obj, title = 'BC')
p_217517_sparse_pca_variance_explained <- plot_pca_variance_explained(d_217517_sparse_pca$normalised$pca_obj, title = 'OV')
p_226163_sparse_pca_variance_explained <- plot_pca_variance_explained(d_226163_sparse_pca$normalised$pca_obj, title = 'iPSC')


variance_explained_plots <- (
  (p_217517_sparse_pca_variance_explained + p_176078_sparse_pca_variance_explained + p_226163_sparse_pca_variance_explained) 
) +
  plot_layout(guides = "collect", axis = "collect_y") +  # Align y-axis (left)
  plot_annotation() &
  theme(legend.position = "top")

ggsave('supp_materials_variance_explained.pdf', variance_explained_plots, width = 18, height = 6)

##################################


#### iPSC ONLY PLOTS #############
plots_226163 <- list(
  correlation_boxplot = p_corr_226163,
  kde_density = p_kde_226163,
  sparse_pca = p_spca_226163,
  som_portraits = p_226163_opossom,
  go_plot = p_226163_sparse_pca_go,
  volcano_plot = p_226163_volcano
)

save_named_plot_list(
  plots_226163,
  folder_path = "plots/ipsc",
  width = 25 /scale_factor, height = 15 / scale_factor, dpi = 300,
  format = 'pdf',
  custom_sizes = list(
    som_portraits = c(25 / scale_factor,15 / scale_factor),  # 3× default width
    go_plot = c(25 / scale_factor, 25 / scale_factor)        # taller PCA
  )
)
###################################

## T TEST OF SPCA SEPERATION
d_176078_sparse_pca_test <- test_spca_separation(d_176078_sparse_pca$normalised)
d_217517_sparse_pca_test <- test_spca_separation(d_217517_sparse_pca$normalised)
d_226163_sparse_pca_test <- test_spca_separation(d_226163_sparse_pca$normalised)

## OPOSSOMS
plot_dynamic_subtype_som_portraits(d_226163_opossom, save_path = 'supp_plot_opossom_226163.png',width = 25/scale_factor, height = 20/scale_factor, force_rename =  T, force_rename_suffix = "Replicate ")
plot_dynamic_subtype_som_portraits(d_176078_opossom, save_path = 'supp_plot_opossom_176078.png', width = 25/scale_factor, height = 30/scale_factor, force_rename =  T, force_rename_suffix = "P")
plot_som_portrait_by_pos(d_217517_opossom, save_path = 'supp_1_sample_opossom_217517.png', width = 8, height = 10, sample_pos = c(1,  8), font_size = 28)


#  ============================
#  =         Saving all       =
#  =          objects         =
#  ============================
all_objects <- list(
  d_176078 = d_176078,
  d_217517 = d_217517,
  d_226163 = d_226163,
  subtype_d_176078 = subtype_d_176078,
  subtype_d_217517 = subtype_d_217517,
  subtype_d_226163 = subtype_d_226163,
  
  d_176078_kde = d_176078_kde,
  d_217517_kde = d_217517_kde,
  d_226163_kde = d_226163_kde,
  
  d_176078_sparse_pca = d_176078_sparse_pca,
  d_217517_sparse_pca = d_217517_sparse_pca,
  d_226163_sparse_pca = d_226163_sparse_pca,
  
  d_176078_opossom = d_176078_opossom,
  d_217517_opossom = d_217517_opossom,
  d_226163_opossom = d_226163_opossom,
  
  d_176078_de_go = d_176078_de_go,
  d_217517_de_go = d_217517_de_go,
  d_226163_de_go = d_226163_de_go,
  
  d_176078_sparse_pca_go = d_176078_sparse_pca_go, 
  d_217517_sparse_pca_go = d_217517_sparse_pca_go,
  d_226163_sparse_pca_go = d_226163_sparse_pca_go,

  d_176078_corr = d_176078_corr,
  d_217517_corr = d_217517_corr,
  d_226163_corr = d_226163_corr,

  sparse_pca_results_list = sparse_pca_results_list,
  de_go_results_list = de_go_results_list,
  
  deseq_de_comparison = deseq_de_comparison,
  
  sparse_pca_go_results_list = sparse_pca_go_results_list
)

if (need_to_save_objects) {
  saveRDS(all_objects, file = "3_full_analysis_objects.rds")
  message("Saved 3_full_analysis_objects.rds")
}