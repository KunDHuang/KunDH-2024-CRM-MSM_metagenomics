

###### Testing code  ######
mat <- "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/prevotellaceae_matrix_4ComplexHeatmap.tsv"
row_md <- "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/prevotellaceae_matrix_4ComplexHeatmap_species_md.txt"
col_md <- "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/prevotellaceae_matrix_4ComplexHeatmap_sample_md.txt"



# source(file = "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/scripts/functions/beta_diversity_funcs.R")

# coor_df <- generate_coordis_df(mat, md, "euclidean")
# View(coor_df)

# pcoa_plot(mat, md, "bray", "condom_use", 20, 4, to_rm = c("no_receptive_anal_intercourse"))
# est_permanova(mat, md, "condom_use", c("Antibiotics_6mo", "HIV_status", "inflammatory_bowel_disease", "BMI_kg_m2_WHO", "diet"))
# est_permanova(mat = mat, 
#               md = md, 
#               variable = "condom_use", 
#               covariables = c("Antibiotics_6mo", "HIV_status", "inflammatory_bowel_disease", "BMI_kg_m2_WHO", "diet"),
#               nper = 999, 
#               to_rm = c("no_receptive_anal_intercourse"),
#               by_method = "margin")

# pcoa_df <- data.frame(read.csv("/Users/kunhuang/R_analysis_mirror/msm_analysis/manuscript_05012023/figure1/pcoa_df.tsv",
#                       header = TRUE,
#                       sep = "\t"))

# length(unique(pcoa_df[, "sexual_orientation"]))
# pcoa_sideplot(coordinate_df = pcoa_df,
#               coordinate_1 = "PC1",
#               coordinate_2 = "PC2",
#               variable = "sexual_orientation")

source(file = "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/scripts/functions/complexheatmap_plotting_funcs.R")
col_func <- viridis::viridis(100)
svglite::svglite("/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/images/prevotellaceae_relab.svg", width = 14.5, height = 7.5)
plot_complex_heatmap(mat,
                     color_bar_name = "relative abundance (log10)",
                     row_md = row_md,
                     column_md = col_md,
                     show_col_names = FALSE,
                     show_row_names = TRUE,
                     width = 2,
                     height = 3.5,
                     row_names_side = "left",
                     cluster_columns = T,
                     cluster_column_slices = F,
                     cluster_rows = F,
                     cluster_row_slices = F,
                     border = F,
                     row_gap = 1,
                     column_gap = 1,
                     color_func = col_func,
                     transformation = "log10")
dev.off()

global_mat <- "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/global_enrichment_matrix.tsv"
global_row_md <- "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/global_enrichment_matrix_rownames.tsv"
global_col_md <- "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/global_enrichment_matrix_colnames.tsv"

col_func <- circlize::colorRamp2(c(0, 1), hcl_palette = "Blues", reverse = T)
svglite::svglite("/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/images/global_presence.svg", width = 17.5, height = 11.5)
plot_complex_heatmap(global_mat,
                     row_md = global_row_md,
                     column_md = global_col_md,
                     show_col_names = F,
                     show_row_names = TRUE,
                     width = 0.3,
                     height = 3.5,
                     row_names_side = "left",
                     column_names_side = "top", 
                     cluster_columns = F,
                     cluster_column_slices = F,
                     cluster_rows = F,
                     cluster_row_slices = F,
                     border = T,
                     row_gap = 1,
                     column_gap = 1,
                     color_func = col_func,
                     transformation = "binary")
dev.off()


source(file = "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/scripts/functions/rocauc_stdv_funcs.R")
roc_auc_merged <- data.frame(read.csv("/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/roc_auc_merged.tsv",
                                      header = TRUE,
                                      sep = "\t"))


std <- data_summary(roc_auc_merged, "roc.auc", "sexual.practice")

dev_plot <- std_deviation_plot(std, "sexual.practice", "roc.auc", "sd",
                               order = c("Receptive anal intercourse", "Number of partners",
                    "Oral sex", "Sex transmitted infection", "Condom use"))
dev_plot + ggplot2::ylim(0, 1) + ggpubr::rotate_x_text(45)
