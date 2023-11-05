

###### Testing code  ######
mat <- read.csv("/Users/kunhuang/R_analysis_mirror/msm_analysis/sexual_practice_analysis/msm_mpa4_run2_matrix.tsv",
                    header = TRUE,
                    sep = "\t")
md <- read.csv("/Users/kunhuang/R_analysis_mirror/msm_analysis/sexual_practice_analysis/msm_mpa4_run2_metadata.tsv",
                    header = TRUE,
                    sep = "\t")



source(file = "/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/scripts/functions/beta_diversity_funcs.R")

coor_df <- generate_coordis_df(mat, md, "euclidean")
View(coor_df)

pcoa_plot(mat, md, "bray", "condom_use", 20, 4, to_rm = c("no_receptive_anal_intercourse"))
est_permanova(mat, md, "condom_use", c("Antibiotics_6mo", "HIV_status", "inflammatory_bowel_disease", "BMI_kg_m2_WHO", "diet"))
est_permanova(mat = mat, 
              md = md, 
              variable = "condom_use", 
              covariables = c("Antibiotics_6mo", "HIV_status", "inflammatory_bowel_disease", "BMI_kg_m2_WHO", "diet"),
              nper = 999, 
              to_rm = c("no_receptive_anal_intercourse"),
              by_method = "margin")

pcoa_df <- data.frame(read.csv("/Users/kunhuang/R_analysis_mirror/msm_analysis/manuscript_05012023/figure1/pcoa_df.tsv",
                      header = TRUE,
                      sep = "\t"))

length(unique(pcoa_df[, "sexual_orientation"]))
pcoa_sideplot(coordinate_df = pcoa_df,
              coordinate_1 = "PC1",
              coordinate_2 = "PC2",
              variable = "sexual_orientation")
