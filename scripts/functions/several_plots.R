#"""
#1. The input should be a matrix compatible with UpSet plot.
#2. mode should be optional, more than just intersect.
#3. combination_column, background_column, and backgroud_patch_column should be given.
#4. set_order can be adjusted by a vector of names.
#5. ylim, and height also should be adjustable.
#6. annotation_name_side, rotation angle.
#7. parameters for decoration also should be adjustable.
#"""





library("ComplexHeatmap")
upset_matrix <- data.frame(read.csv("/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/UpSet_matrix1.tsv",
                                    header = TRUE,
                                    sep = "\t"))
rownames(upset_matrix) <- upset_matrix[, colnames(upset_matrix)[[1]]]
upset_matrix[, colnames(upset_matrix)[[1]]] <- NULL
upset_matrix <- upset_matrix[as.logical(rowSums(upset_matrix != 0)),]
comb <- ComplexHeatmap::make_comb_mat(upset_matrix, mode = "intersect")
c_size <- ComplexHeatmap::comb_size(comb)
sets = c("RAI.yes", "X.partners.3", "Oral.yes",  "STI.positive", "condom.no")
upset_plot <- ComplexHeatmap::UpSet(comb,
                                    comb_col = "#fb5238",
                                    bg_col = "#ffbeab",
                                    bg_pt_col = "#ffdfd5",
                                    set_order = sets,
                                    comb_order = order(c_size),
                                    top_annotation = HeatmapAnnotation(
                                        "# shared taxonomic biomarkers" = anno_barplot(c_size,
                                                                                    ylim = c(0, max(c_size)*1.1),
                                                                                    border = FALSE,
                                                                                    gp = gpar(fill = "#fb5238", col = "#fb5238"),
                                                                                    height = unit(8, "cm")),
                                    annotation_name_side = "left",
                                    annotation_name_rot = 90),
                                    right_annotation = NULL
                                    )

upset_plot <- draw(upset_plot)
col_order <- column_order(upset_plot)

decorate_annotation("# shared taxonomic biomarkers", {
  grid.text(c_size[col_order], x = seq_along(c_size), y = unit(c_size[col_order], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
