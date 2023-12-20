# Visualizations for microbiome biomarkers 
This tutorial is to visualize microbiome biomarkers with different plots.

### Deviation plot

#### R packages required

* [ggpubr](https://rpkgs.datanovia.com/ggpubr/)

In this tutorial, we will use function `ggdotchart` from R package `ggpubr` to visualize LefSe [biomarkers](../example_data/npartners_lefse_deviation_plot.tsv) associated with the number of partners in MSM individuals.   

Open a new working R script, and load our example [biomarkers](../example_data/npartners_lefse_deviation_plot.tsv) from `path_to_the_package/KunDH-2023-CRM-MSM_metagenomics/examples/`.

```{r}
>npartner_lefse_df <- data.frame(read.csv("path_to_the_package/KunDH-2023-CRM-MSM_metagenomics/examples/npartners_lefse_deviation_plot.tsv",
                                           header = TRUE,
                                           sep = "\t"))
```

Use [ggdotchart](https://rpkgs.datanovia.com/ggpubr/reference/ggdotchart.html) implemented in `ggpubr` for visualization.
```{r}
>library(ggpubr)
>ggdotchart(npartner_lefse_df, x = "feature", y = "lda_score",
           color = "class",
           palette = c("#0073C2FF", "#0073C2FF")
           sorting = "descending",                       
           add = "segments",                            
           add.params = list(color = "lightgray", size = 1.5),
           group = "class",             
           rotate = T,
           dot.size = 4,         
           shape = "class",   
           ggtheme = theme_pubr()     
) +  theme(text = element_text(size = 13, family = "Arial")) +  scale_x_discrete(position = "top")
```

![LefSe biomarkers linked with #partners](../images/npartner_lefse_biomarkers.png)

### UpSet plot

#### R packages required

* [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

In this section, we will show you how to visualize LefSe biomarkers associated with multiple groups using UpSet plot.
Our [example data](../example_data/UpSet_matrix1.tsv) comprises LefSe biomarkers associated with sexual practices including RAI: Yes (receiving anal intercourse), having >3 sexual partners (# partners: >3), practicing oral sex (Oral sex: Yes), diagnosed with sexually transmitted infection (STI: Positive), condomless during RAI (Condom use (during RAI): No).
First of all, open a new R working script, and load our [example data](../example_data/UpSet_matrix1.tsv) from `path_to_the_package/KunDH-2023-CRM-MSM_metagenomics/examples/`.

```{r}
>library("ComplexHeatmap")
>upset_matrix <- data.frame(read.csv("/Users/kunhuang/repos/KunDH-2023-CRM-MSM_metagenomics/example_data/UpSet_matrix1.tsv",
                                    header = TRUE,
                                    sep = "\t"))
>rownames(upset_matrix) <- upset_matrix[, colnames(upset_matrix)[[1]]]
>upset_matrix[, colnames(upset_matrix)[[1]]] <- NULL
>upset_matrix <- upset_matrix[as.logical(rowSums(upset_matrix != 0)),] # This step is optional.
```

Once the data is loaded, we use function `UpSet()` implemented in the package [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) to draw an UpSet plot.

```{r}
comb <- make_comb_mat(upset_matrix, mode = "intersect") # generate combination data
c_size <- comb_size(comb) # find combination sizes for setting order later
sets = c("RAI.yes", "X.partners.3", "Oral.yes",  "STI.positive", "condom.no") # manually set the set order

upset_plot <- ComplexHeatmap::UpSet(comb,
                                    comb_col = "#fb5238", # the color for combination columns
                                    bg_col = "#ffbeab", # the color for background of columns 
                                    bg_pt_col = "#ffdfd5", # the color for background of column patches
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
```

Once the backbone is generated, next we will display the values on the barplot using function `decorate_annotation`:

```{r}
upset_plot <- draw(upset_plot)
col_order <- column_order(upset_plot)

decorate_annotation("# shared taxonomic biomarkers", {
  grid.text(c_size[col_order], x = seq_along(c_size), y = unit(c_size[col_order], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
```

There are way more features to explore around for making an UpSet plot, please visit the [manual](https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html).


Similar codes can be used for another set of LefSe biomarkers associated with: RAI: No (not receiving anal intercourse), having 0-3 sexual partners (# partners: 0-3), not practicing oral sex (Oral sex: No), free from sexually transmitted infection (STI: Negative), use condom during RAI (Condom use (during RAI): Yes).

We can generate a nice plot showing the mutual biomarkers shared between different sexual practices by combining these plots.
![UpSet plot](../images/upset_plot.png)

