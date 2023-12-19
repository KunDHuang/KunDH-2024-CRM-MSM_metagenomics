# Visualizations for microbiome biomarkers 
This tutorial is to visualize microbiome biomarkers with different plots.

#### R packages required

* [ggpubr](https://rpkgs.datanovia.com/ggpubr/)

#### Deviation plot

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