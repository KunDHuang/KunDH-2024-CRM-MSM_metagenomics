
SE_converter <- function(md_rows, tax_starting_row, mpa_md) {
    # SE_converter function is to convery metadata-wedged mpa table into SummarisedExperiment structure.
    # md_rows: a vector specifying the range of rows indicating metadata.
    # tax_starting_row: an interger corresponding to the row where taxonomic abundances start.
    # mpa_md: a metaphlan table wedged with metadata, in the form of dataframe.
    md_df <- mpa_md[md_rows,] # extract metadata part from mpa_md table
    tax_df <- mpa_md[tax_starting_row: nrow(mpa_md),] # extract taxonomic abundances part from mpa_md table
    
    ### convert md_df to a form compatible with SummarisedExperiment ### 
    SE_md_df <- md_df[, -1]
    rownames(SE_md_df) <- md_df[, 1]
    SE_md_df <- t(SE_md_df)
    ### convert md_df to a form compatible with SummarisedExperiment ###
    
    ### prep relab values in a form compatible with SummarisedExperiment ###
    SE_relab_df <- tax_df[, -1]
    rownames(SE_relab_df) <- tax_df[, 1]
    col_names <- colnames(SE_relab_df)
    SE_relab_df[, col_names] <- apply(SE_relab_df[, col_names], 2, function(x) as.numeric(as.character(x)))
    ### prep relab values in a form compatible with SummarisedExperiment ###

    SE_tax_df <- tax_df[, 1:2]
    rownames(SE_tax_df) <- tax_df[, 1]
    SE_tax_df <- SE_tax_df[-2]
    colnames(SE_tax_df) <- c("species")
    
    SE_data <- SummarizedExperiment::SummarizedExperiment(assays = list(relative_abundance = SE_relab_df),
                                         colData = SE_md_df,
                                         rowData = SE_tax_df)

    SE_data
}

est_alpha_diversity <- function(se_data) {
    # This function is to estimate alpha diversity (shannon index and richness) of a microbiome and output results in dataframe.
    # se_data: the SummarizedExperiment data structure containing metadata and abundance values.
    se_data <- se_data |>
        mia::estimateRichness(abund_values = "relative_abundance", index = "observed")
    se_data <- se_data |>
        mia::estimateDiversity(abund_values = "relative_abundance", index = "shannon")
    se_alpha_div <- data.frame(SummarizedExperiment::colData(se_data))
    se_alpha_div
}

make_boxplot <- function(df, xlabel, ylabel, font_size = 11, jitter_width = 0.2, dot_size = 1, font_style = "Arial", stats = TRUE, pal = NULL) {
    # This function is to create a boxplot using categorical data.
    # df: The dataframe containing microbiome alpha diversities, e.g. `shannon` and `observed` with categorical metadata.
    # xlabel: the column name one will put along x-axis.
    # ylabel: the index estimate one will put along y-axis.
    # font_size: the font size, default: [11]
    # jitter_width: the jitter width, default: [0.2]
    # dot_size: the dot size inside the boxplot, default: [1]
    # font_style: the font style, default: `Arial`
    # pal: a list of color codes for pallete, e.g. c(#888888, #eb2525). The order corresponds the column order of boxplot.
    # stats: wilcox rank-sum test. default: TRUE
    if (stats) {
          nr_group = length(unique(df[, xlabel])) # get the number of groups
          if (nr_group == 2) {
              group_pair = list(unique(df[, xlabel]))
              ggpubr::ggboxplot(data = df, x = xlabel, y = ylabel, color = xlabel,
              palette = pal, ylab = ylabel, xlab = xlabel,
              add = "jitter", add.params = list(size = dot_size, jitter = jitter_width)) +
              ggpubr::stat_compare_means(comparisons = group_pair, exact = T, alternative = "less") +
              ggplot2::stat_summary(fun.data = function(x) data.frame(y = max(df[, ylabel]), label = paste("Mean=",mean(x))), geom="text") +
              ggplot2::theme(text = ggplot2::element_text(size = font_size, family = font_style))
           }
          else {
              group_pairs = my_combn(unique((df[, xlabel])))
              ggpubr::ggboxplot(data = df, x = xlabel, y = ylabel, color = xlabel,
              palette = pal, ylab = ylabel, xlab = xlabel,
              add = "jitter", add.params = list(size = dot_size, jitter = jitter_width)) +
              ggpubr::stat_compare_means() + ggpubr::stat_compare_means(comparisons = group_pairs, exact = T, alternative = "greater") +
              ggplot2::stat_summary(fun.data = function(x) data.frame(y= max(df[, ylabel]), label = paste("Mean=",mean(x))), geom="text") +
              ggplot2::theme(text = ggplot2::element_text(size = font_size, family = font_style))
           }
    }
    else {
        ggpubr::ggboxplot(data = df, x = xlabel, y = ylabel, color = xlabel,
        palette = pal, ylab = ylabel, xlab = xlabel,
        add = "jitter", add.params = list(size = dot_size, jitter = jitter_width)) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size, family = font_style))
    }
}

my_combn <- function(x) {
  combs <- list()
  comb_matrix <- combn(x, 2)
  for (i in 1: ncol(comb_matrix)) {
    combs[[i]]  <- comb_matrix[,i]
  }
  combs
}