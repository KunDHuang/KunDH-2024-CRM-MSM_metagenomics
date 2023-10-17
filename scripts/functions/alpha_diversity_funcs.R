
SE_converter <- function(md_rows, tax_starting_row, mpa_md) {
    # SE_converter function is to convery metadata-wedged mpa table into SummarisedExperiment structure.
    # md_rows: a vector specifying the range of rows indicating metadata.
    # tax_starting_row: an interger corresponding to the row where taxonomic abundances start.
    # mpa_md: a metaphlan table wedged with metadata, in the form of dataframe.
    md_df <- mpa_md[md_rows,] # extract metadata part from mpa_md table
    tax_df <- mpa_md[tax_starting_row: nrow(mpa_df),] # extract taxonomic abundances part from mpa_md table
    
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