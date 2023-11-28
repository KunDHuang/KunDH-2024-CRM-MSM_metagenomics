
plot_complex_heatmap <- function(mat_file,
                                 column_md = NULL,
                                 row_md = NULL,
                                 color_bar_name = NULL,
                                 transformation = NULL,
                                 font_style = "Arial",
                                 font_size = 11,
                                 show_col_names = TRUE,
                                 show_row_names = TRUE,
                                 row_names_side = "left",
                                 column_names_side = "bottom",
                                 cluster_columns = FALSE,
                                 cluster_rows = FALSE,
                                 cluster_row_slices = FALSE,
                                 cluster_column_slices = FALSE,
                                 color_func = NULL,
                                 border = FALSE,
                                 row_gap = 1,
                                 column_gap = 1,
                                 width = 1,
                                 height = 1){
    mat_df <- data.frame(read.csv(mat_file, header = TRUE, sep = "\t"))
    rownames(mat_df) <- make.names(mat_df$X, unique = TRUE)
    mat_df$X <- NULL
    mat <- data.matrix(mat_df)
    if (is.null(column_md)) {
        column_md <- NULL
    } else {
        column_md <- scan(column_md, character(), quote = "", sep = "\n")
        column_md <- factor(column_md, levels = unique(column_md))
    }

    if (is.null(row_md)) {
        row_md <- NULL
    } else {
        row_md <- scan(row_md, character(), quote = "", sep = "\n")
        row_md <- factor(row_md, levels = unique(row_md))
    }

    if (is.null(transformation)) {
        ch_mat <- mat
    } else if (transformation == "sqrt_asin") {
        ch_mat <- asin(sqrt(mat/100))
    } else if (transformation == "log10") {
        ch_mat <- log10(mat + 0.000000001)
    } else if (transformation == "binary") {
        ch_mat <- mat
        ch_mat[ch_mat > 0] = 1
    } else {
        print("Please choose the function for transforming relative abundances: [sqrt_asin] or [log10] or NULL.")
    }
    
    if (row_names_side == "left") {
        row_title_side <- "right"
    } else if (row_names_side == "right") {
        row_title_side <- "left"
    } else {
        print("[row_names_side] has to be [right] or [left].")
    }

    if (column_names_side == "bottom") {
        column_title_side <- "top"
    } else if (column_names_side == "top") {
        column_title_side <- "bottom"
    } else {
        print("[column_names_side] has to be [bottom] or [top].")
    }

    row_names_gp <- grid::gpar(fontsize = font_size, fontfamily = font_style)
    column_names_gp <- grid::gpar(fontsize = font_size, fontfamily = font_style)
    row_title_gp <- grid::gpar(fontsize = font_size, fontfamily = font_style)
    column_title_gp <- grid::gpar(fontsize = font_size, fontfamily = font_style)
    
    row_order <- 1: nrow(ch_mat)
    width <- ncol(ch_mat)*grid::unit(width, "mm")
    height <- nrow(ch_mat)*grid::unit(height,"mm")
    row_gap <- grid::unit(row_gap, "mm")
    column_gap <- grid::unit(column_gap, 'mm')
    ch_plot <- ComplexHeatmap::Heatmap(ch_mat,
                                       name = color_bar_name,
                                       row_split = row_md,
                                       column_split = column_md,
                                       show_column_names = show_col_names,
                                       show_row_names = show_row_names,
                                       row_names_gp = row_names_gp,
                                       column_names_gp = column_names_gp,
                                       row_title_gp = row_title_gp,
                                       column_title_gp = column_title_gp,
                                       row_order = row_order,
                                       col = color_func,
                                       row_names_side = row_names_side,
                                       row_title_side = row_title_side,
                                       column_names_side = column_names_side,
                                       column_title_side = column_title_side,
                                       cluster_columns = cluster_columns,
                                       cluster_column_slices = cluster_column_slices,
                                       cluster_rows = cluster_rows,
                                       cluster_row_slices = cluster_column_slices,
                                       border = border,
                                       width = width,
                                       height = height,
                                       row_gap = row_gap,
                                       column_gap = column_gap)
    ch_plot
    }