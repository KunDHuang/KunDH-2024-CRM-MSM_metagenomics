
generate_coordis_df <- function(mat, md, dist_method = "bray") {
  # mat: the loaded matrix from mpa-style dataframe.
  # md: the dataframe containing metadata.
  # dist_method: the method for calculating beta diversity. default: ["bray"]. For other method, refer to vegdist()
  # this function is to prepare metadata-added coordinates dataframe.
  bray_dist <- vegan::vegdist(mat, dist_method)
  coordinates <- as.data.frame(ape::pcoa(bray_dist)$vectors)
  coor_df <- cbind(coordinates, md)
  coor_df
}

pcoa_plot <- function(mat,
                      md, 
                      dist_method, 
                      variable, 
                      fsize = 11, 
                      dsize = 1, 
                      fstyle = "Arial", 
                      to_rm = NULL) {
  # mat: the loaded matrix from mpa-style dataframe, [dataframe].
  # md: the dataframe containing metadata, [dataframe].
  # dist_method: the method for calculating beta diversity, [string]. default: ["bray"]. For other method, refer to vegdist(). 
  # fsize: the font size, [int].
  # dsize: the dot size, [int].
  # fstyle: the font style, [string].
  # variable: specify the variable name for separating groups, [string].
  # to_rm: a vector of values in "variable" column where the corresponding rows will be removed first.
  # this function is to draw pcoa plot with confidence ellipse
  coordis_df <- generate_coordis_df(mat, md, dist_method)
  if (is.null(to_rm)) {
    coordis_df <- coordis_df[!(is.na(coordis_df[, variable]) | coordis_df[, variable] == ""), ]
  }
  else {
    coordis_df <- coordis_df[!(is.na(coordis_df[, variable]) | coordis_df[, variable] == "" | coordis_df[, variable] %in% to_rm), ]
  }
  eval(substitute(ggplot(coordis_df, aes(Axis.1, Axis.2, color = c)),list(c = as.name(variable)))) +
    geom_point(size = dsize) + 
    theme_bw() +
    eval(substitute(geom_polygon(stat = "ellipse", aes(fill = c), alpha = 0.1), list(c = as.name(variable)))) +
    labs(x = "PC1", y = "PC2") +
    theme(text = element_text(size = fsize, family = fstyle)) +
    theme(legend.position="bottom") 
}

est_permanova <- function(mat, 
                          md, 
                          variable, 
                          covariables = NULL, 
                          nper = 999, 
                          to_rm = NULL, 
                          by_method = "margin"){
  # mat: the loaded matrix from mpa-style dataframe, [dataframe].
  # md: the dataframe containing metadata, [dataframe].
  # variable: specify the variable for testing, [string].
  # covariables: give a vector of covariables for adjustment, [vector].
  # nper: the number of permutation, [int], default: [999].
  # to_rm: a vector of values in "variable" column where the corresponding rows will be removed first.
  # by_method: "terms" will assess significance for each term, sequentially; "margin" will assess the marginal effects of the terms.
  if (is.null(to_rm)) {
    clean_md <- md[!(is.na(md[, variable]) | md[, variable] == ""), ]
  } else {
    clean_md <- md[!(is.na(md[, variable]) | md[, variable] == "" | md[, variable] %in% to_rm), ]
  }
  clean_idx = rownames(clean_md)
  clean_mat <- mat[rownames(mat) %in% clean_idx, ]
  if (is.null(covariables)) {
    est <- eval(substitute(adonis2(mat ~ cat, data = md, permutations = nper, by = by_method), list(cat = as.name(variable))))
  } else {
    mat_char <- deparse(substitute(mat))
    str1 <- paste0(c(variable, paste0(covariables, collapse = " + ")), collapse = " + ")
    str2 <- paste0(c(mat_char, str1), collapse = " ~ ")
    est <- vegan::adonis2(eval(parse(text = str2)), data = md, permutations = nper, by = by_method)
  }
  est
}