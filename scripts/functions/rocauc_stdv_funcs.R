data_summary <- function(data,
                         var_estimates,
                         groups){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groups, .fun=summary_func,
                  var_estimates)
  data_sum <- rename(data_sum, c("mean" = var_estimates))
  return(data_sum)
}

std_deviation_plot <- function(df,
                               x_axis,
                               y_axis,
                               stdv_column,
                               palette = "jco",
                               y_label = "AUC-ROC",
                               x_label = "",
                               order_x_axis = NULL,
                               font_size = 11,
                               font_family = "Arial"){
    stdv_plot <- ggpubr::ggdotplot(data = df, x = x_axis, y = y_axis, color = x_axis, fill = x_axis,
                           palette = palette, ylab = y_label, xlab = x_label,
                           order = order_x_axis) +
                           ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted", col = 'red') +
                           ggplot2::theme(text = ggplot2::element_text(size = font_size, family = font_family)) +
                           ggplot2::theme(legend.title = ggplot2::element_blank()) +
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = eval(parse(text = y_axis)) - eval(parse(text = stdv_column)),
                                                    ymax =  eval(parse(text = y_axis)) + eval(parse(text = stdv_column)),
                                                    color = eval(parse(text = x_axis))), width = .2,
                                                    position = ggplot2::position_dodge(0.7))
    stdv_plot
}