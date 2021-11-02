#### Plotting Package for OV ####
plot.ov <- function(x, col="color", print_graphic="1", p_contours = c(0.01, 0.05, 0.1), ...) {
  prep_plots = prep_for_plots(x, p_contours)

  #### some defaults ####
  plot_1 = es_plot(prep=prep_plots, col=col)
  plot_2 = add_pvals_plot(prep=prep_plots, col=col)
  plot_3 = es_point_plot(prep=prep_plots,col=col)

  plot_1 <- suppressWarnings(ggplot2::ggplotGrob(plot_1))
  plot_1$layout$clip[plot_1$layout$name == "panel"] <- "off"
  plot_2 <- suppressWarnings(ggplot2::ggplotGrob(plot_2))
  plot_2$layout$clip[plot_2$layout$name == "panel"] <- "off"
  plot_3 <- suppressWarnings(ggplot2::ggplotGrob(plot_3))
  plot_3$layout$clip[plot_3$layout$name == "panel"] <- "off"

  # .pardefault <- graphics::par(no.readonly=TRUE)
  if(print_graphic == "1"){
    suppressWarnings(grid::grid.draw(plot_1))
  } else if(print_graphic == "2"){
    suppressWarnings(grid::grid.draw(plot_2))
    if(prep_plots$text_high!=""){print(prep_plots$text_high)}
    if(prep_plots$text_high_es!=""){print(prep_plots$text_high_es)}
  } else if(print_graphic == "3"){
    suppressWarnings(grid::grid.draw(plot_3))
    if(prep_plots$text_high!=""){print(prep_plots$text_high)}
    if(prep_plots$text_high_es!=""){print(prep_plots$text_high_es)}
  }
  # graphics::par(.pardefault)
  # grDevices::palette("default")
}
