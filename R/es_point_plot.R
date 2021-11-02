#### Treatment effect plot with points ####
es_point_plot = function(prep, col = "color"){
  r1_df = prep$r1_df
  pvals = prep$pvals
  pval_lines = prep$pval_lines
  text_high = prep$text_high
  obs_cors = prep$obs_cors
  p_labels = pvals
  if(length(which(pval_lines=="twodash")) > 1){
    p_labels[which(pval_lines=="twodash")] = "Others"
  }

  # Black and white or color
  color=rep(NA,2)
  color = ifelse(col == "color", list(c("red", "blue")),
                 ifelse(col == "bw", list(c("grey70", "black")),
                        stop("Please specify 'bw' or 'color'.")))
  redregion = "white"

  raw = prep$raw

  v3 <- ggplot2::ggplot(r1_df, ggplot2::aes(x = .data$es_grid, y = .data$rho_grid,
                                            z = .data$trt_effect)) +
    ggplot2::geom_rect(ggplot2::aes(ymin = max(r1_df$rho_grid),
                                    ymax = Inf,
                                    xmin = -Inf,
                                    xmax = Inf),
                       fill = redregion,
                       alpha = 0.01) +
    ggplot2::geom_rect(ggplot2::aes(ymin = -Inf,
                                    ymax = 0,
                                    xmin = -Inf,
                                    xmax = Inf),
                       fill = "snow2",
                       alpha = 0.01) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0,max(c(obs_cors$Cor_Outcome_Actual,
                                                   r1_df$rho_grid))),expand = c(0, 0)) +
    ggplot2::theme_linedraw() +
    ggplot2::geom_hline(yintercept=max(r1_df$rho_grid)) +
    ggplot2::geom_contour(col='black') + ggplot2::xlab("Association with Treatment Indicator\n(effect size scale)") +
    ggplot2::ylab("Absolute Correlation with Outcome (rho)") +
    ggplot2::ggtitle(paste0("Estimated treatment effect contours\n", raw)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
          legend.key = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 10),
          legend.key.size =  grid::unit(0.5, "in"),
          panel.grid.major = ggplot2::element_line(color="snow2"),
          panel.grid.minor = ggplot2::element_line(color="snow2")) +
    ggplot2::geom_vline(xintercept = 0, color = 'darkgrey')
  for(i in 1:length(pvals)){
    v3 = v3 +
      ggplot2::geom_contour(data = r1_df, ggplot2::aes(x = .data$es_grid,
                                                       y = .data$rho_grid,
                                                       z = .data$p_val,
                                                       linetype=!!pval_lines[i]),
                            color = color[[1]][1],
                            breaks=pvals[i],
                            size=1.25)
  }
  v3 = v3 +
    ggplot2::scale_linetype_manual(name = 'P-value Threshold', values = pval_lines,
                          labels = p_labels) +
    ggplot2::geom_point(data = obs_cors, col=color[[1]][2],
                        ggplot2::aes(x = .data$ES,
                                     y = .data$Cor_Outcome_Actual,
                                     z = NULL)) +
    # ggplot2::annotation_custom(grob = grid::textGrob(label = raw,
    #                                                  vjust = 3,
    #                                                  gp = grid::gpar(cex = .75)),
    #                   ymin = .1,
    #                   ymax = .1,
    #                   xmax = max(r1_df$es_grid)+.1,
    #                   xmin=max(r1_df$es_grid)+.05)  +
    ggrepel::geom_text_repel(data = obs_cors,
                             ggplot2::aes(x = .data$ES,
                                          y = .data$Cor_Outcome_Actual,
                                          z = NULL,
                                          label = .data$cov),
                    box.padding = grid::unit(0.45, "lines"),
                    col=color[[1]][2]) +
    metR::geom_text_contour(ggplot2::aes(z=.data$trt_effect),
                            stroke=.2,
                            check_overlap = T)

  if(col == "bw"){
    v3 = v3 + ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.key = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 10),
                     legend.key.size =  grid::unit(0.5, "in"))
  }
  return(v3)
}

