#### summarize fn ####
summary.ov <- function(object, model_results, sig_level=0.05, progress = TRUE, ...){
  temp = prep_for_plots(object, p_contours = stats::coef(summary(model_results$mod_results))[2,4])
  raw_treat = object$trt_effect[which(object$es_grid<.000000001 & object$es_grid>-.000000001),
                            which(object$rho_grid==0)]
  raw_pval = object$p_val[which(object$es_grid<.000000001 & object$es_grid>-.000000001),
                      which(object$rho_grid==0)]

  # Find effects and pvalues
  pvals = rep(NA, nrow(temp$obs_cors))
  trt_effect = rep(NA, nrow(temp$obs_cors))
  # set status bar
  pb <- progress::progress_bar$new(
    format = "  running simulation [:bar] :percent completed in :elapsed",
    total = nrow(temp$obs_cors), clear = FALSE, width= 60)
  for(i in 1:nrow(temp$obs_cors)){
    calculate_exact = ov_sim(model_results = model_results,
                                     plot_covariates = object$cov,
                                     rho_grid = temp$obs_cors$Cor_Outcome[i],
                                     es_grid=temp$obs_cors$ES[i],
                                     n_reps = object$n_reps)
    trt_effect[i] = calculate_exact$trt_effect[[1]]
    pvals[i] = calculate_exact$p_val[[1]]
    if(progress == TRUE){
      pb$tick()
      Sys.sleep(1 / nrow(temp$obs_cors))
    }
  }

  # Effect Size cases
  effect_size_text = dplyr::case_when((raw_treat < 0 & all(trt_effect < 0)) |
                                 raw_treat >= 0 & all(trt_effect >= 0) ~ "no sign changes",
                               raw_treat < 0 & all(trt_effect >= 0) |
                                 raw_treat >=0 & all(trt_effect < 0) ~ "all sign changes",
                               TRUE ~ "some sign changes")

  if(effect_size_text == "no sign changes"){
    diff_which.max = which.max(abs(trt_effect - raw_treat))
    most_extreme = trt_effect[diff_which.max]
    text = paste0("The sign of the estimated effect is expected to remain consistent when simulated unobserved confounders have the same strength of association with the treatment indicator and outcome that are seen in the observed confounders. In the most extreme observed case, the estimated effect size shifts from ", round(raw_treat,3), " to ", round(most_extreme,3), ".")

  } else if(effect_size_text == "all sign changes"){
    diff_which.max = which.max(abs(trt_effect - raw_treat))
    most_extreme = trt_effect[diff_which.max]
    text = paste0("The sign of the estimated effect is *not* expected to be robust to unobserved confounders that have the same strength of association with the treatment indicator and outcome that are seen in any of the observed confounders. In the most extreme observed case in which the sign changes, the estimated effect size shifts from ", round(raw_treat,3), " to ", round(most_extreme,3), ".")

  } else if(effect_size_text == "some sign changes"){
    # calculate sign changes
    if(raw_treat < 0){
      change = length(which(trt_effect >= 0))
      which.change = which(trt_effect >= 0)
      nochange = length(which(trt_effect < 0))
      total = nrow(temp$obs_cors)
      changes = temp$obs_cors$cov[which(trt_effect >=0)]
    } else{
      change = length(which(trt_effect < 0))
      which.change = which(trt_effect < 0)
      nochange = length(which(trt_effect >= 0))
      total = nrow(temp$obs_cors)
      changes = temp$obs_cors$cov[which(trt_effect < 0)]
    }
    diff_which.max = which.max(abs(trt_effect[which.change] - raw_treat))
    most_extreme = trt_effect[which.change][diff_which.max]
    text = paste0("The sign of the estimated effect is expected to remain consistent when simulated unobserved confounders have the same strength of associations with the treatment indicator and outcome that are seen in ", nochange, " of the ", total, " observed confounders. In the most extreme observed case in which the sign changes, the estimated effect size shifts from ", round(raw_treat,3), " to ", round(most_extreme,3), ". The sign of the estimate would not be expected to be preserved for unobserved confounders that have the same strength of association with the treatment indicator and outcome as ", paste(changes, collapse=", "), ".")
   }

  # Pvalue Size cases

  if(raw_pval < sig_level & all(pvals < sig_level)){
    text_p = paste0("Statistical significance at the ", sig_level, " level is expected to be robust to unobserved confounders with strengths of associations with the treatment indicator and outcome that are seen in the observed confounders. In the most extreme observed case, the p-value would be expected to increase from ", format(round(raw_pval,3), nsmall=3), " to ", format(round(max(pvals),3), nsmall=3), ".")
  } else if(raw_pval < sig_level & all(pvals > sig_level)){
    text_p = paste0("Statistical significance at the ", sig_level, " level is *not* expected to be robust to unobserved confounders with strengths of associations with the treatment indicator and outcome that are seen in any of the observed confounders. In the most extreme observed case, the p-value would be expected to increase from  ", format(round(raw_pval,3), nsmall=3), " to ", format(round(max(pvals),3), nsmall=3), ".")
  } else if(raw_pval < sig_level & !(all(pvals<sig_level) | all(pvals >sig_level))){
    nonsig = temp$obs_cors$cov[which(pvals>sig_level)]
    total = nrow(temp$obs_cors)
    sig_count = total - length(nonsig)
    text_p = paste0("Statistical significance at the ", sig_level, " level is expected to be robust to unobserved confounders with strengths of associations with the treatment indicator and outcome that are seen in ", sig_count, " of the ", total, " observed confounders. In the most extreme observed case, the p-value would be expected to increase from ", format(round(raw_pval,3), nsmall=3), " to ", format(round(max(pvals),3), nsmall=3), ". Significance at the ", sig_level, " level would not be expected to be preserved for unobserved confounders that have the same strength of association with the treatment indicator and outcome as ",  paste(nonsig, collapse=", "),".")
  }
  print("Recommendation for reporting the sensitivity analyses")

  if(raw_pval < 0.05){
    print(text)
    print(text_p)
  } else{
    print(text)
  }
}
