add_reps <- function(OVtool_results, model_results, more_reps){
  OVtool_results_new = ov_sim(model_results = model_results,
                              plot_covariates = OVtool_results$cov,
                              es_grid = OVtool_results$es_grid,
                              rho_grid = OVtool_results$rho_grid,
                              n_reps = more_reps,
                              add = TRUE,
                              progress = TRUE,
                              sim_archive = OVtool_results$es_se_raw)
  return(OVtool_results_new)
}
