#### ov_sim fn ####
ov_sim <- function(model_results, plot_covariates,
                   es_grid = seq(-.4, .4, by = 0.05),
                   rho_grid = seq(0, .4, by = 0.05), n_reps = 50,
                   progress = TRUE,
                   add = FALSE,
                   sim_archive = NULL){
  tx = model_results$tx
  y = model_results$y
  dta = model_results$data
  estimand = model_results$estimand

  # covariates
  if(typeof(plot_covariates) == "language"){
    cov = all.vars(plot_covariates)
  } else{
    cov = plot_covariates
  }

  # create formula
  formula = model_results$outcome_mod_fmla

  # checks
  if(add != TRUE){
    if(n_reps <= 1){stop("Please specify at least two (2) n_reps. For example: n_reps = 50")}
    if(!is.null(rho_grid)){
      if(any(rho_grid <0) | any(rho_grid > 1)){
        stop("Please specify a rho grid between 0 and 1. For example: rho_grid = seq(0, .45, by = 0.05)")
      }
    }
    # if(!is.null(es_grid)){
    #   if(any(abs(es_grid) > 1)){
    #     stop("This tool only considers small to moderate correlations. An effect size of 1 corresponds to a correlation of 0.44. Please reduce the size of your grid. For example: es_grid = seq(-.4, .4, by = 0.05)")
    #   }
    # }

    if(!all(dta[,tx] %in% c(0,1))) stop("Treatment variable `tx` must be only 0/1 values.")

    # determine reasonable grid to simulate over
    if(length(es_grid) != 1){
      jdp_test=find_esgrid(my_data = dta,
                           my_cov = cov,
                           treatment = tx,
                           outcome = y,
                           my_estimand = estimand)
      # SET final es_grid.
      if(!is.null(es_grid)){
        es_drop = jdp_test[which(jdp_test$ES > max(abs(es_grid))),'cov']
        if(max(abs(es_grid)) < max(abs(jdp_test$ES))){
          warning(paste0("You specified an effect size grid whose maximum in absolute value is less than the maximum absolute association at least one observed covariate has with the treatment indicator. The effect size grid was automatically expanded to include all plot_covariates specified in the relevant graphics. If you want the effect size grid range to remain from ", min(es_grid), " to ", max(es_grid), " then you must exclude the following variables from the plot_covariates argument: ", paste(es_drop, collapse=", "), "."))
          es_upper = .05 * ceiling(max(jdp_test$ES)/.05)
          es_lower = -es_upper
          es_grid = seq(es_lower, es_upper, by=0.05)
        }
      } else{
        es_upper = .05 * ceiling(max(jdp_test$ES)/.05)
        es_lower = -es_upper
        es_grid = seq(es_lower, es_upper, by=0.05)
      }

      if(!is.null(rho_grid)){
        rho_drop = jdp_test[which(jdp_test$Cor_Outcome > max(abs(rho_grid))),'cov']
        if(max(rho_grid) < max(jdp_test$Cor_Outcome)){
          warning(paste0("You specified a rho grid whose maximum value is less than the maximum absolute correlation at least one observed covariate has with the outcome. The rho grid was automatically expanded to include all plot_covariates specified in the relevant graphics. If you want the rho grid range to remain from ", min(rho_grid), " to ", max(rho_grid), " then you must exclude the following variables from the plot_covariates argument: ", paste(rho_drop, collapse=", "), "."))
          rho_upper = .05 * ceiling(max(jdp_test$Cor_Outcome)/.05)
          rho_grid = seq(0, rho_upper, by=0.05)
        }
      } else{
        rho_upper = .05 * ceiling(max(jdp_test$Cor_Outcome)/.05)
        rho_grid = seq(0, rho_upper, by=0.05)
      }
    }

    # Setup simulation:
    trt_effect_nodr <- matrix(0,length(es_grid), length(rho_grid))
    p_val_nodr <- matrix(0,length(es_grid), length(rho_grid))
    std_error_nodr <- matrix(0, length(es_grid), length(rho_grid))
    pValHd <- esHd <- StdError <- rep(NA, n_reps)
    # save results:
    temp_store = expand.grid(es=es_grid, rho=rho_grid); store_results = data.frame()
    for(i in 1:n_reps){
      store_results = dplyr::bind_rows(store_results, temp_store)
    }
    store_results$esHd = NA; store_results$StdError = NA
    # create w_new and set it to the original weights for now
    dta$w_new = dta$w_orig
    # residuals
    my_res = as.vector(stats::lm(formula = model_results$outcome_mod_fmla, data = dta)$residuals)

    # set status bar:
    pb <- progress::progress_bar$new(
      format = "  running [:bar] :percent completed in :elapsed",
      total = length(es_grid), clear = FALSE, width= 60)
    for(i in 1:length(es_grid)){
      for(j in 1:length(rho_grid)){
        for(k in 1:n_reps){
          a_prep <- gen_a_start(y=dta[,y],
                                tx = dta[,tx],
                                residuals = my_res,
                                es = es_grid[i],
                                rho = rho_grid[j],
                                my_estimand = estimand)
          a <- gen_a_finish(a_res=a_prep, my_estimand=estimand)
          dta$w_new <- dta$w_orig * a
          design_u <- survey::svydesign(ids=~1, weights=~w_new, data=dta)
          glm0_u_nodr <- survey::svyglm(formula, design=design_u)
          esHd[k] <- summary(glm0_u_nodr)$coefficients[tx,1]
          StdError[k] <- summary(glm0_u_nodr)$coefficients[tx,2]
        }
        store_results[which(store_results$es == es_grid[i] &
                              store_results$rho == rho_grid[j]),c('esHd','StdError')] = data.frame(esHd, StdError)
        combine = Amelia::mi.meld(q=data.frame(esHd), se=data.frame(StdError))
        melded_summary <- as.data.frame(cbind(t(combine$q.mi),
                                              t(combine$se.mi))) %>%
          purrr::set_names(c("estimate", "std.error")) %>%
          data.frame() %>%
          tibble::rownames_to_column(var = "term") %>%
          tibble::as_tibble() %>%
          dplyr::select(.data$term, tidyselect::everything()) %>%
          dplyr::mutate(statistic = .data$estimate / .data$std.error,
                        p.value = 2 * stats::pnorm(abs(.data$statistic),lower.tail = FALSE))
        p_val_nodr[i,j] <- melded_summary$p.value
        trt_effect_nodr[i,j] <- melded_summary$estimate
        std_error_nodr[i,j] <- melded_summary$std.error
      }
      if((length(es_grid) >1) & progress == TRUE){
        # print(paste0(round(i/length(es_grid)*100,0), "% Done!"))
          pb$tick()
          Sys.sleep(1 / length(es_grid))
      }
    }
  } else{
    if(n_reps <= 1){stop("Please specify at least one (1) more_reps. For example: more_reps = 10")}

    # Setup simulation:
    trt_effect_nodr <- matrix(0,length(es_grid), length(rho_grid))
    p_val_nodr <- matrix(0,length(es_grid), length(rho_grid))
    std_error_nodr <- matrix(0, length(es_grid), length(rho_grid))
    pValHd <- esHd <- StdError <- rep(NA, n_reps)
    temp_store = expand.grid(es=es_grid, rho=rho_grid); store_results = data.frame()
    for(i in 1:n_reps){
      store_results = dplyr::bind_rows(store_results, temp_store)
    }
    store_results$esHd = NA; store_results$StdError = NA
    dta$w_new = dta$w_orig
    # residuals:
    my_res = as.vector(stats::lm(formula = model_results$outcome_mod_fmla, data = dta)$residuals)

    # set status bar
    pb <- progress::progress_bar$new(
      format = "  running simulation [:bar] :percent completed in :elapsed",
      total = length(es_grid), clear = FALSE, width= 60)
    for(i in 1:length(es_grid)){
      for(j in 1:length(rho_grid)){
        for(k in 1:n_reps){
          a_prep <- gen_a_start(y=dta[,y],
                                tx = dta[,tx],
                                residuals = my_res,
                                es = es_grid[i],
                                rho = rho_grid[j],
                                my_estimand = estimand)
          a <- gen_a_finish(a_res=a_prep, my_estimand = estimand)
          dta$w_new <- dta$w_orig * a
          design_u <- survey::svydesign(ids=~1, weights=~w_new, data=dta)
          glm0_u_nodr <- survey::svyglm(formula, design=design_u)
          esHd[k] <- summary(glm0_u_nodr)$coefficients[tx,1]
          StdError[k] <- summary(glm0_u_nodr)$coefficients[tx,2]
        }
        pos = which(sim_archive$es == es_grid[i] & sim_archive$rho == rho_grid[j])
        esHd_new = c(sim_archive[pos,'esHd'], esHd)
        StdError_new = c(sim_archive[pos,'StdError'], StdError)
        store_results[which(store_results$es == es_grid[i] &
                              store_results$rho == rho_grid[j]),c('esHd','StdError')] = data.frame(esHd, StdError)
        combine = Amelia::mi.meld(q=data.frame(esHd_new), se=data.frame(StdError_new))
        melded_summary <- as.data.frame(cbind(t(combine$q.mi),
                                              t(combine$se.mi))) %>%
          purrr::set_names(c("estimate", "std.error")) %>%
          data.frame() %>%
          tibble::rownames_to_column(var = "term") %>%
          tibble::as_tibble() %>%
          dplyr::select(.data$term, tidyselect::everything()) %>%
          dplyr::mutate(statistic = .data$estimate / .data$std.error,
                        p.value = 2 * stats::pnorm(abs(.data$statistic),lower.tail = FALSE))
        p_val_nodr[i,j] <- melded_summary$p.value
        trt_effect_nodr[i,j] <- melded_summary$estimate
        std_error_nodr[i,j] <- melded_summary$std.error
      }
      if((length(es_grid) >1) & progress == TRUE){
        # print(paste0(round(i/length(es_grid)*100,0), "% Done!"))
        pb$tick()
        Sys.sleep(1 / length(es_grid))
      }
    }
    store_results = dplyr::bind_rows(sim_archive, store_results)
  }
  results = list(p_val = p_val_nodr, trt_effect = trt_effect_nodr,
                 es_grid = es_grid, rho_grid = rho_grid, cov = cov,
                 data = dta, tx = tx, y = y, estimand = estimand,
                 n_reps = n_reps, std.error = std_error_nodr,
                 es_se_raw = store_results)
  class(results) <- "ov"
  return(results)
}
