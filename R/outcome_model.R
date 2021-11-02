outcome_model <- function(ps_object=NULL, stop.method=NULL, data, weights=NULL,
                          treatment, outcome, model_covariates,
                          estimand = "ATE"){
  data$w_orig = NA
  if(!(estimand %in% c("ATT", "ATE"))){
    stop("At the moment, OVtool only handles ATT and ATE.")
  }
  if(class(ps_object)!="ps"){
    if(missing(data) | missing(weights) | missing(treatment) | missing(outcome) | missing(model_covariates)){
      stop("Please supply either a ps class object, stop.method, and relevant column names for:
      outcome, and covariates OR the data and relevant column names for: weights, treatment, outcome, & covariates")
    } else{
      # weights
      data$w_orig = data[,weights]
      # outcome
      y = outcome
      # treamtment
      tx = treatment
      # covariates
      if(typeof(model_covariates) == "language"){
        cov = all.vars(model_covariates)
      } else{
        cov = model_covariates
      }
      # survey design
      design_u <- survey::svydesign(ids=~1, weights=~w_orig, data=data)
    }
  } else{
    if(is.null(stop.method) | is.null(outcome) | is.null(model_covariates)){
      stop("Please supply a stop.method to generate the weights (e.g. \"ks.max/")
    }
    # outcome
    y = outcome
    # treatment
    tx = ps_object$treat.var
    # covariates
    if(typeof(model_covariates) == "language"){
      cov = all.vars(model_covariates)
    } else{
      cov = model_covariates
    }
    # data
    data = ps_object$data
    # weights
    data$w_orig = twang::get.weights(ps_object, stop.method = stop.method, estimand = estimand)
    # survey design
    design_u <- survey::svydesign(ids=~1, weights=~w_orig, data=data)
  }
  # create formula
  formula = formula(paste0(y, " ~ ", tx, " + ", glue::glue_collapse(cov, sep=" + ")))
  # run outcomes model:
  if(length(unique(data[, y ])) <=1){
    stop("No variation in outcome.")
  } else {
    outcome_mod_results = survey::svyglm(formula, design=design_u)
  }
  # print(summary(model_results))

  return(list(ps_object = ps_object, stop.method = stop.method, data = data,
              weights = data$w_orig, tx = tx, y = y, outcome_mod_fmla = formula,
              estimand = estimand, mod_results = outcome_mod_results))
}

