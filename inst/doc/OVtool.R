## ----setup, include=F, linewidth=60-------------------------------------------
# github_document
# html_document
# word_document
# pdf_document
library(knitr)
knitr::opts_chunk$set(echo = TRUE,
                      fig.path = "inst/",
                      comment = "#>")
hook_output = knit_hooks$get('output')
knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x = knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    if (any(nchar(x) > n)) x = strwrap(x, width = n)
    x = paste(x, collapse = '\n')
  }
  hook_output(x, options)
})
time = Sys.time()

## ---- begin, message=F, linewidth=60------------------------------------------
library(OVtool)

## ---- data--------------------------------------------------------------------
set.seed(24)
data(sud); sud = data.frame(sud)
sud$treat = ifelse(sud$treat == "A", 1, 0)

## ---- estTWANG, message=F, linewidth=60---------------------------------------
## Create Formula
my_formula = as.formula(treat ~ eps7p_0 + sfs8p_0 + sati_0 + ada_0 + recov_0 + 
                          tss_0 + dss9_0)

## Get weights
# library(twang)
ps.twang <- twang::ps(my_formula, data = sud, estimand = 'ATE', booster = "gbm",
                      stop.method = "ks.max", verbose=F, ks.exact = T)

## ---- balance, message=F------------------------------------------------------
# Check Balance
twang::bal.table(ps.twang)

## ---- getweights, message=F, linewidth=60-------------------------------------
# Get weights (not needed if user inserts a ps object in OVtool::outcome_model)
sud$w_twang = ps.twang$w$ks.max.ATE

## ---- outcomemodel------------------------------------------------------------
# Run Models -- first standardize outcome
sud$eps7p_3_std = sud$eps7p_3/sd(sud$eps7p_3) 

# Run outcome model (function in OVtool that calls svyglm from survey)
results = outcome_model(ps_object = NULL,
                        stop.method = NULL, 
                        data = sud,
                        weights = "w_twang", 
                        treatment = "treat",
                        outcome = "eps7p_3_std", 
                        model_covariates = c("eps7p_0", "sfs8p_0",
                                             "sati_0", "ada_0",
                                             "recov_0", "tss_0",
                                             "dss9_0"),
                        estimand = "ATE")

summary(results$mod_results)

## ---- ov_sim, linewidth=60, eval=F--------------------------------------------
#  # Run OVtool (with weights (not a ps object))
#  ovtool_results_twang = ov_sim(model_results=results,
#                                plot_covariates=c("eps7p_0", "sfs8p_0",
#                                                  "sati_0", "ada_0",
#                                                  "recov_0", "tss_0",
#                                                  "dss9_0"),
#                                es_grid = NULL,
#                                rho_grid = seq(0, 0.40, by = 0.05),
#                                n_reps = 50,
#                                progress = TRUE,
#                                add = FALSE,
#                                sim_archive = NULL)

## ---- ov_sim_output, linewidth=60, echo=F-------------------------------------
# save(ovtool_results_twang, file = "vignettes/ATE_ovtool.rda", version = 2)
load("ATE_ovtool.rda")

## ---- fig1, fig.width=10, fig.height=7, fig.align='center', warning=F, linewidth=60----
plot.ov(ovtool_results_twang, print_graphic = "1", col = "bw")

## ---- fig2, fig.width=10, fig.height=7, fig.align='center', linewidth=60------
plot.ov(ovtool_results_twang, print_graphic = "2", col = "color")

## ---- fig3, fig.width=10, fig.height=7, fig.align='center', linewidth=60------
plot.ov(ovtool_results_twang, print_graphic = "3", col = "color")

## ---- addreps, linewidth=60---------------------------------------------------
# If you want to add repetitions, run the following line.
# ovtool_results_twang = add_reps(OVtool_results = ovtool_results_twang,
#                                 model_results = results,
#                                 more_reps = 30)
#
# Recreate Graphic
# plot.ov(ovtool_results_twang, print_graphic = "1", col = "bw")

## ---- summary`, linewidth=60--------------------------------------------------
summary.ov(object = ovtool_results_twang, model_results = results)

## ---- ATT---------------------------------------------------------------------
## Create formula - here we use same formula as in the ATE case.
my_formula = as.formula(treat ~ eps7p_0 + sfs8p_0 + sati_0 + ada_0 + recov_0 + 
                          tss_0 + dss9_0)

# Propensity score weights:
ps.twang_att <- twang::ps(my_formula, data = sud, estimand = 'ATT',
                          booster = "gbm", stop.method = "ks.max", 
                          verbose=F, ks.exact = T)
twang::bal.table(ps.twang_att)

## ---- outcome_att-------------------------------------------------------------
results_att = outcome_model(ps_object = ps.twang_att,
                            stop.method = "ks.max",
                            data = sud,
                            weights = NULL, 
                            treatment = "treat",
                            outcome = "eps7p_3_std",
                            model_covariates = c("eps7p_0", "sfs8p_0",
                                                 "sati_0", "ada_0",
                                                 "recov_0", "tss_0",
                                                 "mhtrt_0", "dss9_0"),
                            estimand = "ATT")
summary(results_att$mod_results)

## ---- ovsim_att, linewidth=60, eval=F-----------------------------------------
#  ovtool_results_twang_att = ov_sim(model_results=results_att,
#                                    plot_covariates=c("eps7p_0", "sfs8p_0",
#                                                      "sati_0", "ada_0",
#                                                      "recov_0", "tss_0",
#                                                      "dss9_0"),
#                                    es_grid = NULL,
#                                    rho_grid = seq(0, 0.40, by = 0.05),
#                                    n_reps = 50,
#                                    progress = TRUE)

## ---- ovsim_att_output, linewidth=60, echo=F----------------------------------
# save(ovtool_results_twang_att, file = "vignettes/ATT_ovtool.rda", version = 2)
load("ATT_ovtool.rda")

## ---- fig1_att, fig.width=10, fig.height=7, fig.align='center', warning=F, linewidth=60----
plot.ov(ovtool_results_twang_att, print_graphic = "1", col = "bw")

## ---- addreps_att, linewidth=60, eval=F---------------------------------------
#  # If you want to add repetitions, run the following line.
#  ovtool_results_twang_att = add_reps(OVtool_results = ovtool_results_twang_att,
#                                      model_results = results_att,
#                                      more_reps = 10)

## ---- addreps_att_output, linewidth=60, echo=F--------------------------------
# save(ovtool_results_twang_att, file = "vignettes/ATT_ovtool2.rda", version = 2)
load("ATT_ovtool2.rda")

## ---- fig2_att_new, fig.width=10, fig.height=7, fig.align='center', warning=F, linewidth=60----
plot.ov(ovtool_results_twang_att, print_graphic = "1", col = "bw")

## ---- fig2_att, fig.width=10, fig.height=7, fig.align='center', warning=F, linewidth=60----
plot.ov(ovtool_results_twang_att, print_graphic = "2", col = "color")

## ---- fig3_att, fig.width=10, fig.height=7, fig.align='center', warning=F, linewidth=60----
plot.ov(ovtool_results_twang_att, print_graphic = "3", col = "color")

## ---- summary_att, linewidth=60-----------------------------------------------
summary.ov(object = ovtool_results_twang_att, model_results = results_att)

