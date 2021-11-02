\name{summary.ov}
\alias{summary.ov}
\title{summary.ov
}
\description{
Produces summary information that contains a recommendation for reporting the sensitivity analyses
}
\usage{
\method{summary}{ov}(object, model_results, sig_level = 0.05, progress = TRUE, \dots)
}
\arguments{
  \item{object}{
The object returned from OVtool::ov_simgrid()
}
  \item{model_results}{
The object returned from OVtool::outcome_model()
}
  \item{sig_level}{
The alpha level with default 0.05
}
  \item{progress}{
Whether or not the function progress should print to screen. The default value is TRUE. If the user does not want the output to print to screen, they should set to FALSE.
}
  \item{...}{
Additional arguments.
}
}
\value{
This function will print a recommendation for reporting the sensitivity analyses.
}
\examples{
data(sud)
sud = data.frame(sud)
sud$treat = ifelse(sud$treat == "A", 1, 0)
sud$wts = sample(seq(1, 10, by=.01), size=nrow(sud), replace = TRUE)
outcome_mod = outcome_model(data = sud,
                            weights = "wts",
                            treatment = "treat",
                            outcome = "eps7p_6",
                            model_covariates = c("sfs8p_0", "eps7p_0"),
                            estimand = "ATE")

ovtool_results = ov_sim(model_results=outcome_mod,
                        plot_covariates=c("sfs8p_0"),
                        es_grid = NULL,
                        rho_grid = NULL,
                        n_reps = 2,
                        progress=FALSE)
summary = summary.ov(object = ovtool_results,
                     model_results = outcome_mod,
                     sig_level=0.05,
                     progress = FALSE)
}

