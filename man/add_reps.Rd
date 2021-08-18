\name{add_reps}
\alias{add_reps}
\title{add_reps
}
\description{
This function will run additional simulations of the unobserved confounder and add the results to the object returned from the previous call to OVtool::ov_sim
}
\usage{
add_reps(OVtool_results, model_results, more_reps)
}
\arguments{
  \item{OVtool_results}{
The object returned from OVtool::ov_sim()
}
  \item{model_results}{
The object returned from OVtool::outcome_model()
}
  \item{more_reps}{
The number of additional repetitions the user wants to simulate the unobserved confounder
}
}

\value{
add_reps returns an updated object returned from OVtool::ov_sim()
}
\examples{
data(sud)
sud = data.frame(sud[sample(1:nrow(sud),100),])
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
                        es_grid = .1,
                        rho_grid = .1,
                        n_reps = 2,
                        progress=FALSE)

additional = add_reps(OVtool_results = ovtool_results,
                      model_results = outcome_mod,
                      more_reps = 2)
}


