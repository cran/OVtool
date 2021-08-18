\name{ov_sim}
\alias{ov_sim}
\title{ov_sim
}
\description{
This function will create the simulation grid. The simulation will iterate over
effects sizes and absolute correlations with the outcome (rho) and see how the
treatment effect and relevant p-value changes
}
\usage{
ov_sim(model_results, plot_covariates, es_grid = seq(-.4, .4, by = 0.05),
rho_grid = seq(0, .4, by = 0.05), n_reps = 50, progress = TRUE, add = FALSE,
sim_archive = NULL)
}
\arguments{
  \item{model_results}{
  object returned from outcome_model
}
  \item{plot_covariates}{
  vector of column names representing the covariates that will be plotted on the graphic as observed covariates. Most users will include the variables on the right-hand side of the propensity score model
}
  \item{es_grid}{
  Not required. A grid of effect sizes to simulate over
}
  \item{rho_grid}{
  Not required. A grid of correlations to simulate over; rho relates the correlation to the effect size.
}
  \item{n_reps}{
  Number of repetitions to simulate over
}
  \item{progress}{
  Whether or not the function progress should print to screen. The default value is TRUE. If the user does not want the output to print to screen, they should set to FALSE.
}
\item{add}{
Default is FALSE. This is set to true if the user is running additional repetitions after the first call to ov_sim
}
\item{sim_archive}{
Default is NULL
}
}
\value{
ov_sim returns a list containing the following components:

\item{p_val}{matrix of pvalues for each grid point}
\item{trt_effect}{matrix of effect sizes for each grid point}
\item{es_grid}{vector of the effect size grid}
\item{rho_grid}{vector of the rho grid}
\item{cov}{vector of covariates used to estimate propensity score weights}
\item{data}{the initial data frame containing data with new weights}
\item{tx}{column name in data representing the treatment indicator}
\item{y}{column name in data representing the outcome}
\item{estimand}{estimand used}
\item{n_reps}{number of repetitions to simulate over}
\item{std.error}{matrix of standard errors for each grid point}
\item{es_se_raw}{matrix that stores each repetitions results at every grid point}
}

\examples{
data(sud)
sud = data.frame(sud)
sud$treat = ifelse(sud$treat == "A", 1, 0)
sud$wts = sample(seq(1, 10, by=.01), size=nrow(sud), replace = TRUE)
outcome_mod = outcome_model(data = sud,
                            weights = "wts",
                            treatment = "treat",
                            outcome = "eps7p_3",
                            model_covariates = c("sfs8p_0"),
                            estimand = "ATE")

ovtool_results = ov_sim(model_results=outcome_mod,
                        plot_covariates=c("sfs8p_0"),
                        es_grid = NULL,
                        rho_grid = NULL,
                        n_reps = 2,
                        progress=FALSE)
}


