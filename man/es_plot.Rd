\name{es_plot}
\alias{es_plot}
\title{es_plot
}
\description{
Plot of the effect size contours
}
\usage{
es_plot(prep, col="bw")
}
\arguments{
  \item{prep}{
Input from prep_for_plots
}
  \item{col}{
If user wants color (a heat map layered with contours) or black and white (contours only). Specify the heat map with "color" or black and white "bw".
}
}
\value{
a list of class gg and ggplot
}
\examples{
data(sud)
sud = data.frame(sud[sample(1:nrow(sud),100),])
sud$treat = ifelse(sud$treat == "A", 1, 0)
sud$wts = sample(seq(0.01, .1,by=.01), size=nrow(sud), replace = TRUE)
outcome_mod = outcome_model(data = sud,
                            weights = "wts",
                            treatment = "treat",
                            outcome = "eps7p_6",
                            model_covariates = c("sfs8p_0", "eps7p_0",
                                                 "ada_0"),
                            estimand = "ATE")

ovtool_results = ov_sim(model_results=outcome_mod,
                        plot_covariates=c("sfs8p_0", "ada_0"),
                        es_grid = 0,
                        rho_grid = 0,
                        n_reps = 2,
                        progress=FALSE)
prep = prep_for_plots(ovtool_results, p_contours=.05)
plot = es_plot(prep = prep)
}
