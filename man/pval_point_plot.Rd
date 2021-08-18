\name{add_pvals_plot}
\alias{add_pvals_plot}
\title{add_pvals_plot
}
\description{
Plot of effect size contours with pvalue overlay
}
\usage{
add_pvals_plot(prep, col = "color")
}
\arguments{
  \item{prep}{
Input from prep_for_plots
}
  \item{col}{
If user wants color or black and white. Specify color with "color" or black and white "bw"
}
}
\value{
a list of class gg and ggplot
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
plot = add_pvals_plot(prep = prep)
}
