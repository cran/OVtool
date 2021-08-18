\name{prep_for_plots}
\alias{prep_for_plots}
\title{prep_for_plots
}
\description{
Data preparation for producing the graphics and summary results.
}
\usage{
prep_for_plots(r1, p_contours)
}
\arguments{
  \item{r1}{
An object returned from ov_sim
}
  \item{p_contours}{
P-value countours to plot. The default plots: 0.01, 0.05, and 0.1. We only recommend changing this if the raw effect p-value is very close to one of these values. Do not specify more than four p-value contours.
}
}
\value{
prep_for_plots returns a list containing the following components:
\item{r1}{a list with the components returned from ov_simgrid}
\item{r1_df}{a data frame with components used to create the contour graphic}
\item{obs_cors}{a data frame with components used to plot the observed covariates on plot_graphic = "2" and plot_graphic = "3"}
\item{text_high}{a character noting the covariates whose absolute correlation with the outcome is greater than the grid allows}
\item{text_high_es}{a character noting the covariates with effect sizes greater than the maximum the plot will allow}
\item{pvals}{a vector of p-value thresholds to be plotted on the graphics}
\item{pval_lines}{a vector of line types to represent pvals}
\item{raw}{a character with the raw effect and pvalue from the outcome model}

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
}
