\name{plot.ov}
\alias{plot.ov}
\title{plot
}
\description{
Plots the user specified graphic(s)
}
\usage{
\method{plot}{ov}(x, col="color", print_graphic="1", p_contours = c(0.01, 0.05, 0.1), \dots)
}

\arguments{
  \item{x}{
Object returned from the call to ov_sim
}
  \item{col}{
If user wants color or black and white. Specify color with "color" or black and white "bw"
}
  \item{print_graphic}{
Takes values "1", "2", or "3", depending what graphics the user wants
}
  \item{p_contours}{
P-value countours to plot. The default plots: 0.01, 0.05, and 0.1. We only recommend changing this if the raw effect p-value is very close to one of these values. Do not specify more than four p-value contours.
}
\item{...}{
Additional arguments.
}
}

\value{
This function will print the plot to screen that the use specifies with print_graphic.
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
plot = plot.ov(ovtool_results, print_graphic="3", p_contours=.05)
}
