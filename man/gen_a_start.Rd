\name{gen_a_start}
\alias{gen_a_start}
\title{gen_a_start
}
\description{
This function is a wrapper to ov_simgrid. It generates the a. a is used to control the strength of the relationshp between the unobserved counfounder, U, and the treatment indicator
}
\usage{
gen_a_start(y, tx, residuals, es, rho, my_estimand)
}
\arguments{
  \item{y}{
A vector that represents the outcome.
}
  \item{tx}{
A vector for the treatment indicator (must be 0s and 1s).
}
  \item{residuals}{
A vector of residuals from regressing Y on X and controlling for treatment.
}
  \item{es}{
An effect size value to simulate over.
}
  \item{rho}{
A rho (correlation) value to simulate over.
}
  \item{my_estimand}{
"ATE" or "ATT"
}
}

\value{
gen_a_start returns a list containing the following components:
\item{n1}{scalar representing sample size of treatment group (treat == 1)}
\item{ve1}{1 - b1^2 multiplied by the variance of Ystar1}
\item{b1}{bounded parameter for treatment group (treat == 1) that it with b0 are selected to set the correlation of the omitted variable and the outcome equal to rho}
\item{es}{}
\item{pi}{proportion of population that is in the treatment group (treat == 1)}
\item{n0}{scalar represnting sample size of control group (treat == 0)}
\item{ve0}{1 - b0^2 multiplied by the variance of Ystar0}
\item{b0}{bounded parameter for control group (treat == 0) that it with b1 are selected to set the correlation of the omitted variable and the outcome equal to rho}
\item{n}{scalar representing the total sample size}
\item{ind}{vector of positions in data that represent treatment group (treat == 1)}
\item{Rstar_1}{Residuals in treatment group}
\item{Rstar_0}{Residuals in control group}
}
\examples{
data(sud)
sud = data.frame(sud)
sud$treat = ifelse(sud$treat == "A", 1, 0)
sud$wts = sample(seq(0.01, 1,by=.01), size=nrow(sud), replace = TRUE)
outcome_mod = outcome_model(data = sud,
                            weights = "wts",
                            treatment = "treat",
                            outcome = "eps7p_3",
                            model_covariates = c("sfs8p_0"),
                            estimand = "ATE")

start = gen_a_start(y=sud$eps7p_3, tx=sud$treat,
                    residuals=residuals(outcome_mod$mod_results),
                    es = .01,
                    rho = .01,
                    my_estimand = "ATE")
}

