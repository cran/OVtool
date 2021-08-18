\name{find_esgrid}
\alias{find_esgrid}
\title{find_esgrid
}
\description{
Finds a reasonable effect size grid to simulate over.
}
\usage{
find_esgrid(my_data, my_cov, treatment, outcome, my_estimand)
}
\arguments{
  \item{my_data}{
Data
}
  \item{my_cov}{
vector of covariates
}
  \item{treatment}{
column name of treatment indicator in my_data
}
  \item{outcome}{
column name of outcome in my_data
}
  \item{my_estimand}{
Relevant estimand ("ATE" or "ATT")
}
}
\value{
a data frame with three columns, "Cor_Outcome", "es", and "cov". "Cor_Outcome" represents rho grid values, "ES" represents the range of grid values to represent the association between the unobserved confounder and the treatment indicator on the effect size scale, and "cov" is a vector of all the covariates used in the propensity score model
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

find_es = find_esgrid(sud, my_cov="sfs8p_0", treatment="treat",
                      outcome="eps7p_3",
                      my_estimand="ATE")
}

