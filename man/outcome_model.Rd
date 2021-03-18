\name{outcome_model}
\alias{outcome_model}
\title{outcome_model
}
\description{
This function will run the outcomes model for your analysis. Upon completeion,
use the model object returned from this function and call ov_simgrid to check
the sensitivity of your findings.
}
\usage{
outcome_model(ps_object = NULL, stop.method=NULL, data, weights=NULL, treatment,
outcome, model_covariates, estimand = "ATE")
}
\arguments{
  \item{ps_object}{
  A ps object exported from TWANG
}
  \item{stop.method}{
  If the user specifies ps_object, stop.method should be used to export the weights (e.g "ks.max")
}
  \item{data}{
  A data frame containing the data
}
  \item{weights}{
  A column name in data that represents the relevant weights
}
  \item{treatment}{
  A column name in data for the treatment indicator
}
  \item{outcome}{
  A column name in data indicating the outcome vector
}
  \item{model_covariates}{
  A vector of column names representing the covariates in your final outcome's model
}
  \item{estimand}{
  "ATE" or "ATT"
}
}
\value{
outcome_model returns a list containing the following components:
\item{ps_object}{The ps_object from TWANG specified in the function call. If ignored, this component will be NULL}
\item{stop.method}{The stop method, if applicable, specified in the function call}
\item{data}{the updated data frame}
\item{weights}{the original vector of weights}
\item{tx}{a character name in data indicating the treatment indicator}
\item{y}{a character name in data indicating the outcome}
\item{outcome_mod_fmla}{the final outcome model formula}
\item{estimand}{The estimand specified in the function call}
\item{mod_results}{an object of class "svyglm"}

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
}

\references{
Lumley T (2020). “survey: analysis of complex survey samples.” R package version 4.0.
}

