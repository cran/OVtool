\name{gen_a_finish}
\alias{gen_a_finish}
\title{gen_a_finish
}
\description{
This function will generate a at finish.
}
\usage{
gen_a_finish(a_res, my_estimand)
}
\arguments{
  \item{a_res}{
A list of values returned by gen_a_start
}
  \item{my_estimand}{
"ATE" or "ATT"
}
}
\value{
\item{a}{used to control the strength of the relationship between the omitted variable and the treatment}
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

finish = gen_a_finish(a_res = start, my_estimand = "ATE")
}
