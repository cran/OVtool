\name{sud}
\alias{sud}
\docType{data}
\title{
Longitudinal observational data from adolescents receiving SUD treatment.
}
\description{
A dataset containing substance use disorder and mental health measures for adolescents who had one of two substance use treatments.
}

\usage{data("sud")}
\format{
  A data frame with 4000 observations on the following 29 variables.
  \describe{
  \item{\code{treat}}{treatment indicator}
  \item{\code{tss_0}}{Traumatic Stress Scale, baseline}
  \item{\code{tss_3}}{Traumatic Stress Scale, recorded at 3-months}
  \item{\code{tss_6}}{Traumatic Stress Scale, recorded at 6-months}
  \item{\code{sfs8p_0}}{Substance Frequency Scale, baseline}
  \item{\code{sfs8p_3}}{Substance Frequency Scale, recorded at 3-months}
  \item{\code{sfs8p_6}}{Substance Frequency Scale, recorded at 6-months}
  \item{\code{eps7p_0}}{Emotional Problems Scale, baseline}
  \item{\code{eps7p_3}}{Emotional Problems Scale, recorded at 3-months}
  \item{\code{eps7p_6}}{Emotional Problems Scale, recorded at 6-months}
  \item{\code{ias5p_0}}{Illegal Activity Scale, baseline}
  \item{\code{dss9_0}}{Depressive Symptom Scale-9 Item, baseline}
  \item{\code{mhtrt_0}}{MH treatment, past 90 days, baseline}
  \item{\code{sati_0}}{Substance Abuse Tx Index, baseline}
  \item{\code{sp_sm_0}}{Substance Problem Scale, Past Month, baseline}
  \item{\code{sp_sm_3}}{Substance Problem Scale, Past Month, recorded at 3-months}
  \item{\code{sp_sm_6}}{Substance Problem Scale, Past Month, recorded at 6-months}
  \item{\code{gvs}}{General Victimization Scale}
  \item{\code{ers21_0}}{Environment Risk Scale, baseline}
  \item{\code{nproc}}{Count of Treatment A procedures delivered to client}
  \item{\code{ada_0}}{Adjusted Days Abstinent-Any, baseline}
  \item{\code{ada_3}}{Adjusted Days Abstinent-Any, recorded at 3-months}
  \item{\code{ada_6}}{Adjusted Days Abstinent-Any, recorded at 6-months}
  \item{\code{recov_0}}{Binary indicator indicating if in recovery, baseline}
  \item{\code{recov_3}}{Binary indicator indicating if in recovery, recorded at 3-months}
  \item{\code{recov_6}}{Binary indicator indicating if in recovery, recorded at 6-months}
  \item{\code{subsgrps_n}}{Categorical variable where: 1="Alcohol and/or marijuana disorder/weekly use; 2="Other drugs"; 3="Opiate disorder/weekly use"}
  \item{\code{sncnt}}{Total number of sessions for Treatment A}
  \item{\code{engage}}{Binary indicator indicating initiated treatment and had 4+ sesssions within 45 days for Treatment A}
  }
}
\source{
Global Appraisal of Individual Needs biopsychosocial assessment instrument - GAIN - Dennis, Titus et al. 2003
}
\examples{
data(sud)
}
\keyword{datasets}
