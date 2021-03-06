% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{adaptive_invnorm_ttest_os}
\alias{adaptive_invnorm_ttest_os}
\title{Inverse normal adaptive t-test}
\usage{
adaptive_invnorm_ttest_os(x, n1, n, ne, alpha = 0.025)
}
\arguments{
\item{x}{Observations}

\item{n1}{First stage sample size}

\item{n}{Pre-planned total sample size}

\item{ne}{Extended sample size}
}
\value{
p-value
}
\description{
Adaptive combination test of stage-wise t-tests using the inverse normal combination function.
}
\author{
Florian Klinglmueller
}

