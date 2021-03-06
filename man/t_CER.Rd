% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive.R
\name{t_CER}
\alias{t_CER}
\title{t - Test Conditional Error Rate}
\usage{
t_CER(x1, g1, n, alpha = 0.025, one_sample = FALSE, ...)
}
\arguments{
\item{x1}{First stage observations}

\item{g1}{First stage treatment assignments (1 for control, -1 for treatment)}

\item{n}{Preplanned overal sample size}

\item{alpha}{pre-fixed significance level}
}
\value{
numeric value of the conditional error rate
}
\description{
Computes the conditional type I error rate of a pre-planned t-test (using inverse normal combination of stage-wise t-tests) in a two-stage parallel group adaptive design. We condition on the observed first stage data and treatment assignment as well as the pre-planned second stage sample size.
}
\author{
Florian Klinglmueller
}

