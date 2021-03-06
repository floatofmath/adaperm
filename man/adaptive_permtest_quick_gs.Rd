% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive_permtest.R
\name{adaptive_permtest_quick_gs}
\alias{adaptive_permtest_quick_gs}
\title{Permutation test for adaptive group sequential trials}
\usage{
adaptive_permtest_quick_gs(x1, x2, xE, g1, g2, gE, a0, stat,
  permutations = 10000, restricted = T, stratified = T)
}
\arguments{
\item{x1}{first stage observations}

\item{x2}{second stage observations}

\item{xE}{extended stage observations}

\item{g1}{frist stage treatment assignments}

\item{g2}{second stage treatment assignments}

\item{gE}{extended stage treatment assignments}

\item{a0}{early rejection boundary}

\item{stat}{test statistic to use}

\item{permutations}{number of permutations to use}

\item{restricted}{should the treatment group sizes be fixed}

\item{stratified}{}
}
\value{
pvalue
}
\description{
Quick no-nonsense permutation test for two-stage adaptive group-sequential trials that computes only the p-value
}
\details{
We assume that \code{stat} may also be used to compute the first stage statistic
}
\author{
Florian Klinglmueller
}

