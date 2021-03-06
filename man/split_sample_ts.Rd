% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{split_sample_ts}
\alias{split_sample_ts}
\title{Split two-group sample}
\usage{
split_sample_ts(x, y, n1, n, m1, m)
}
\arguments{
\item{x}{Observations control group}

\item{y}{Observations treatment group}

\item{n1}{First stage sample size}

\item{n}{Pre-planned total sample size}

\item{m1}{First stage sample size}

\item{m}{Pre-planned total sample size}
}
\value{
list with stagewise observations and stagewise sign-indicators
}
\description{
Split sample of a two group experiment into stages
}
\author{
Florian Klinglmueller
}

