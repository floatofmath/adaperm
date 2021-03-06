% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{cond_power_rule_u_ts}
\alias{cond_power_rule_u_ts}
\title{Conditional power sample size reassessment rule (Mann-Whitney U Test)}
\usage{
cond_power_rule_u_ts(x1, y1, target = 0.9, alpha = 0.025,
  maxN = length(x1) * 6)
}
\arguments{
\item{x1}{First stage control group observations}

\item{y1}{First stage treatment group observations}

\item{target}{Target power}

\item{alpha}{Significance level}

\item{maxN}{Maximum sample size}
}
\value{
Total sample size (in the control group) required to achieve the target power
}
\description{
Conditional power rule for the Mann-Whitney U test, where we estimate the relative effect size based on the observed first stage outcomes by simply counting the number of treatment group samples that have larger outcome values than control group samples.
}
\author{
float
}

