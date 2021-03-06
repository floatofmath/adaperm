% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation.R
\name{p2t}
\alias{p2t}
\title{p to t}
\usage{
p2t(p, dist)
}
\arguments{
\item{p}{p-value/alpha level}

\item{dist}{permutation distribution}
}
\value{
critical value
}
\description{
Level alpha critical boundary of permutation test. Returns the smallest unique value of the permutation distribution such that the proportion of statistics larger is smaller than alpha.
}
\author{
Florian Klinglmueller
}

