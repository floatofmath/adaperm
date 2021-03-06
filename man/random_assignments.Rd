% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutation.R
\name{random_assignments}
\alias{random_assignments}
\title{Random combinations of two group assignments}
\usage{
random_assignments(n, B)
}
\arguments{
\item{n}{number of observations}

\item{B}{number of random combinations}
}
\value{
integer matrix of size \code{n} x \code{B}
}
\description{
Returns a logical matrix with B random assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control)
}
\author{
Florian Klinglmueller
}

