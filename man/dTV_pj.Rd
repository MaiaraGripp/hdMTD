% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dTV_pj.R
\name{dTV_pj}
\alias{dTV_pj}
\title{Total variation distance for p_j distribution.}
\usage{
dTV_pj(p_j, rows, which_max = FALSE)
}
\arguments{
\item{p_j}{A stochastic matrix.}

\item{rows}{A matrix with two columns, the pair of numbers in each row represents the rows of p_j that
will be compared, that is, the distributions in p_j to be compared.}

\item{which_max}{If TRUE, the function will also return the rows of the p_j matrices with the distributions
that have the greatest total variation distance.}
}
\value{
Returns the total variation distance of a set of distributions in p_j.
}
\description{
Calculates the total variation distance given a stochastic matrix p_j.
}
