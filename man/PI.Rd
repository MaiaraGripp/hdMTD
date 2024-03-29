% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PI.R
\name{PI}
\alias{PI}
\title{Estimated stationary distribution.}
\usage{
PI(S, freqTabSj, x_S, lenX, d)
}
\arguments{
\item{S}{A subset of 1:d}

\item{freqTabSj}{A table with frequencies of sequences in a MTD Markov chain sample.}

\item{x_S}{A specific sequence indexed by S.}

\item{lenX}{The length of the sample vector.}

\item{d}{The order of the chain.}
}
\value{
Returns a vector of estimated stationary distributions of a sequence x_S
with the element x_j ( for j in S)varying across the entire state space.
}
\description{
Estimated stationary distribution vector for each x_j and past x_S.
}
