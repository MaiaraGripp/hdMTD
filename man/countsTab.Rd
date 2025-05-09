% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countsTab.R
\name{countsTab}
\alias{countsTab}
\title{Counts sequences of length d+1 in a sample}
\usage{
countsTab(X, d)
}
\arguments{
\item{X}{A numeric vector, a single-column data frame, or a list with a
sample from a Markov chain. The first element must be the most recent observation.}

\item{d}{A positive integer specifying the number of elements in each sequence,
which will be \code{d+1}. Typically, \code{d} represents the chain order or
serves as an upper limit for it.}
}
\value{
A tibble with all observed sequences of length \code{d+1} and their
absolute frequencies.
}
\description{
Creates a tibble containing all unique sequences of length \code{d+1} found in
the sample, along with their absolute frequencies.
}
\details{
The function generates a tibble with \code{d+2} columns. In the first
\code{d+1} columns, each row displays a unique sequence of size \code{d+1}
observed in the sample. The last column, called \code{Nxa}, contains the
number of times each of these sequences appeared in the sample.

The number of rows in the output varies between \eqn{1} and \eqn{|A|^{d+1}},
where \eqn{|A|} is the number of unique states in \code{X}, since it depends
on the number of unique sequences that appear in the sample.
}
\examples{
countsTab(c(1,2,2,1,2,1,1,2,1,2), 3)

# Using test data.
countsTab(testChains[, 1], 2)

}
