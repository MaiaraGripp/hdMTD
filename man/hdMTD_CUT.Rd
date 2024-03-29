% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdMTD_CUT.R
\name{hdMTD_CUT}
\alias{hdMTD_CUT}
\title{The Cut method.}
\usage{
hdMTD_CUT(
  X,
  d,
  S = 1:d,
  alpha = 0.05,
  mu = 1,
  xi = 0.5,
  A = NULL,
  warning = FALSE,
  ...
)
}
\arguments{
\item{X}{A MTD chain sample.}

\item{d}{An upper bound for the chains order.}

\item{S}{A subset of 1:d that contains the relevant lag set, if empty S=\eqn{1,2,\dots, d}.}

\item{alpha}{A parameter of the CUT algorithm. Alpha is constant of a threshold used in the CUT
step to determine if two distributions are different enough. The larger the alpha, the
greater the distance required between the distributions (to be considered different).}

\item{mu}{A parameter of the CUT algorithm. mu is also a component of the same threshold as alpha.}

\item{xi}{A parameter of the CUT algorithm. xi is also a component of the same threshold as alpha and mu.}

\item{A}{The states space. "A" only needs to be informed if X does not already contain all elements of "A".}

\item{warning}{If TRUE may return warnings.}

\item{...}{Used to accommodate any extra arguments passed by the \code{\link[=hdMTD]{hdMTD()}} function.}
}
\value{
Returns a estimated set of relevant lags.
}
\description{
A function for inference in MTD Markov chains with CUT method. It applies CUT algorithm to estimate a relevant lag set \eqn{\Lambda} of a MTD model.
}
\details{
The "Forward Stepwise and Cut" (FSC)is an algorithm for inference in
Mixture Transition Distribution (MTD) models.
It consists in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
This method was developed by \href{https://arxiv.org/abs/2202.08007}{Ost and Takahashi} and is specially useful for high-order MTD Markov chains.
This function will only apply the CUT step of the algorithm.
}
\examples{
X <- perfectSample(MTDmodel(Lambda=c(1,4),A=c(0,1)),N=1000)
hdMTD_CUT(X,4,alpha=0.02,mu=1,xi=0.4)
hdMTD_CUT(X,d=6,S=c(1,4,6),alpha=0.0065)

}
