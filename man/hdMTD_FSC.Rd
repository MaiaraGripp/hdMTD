% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdMTD_FSC.R
\name{hdMTD_FSC}
\alias{hdMTD_FSC}
\title{Forward Stepwise and Cut method.}
\usage{
hdMTD_FSC(
  X,
  d,
  l,
  alpha = 0.05,
  mu = 1,
  xi = 0.5,
  A = NULL,
  warning = FALSE,
  ...
)
}
\arguments{
\item{X}{A MTD Markov chain sample.}

\item{d}{An upper bound for the chains order.}

\item{l}{Stop point for FS algorithm.}

\item{alpha}{A parameter of CUT step.}

\item{mu}{A parameter of CUT step.}

\item{xi}{A parameter of CUT step.}

\item{A}{The states space.}

\item{warning}{If TRUE may return warnings.}

\item{...}{Used to accommodate any extra arguments passed by the \code{\link[=hdMTD]{hdMTD()}}  function.}
}
\value{
Returns a estimated set of relevant lags.
}
\description{
A function for inference in MTD Markov chains with FSC method. This function estimates the relevant lag set \eqn{\Lambda}
of a MTD model through the FSC algorithm.
}
\details{
The "Forward Stepwise and Cut" (FSC)is an algorithm for inference in
Mixture Transition Distribution (MTD) models.
It consists in the application of the "Forward Stepwise" (FS) step followed by the CUT step.
This method was developed by \href{https://arxiv.org/abs/2202.08007}{Ost and Takahashi} and is specially useful for high order MTD Markov chains.
}
\examples{
X <- perfectSample(MTDmodel(Lambda=c(1,3),A=c(0,1)),N=2000)
hdMTD_FSC(X,4,3,alpha=0.02)
hdMTD_FSC(X,4,2,alpha=0.001)

}
