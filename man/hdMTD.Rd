% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdMTD.R
\name{hdMTD}
\alias{hdMTD}
\title{Inference in MTD models}
\usage{
hdMTD(X, d, method = "FS", ...)
}
\arguments{
\item{X}{A vector or single-column data frame containing a chain sample.}

\item{d}{A positive integer representing an upper bound for the chain order.}

\item{method}{A method for estimation of the relevant lag set. The available methods in this
package are "FS" (default), "FSC", "CUT", and "BIC". Refer to the documentation of
each method's respective function for details.}

\item{...}{Additional arguments relevant to the selected method. Refer to the \emph{Details} section
for more information.}
}
\value{
Returns a vector with the estimated relevant lag set from a sample of an MTD model.
}
\description{
A function for inference in Mixture Transition Distribution (MTD) models. This function
can use a selected \code{method} to perform estimation of the relevant lag set from a sample of
an MTD model. By default \code{method="FS"} (Forward Stepwise) which is especially useful
in high dimension. The other available "methods" are "CUT", "FSC" (Forward Stepwise and Cut)
which is an application of the "FS" method followed by the "CUT" method, and lastly the "BIC"
(Bayesian Information Criterion) method. For more information on these methods see \emph{Details}
and the documentation of their specific functions.
}
\details{
This function allows the user to call any of the "\code{hdMTD_(method)}" functions
by specifying the "method" as an argument.
For example, if the \code{\link[=hdMTD]{hdMTD()}} function is used with \code{method="FSC"} it will call the
\code{\link[=hdMTD_FSC]{hdMTD_FSC()}} function. Note that, in this case, any extra parameters must match those used
by \code{\link[=hdMTD_FSC]{hdMTD_FSC()}}. Each method may use a different set of parameters, and they can be passed
to \code{\link[=hdMTD]{hdMTD()}} through the \code{...} argument. In order to see which parameters can be passed
for each method see the documentation of the "\code{hdMTD_(method)}" function:
\itemize{
\item For "FS" method, extra parameters are listed in the documentation of \code{\link[=hdMTD_FS]{hdMTD_FS()}}.
\item For "FSC" method, extra parameters are listed in the documentation of \code{\link[=hdMTD_FSC]{hdMTD_FSC()}}.
\item For "CUT" method, extra parameters are listed in the documentation of \code{\link[=hdMTD_CUT]{hdMTD_CUT()}}.
\item For "BIC" method, extra parameters are listed in the documentation of \code{\link[=hdMTD_BIC]{hdMTD_BIC()}}.
}
}
\examples{
X <- testChains[,1]
hdMTD(X=X,d=5, method = "FS",l=2)
hdMTD(X=X,d=5, method = "BIC",xi=1, minl=3, maxl=3)
}
