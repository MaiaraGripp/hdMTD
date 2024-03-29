#' Inference in MTD models
#'
#' A function for inference in Mixture Transition Distribution (MTD) Markov chains. This function can use a selected "method" to perform estimation of the relevant lags set of a MTD chain sample.
#' The default method is "FS" (Forward Stepwise) which is specially useful in high dimension. The other available methods are "CUT", FSC" (Forward Stepwise and Cut) which is an application of the
#' "FS" method followed by the "CUT" method, and lastly the "BIC" (Bayesian Information Criterion) method. For more information on these methods see the documentation of their specific functions
#' listed on "details" below.
#'
#'
#' @details This function is simply a way to gather all of the "hdMTD_(method)" functions in a single place.
#' For example, if the [hdMTD()] function is used with method="FSC" it will call the [hdMTD_FSC()] function. Note that, in this case, any extra parameters must match those used by [hdMTD_FSC()].
#' Each method may use a different set of parameters, and they can be passed to [hdMTD()] through the ... argument. In order to see which parameters can be passed for each method
#' seek the documentation of the "hdMTD_(method)" function as follows:
#' \itemize{
#' \item For "FS" method, extra parameters are listed in the documentation of [hdMTD_FS()].
#' \item For "FSC" method, extra parameters are listed in the documentation of [hdMTD_FSC()].
#' \item For "CUT" method, extra parameters are listed in the documentation of [hdMTD_CUT()].
#' \item For "BIC" method, extra parameters are listed in the documentation of [hdMTD_BIC()].
#' }
#'
#' @param X A MTD chain sample.
#' @param d An upper bound for the chains order.
#' @param method A method for estimation of the relevant lag set. The available methods in this package are "FS" (default), "FSC", "CUT", and "BIC". Refer to the documentation of
#' each method's respective function for details.
#' @param ... Additional arguments relevant to the selected method. Refer to the "details" section for more information.
#'
#' @return Returns the estimated relevant lag set for a MTD chain sample.
#' @export
#'
#' @examples
#' X <- perfectSample(MTDmodel(c(1,5),c(0,1)),500)
#' hdMTD(X=X,d=5, method = "FS",l=2)
#' hdMTD(X=X,d=5, method = "BIC",xi=1, minl=3, maxl=3)
hdMTD <- function(X,d,method="FS", ...){

  if( !method %in% c("FSC","FS","CUT","BIC") ){stop("The chosen method is unknown")}
  fmtd <-  match.fun(paste0("hdMTD_",method))
  fmtd_params <- names(formals(fmtd))

  params <- list(...)
  if( method=="FS" | method=="FSC" ){
    while ( is.na(params$l) ||
            !is.numeric(params$l) ||
            params$l%%1 != 0 ||
            params$l>d ) {
      cat("The value of l is not valid or has not been declared for the FS step. It should be a positive integer less than or equal to d.")
      params$l <- readline(prompt = "Please enter a valid l: ")
      params$l <- suppressWarnings(as.numeric(params$l))
    }
  }

  dparams <- list(S=seq(1,d,1),l=d,alpha=0.05,mu=1,xi=0.5,minl=d,maxl=d,
                         A=NULL,byl=FALSE,BICvalue=FALSE,elbowTest=FALSE,warning=FALSE)

  if(length(params)!=0){
    if( !all(names(params) %in% fmtd_params) ){
      stop( paste0("Some parameters do not match those used in hdMTD_.",
                   method," function. Please check hdMTD_",method,"() documentation.") )
    }
    params_names <- names(params)
    dparams[params_names] <- params
  }
  fmtd(X=X,d=d,S=dparams$S,l=dparams$l,
       alpha=dparams$alpha,mu=dparams$mu,
       xi=dparams$xi,minl=dparams$minl,
       maxl=dparams$maxl,A=dparams$A,
       byl=dparams$byl,BICvalue=dparams$BICvalue,
       elbowTest=dparams$elbowTest,
       warning=dparams$warning)
}
