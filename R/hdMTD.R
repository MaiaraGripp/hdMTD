#' Inference in MTD models
#'
#' A function for inference in mixture transition distribution (MTD) Markov chains. This function can use a selected "method" to perform estimation of the relevant lag set of a MTD chain sample.
#' The default method is "FS" (Foward Stepwise) which is specially useful in high dimension. The other available methods are "CUT", FSC" (Foward Stepwise and Cut) which is a application of the
#' "FS" method followed by the "CUT" method, and lastly the "BIC" (Bayesian Information Criterion) method. For more information on these methods see the documentation of their specific functions
#' listed on "details" below.
#'
#'
#' @details This function is simply a way to gather all of the hdMTD_ functions in a single place.
#' So, for example, if the [hdMTD()] function is used with method="FSC" it will call the [hdMTD_FSC()] function. Note that, in this case, any extra parameters must match those used by [hdMTD_FSC()].
#' Each method may use a different set of parameters, and they can be passed to [hdMTD()] trough the ... argument. In other to see with parameters can be passed for each method
#' seek the documentation of the hdMTD_"method" function as follows:
#' \itemize{
#' \item For "FS" method, extra parameters can be cheeked in [hdMTD_FS()] documentation.
#' \item For "FSC" method, extra parameters can be cheeked in [hdMTD_FSC()] documentation.
#' \item For "CUT" method, extra parameters can be cheeked in [hdMTD_CUT()] documentation.
#' \item For "BIC" method, extra parameters can be cheeked in [hdMTD_BIC()] documentation.
#' }
#'
#' @param X A mixture transition distribution (MTD) chain sample.
#' @param d An upper bound for the chains order.
#' @param method A method for estimation of the relevant lag set. The methods available in this package are "FS" (default),"FSC","CUT" and "BIC". See documentation for
#' each method with its respective function as explained in "details".
#' @param ... Extra arguments relevant to chosen method see "details" for more information.
#'
#' @return Returns the estimated relevant lag set for a MTD chain sample.
#' @export
#'
#' @examples
#' X <- perfectSample(MTDmodel(c(1,5),c(0,1)),5000)
#' hdMTD(X=X,d=5, method = "FS",l=2)
#' hdMTD(X=X,d=5, method = "FSC",alpha=0.001,xi=1,l=3)
#' hdMTD(X=X,d=5, method = "BIC",xi=1, minl=3, maxl=3)
#' hdMTD(X=X,d=5, method = "CUT",S=c(1,5,3,2),alpha=0.1,xi=1)
#'
hdMTD <- function(X,d,method="FS", ...){

  if( !method %in% c("FSC","FS","CUT","BIC") ){stop("The chosen method is unkown.")}
  fmtd <-  match.fun(paste0("hdMTD_",method))
  fmtd_params <- names(formals(fmtd))

  params <- list(...)
  if( method=="FS" | method=="FSC" ){
    while ( is.na(params$l) ||
            !is.numeric(params$l) ||
            params$l%%1 != 0 ||
            params$l>d ) {
      cat("l value is not valid or not declared for FS step. l should be a positive integer less than or equal to d.")
      params$l <- readline(prompt = "Please enter a valid l: ")
      params$l <- suppressWarnings(as.numeric(params$l))
    }
  }

  dparams <- list(S=seq(1,d,1),l=d,alpha=0.05,mu=1,xi=0.5,minl=d,maxl=d,
                         A=NULL,byl=FALSE,BICvalue=FALSE,elbowTest=FALSE,warning=FALSE)

  if(length(params)!=0){
    if( !all(names(params) %in% fmtd_params) ){
      stop( paste0("Some parameters do not match the used in hdMTD_",
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
