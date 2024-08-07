#' A tibble with probabilities and frequencies of sample sequences
#'
#' A tibble with the sample sequences, their frequencies and the estimated
#' transition probabilities.
#'
#' @details Within this function, both arguments \code{S} and \code{j} perform the same role of
#' filtering the columns of the \code{countsTab} argument to be kept in the \code{freqTab()} output.
#' In fact, any specific lag \code{j} could be provided to this function via the \code{S} argument,
#' while leaving the \code{j} argument as \code{NULL} as defaulted. The output would be exactly the
#' same as when specifying \code{S} and \code{j} separately. The inclusion of \code{j} as a
#' parameter is merely for organizational purposes and is a useful tool for clarity within
#' the package's algorithms.
#'
#' @param S A numeric vector of positive integers or \code{NULL}. Represents a set
#'  of past lags that must be present within the columns of the \code{countsTab} argument and are
#'  to be considered while estimating the transition probabilities. Both \code{S} and \code{j}
#'  cannot be \code{NULL} at the same time.
#' @param j An integer or \code{NULL}. Typically represents a lag \code{j} in the
#'  \eqn{complement} of \code{S}. Both \code{S} and \code{j} cannot be \code{NULL} at the same time.
#'  See *Details* for further information.
#' @param A A vector with positive integers. Must have at least two different entries. \code{A}
#' represents the state space.
#' @param countsTab A tibble or a data frame with all sequences of length d+1 that appear in the sample,
#'  and their absolute frequency. Such tibble typically is generated by the function [countsTab()].
#'  If any data frames that were not generated by [countsTab()] are used, problems may occur since
#'  its format and column names are important within [freqTab()] function.
#' @param complete Logical. If \code{TRUE} all sequences that didn't appear in the sample
#' will be included in the output with frequency equal to 0 .
#'
#' @return A tibble where each row represents a sequence of elements from \code{A}.
#' The initial columns display each sequence symbol separated into columns corresponding to their time indexes.
#' The remaining columns show the sample frequencies of the sequences and the MLE (Maximum Likelihood Estimator)
#' of the transition probabilities.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom methods is
#' @examples
#' freqTab(S=c(1,4),j=2,A=c(1,2,3),countsTab = countsTab(testChains[,2],d=5))
#' #Equivalent to freqTab(S=c(1,2,4),j=NULL,A=c(1,2,3),countsTab = countsTab(testChains[,2],d=5))
freqTab <- function(S,j=NULL,A,countsTab,complete=TRUE){
  # Cheking inputs
  if(!is.data.frame(countsTab)){stop("countsTab must be a tibble or a dataframe. Try using countsTab() function.")}
  d <- ncol(countsTab)-2
  if(length(S)!=0){
    if(!is.vector(S) || any(S%%1 != 0) || any(S<1) || length(S)<0 )stop("S must be a positive integer vector, a number or NULL.")
  }
  if(length(j)!=0){
    if(length(j)!=1 || j%%1!=0 || j %in% S ){stop("j must be a integer number in the complement of S.")}
  }
  Sj <- sort(c(S,j),decreasing = TRUE)
  if(length(Sj)==0){stop("The set {S}U{j} can't be NULL.")}
  if(Sj[1] > d){stop("The set {S}U{j} can't have an element greater than d which is the upper bound for the chains order.")}
  if(!is.logical(complete)){stop("Complete must be a logical argument.")}
  if( !all(unique(unlist(countsTab[,-(d+2)])) %in% A) ){stop("A must contain all elements that appear in the countsTab sequencies.")}
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")

  A <- sort(A)
  lenSj <- length(Sj)
  filtrs <- c(paste0("x",Sj),"a")
  lenA <- length(A)

## summarising countsTab
  freqTab <- countsTab %>%
                dplyr::group_by_at(filtrs) %>%
                dplyr::summarise(Nxa_Sj=sum(Nxa), .groups="drop")
## If there are sequences that didn't appear in sample and complete=TRUE,
# adds thoses sequences with frequency 0.
  if ( ( nrow(freqTab) < lenA^(lenSj+1) ) && complete ){
## However adding those sequences might create too big a dataframe, so it must test for that.
    Tablexa <- try(expand.grid(rep(list(A),lenSj+1))[,order((lenSj+1):1)],silent = TRUE)
    if(is(Tablexa,"try-error")){stop(paste0("The dataset with all sequences of size length(S) is too large."))}

    list1 <- apply( freqTab[,1:(lenSj+1)],1,paste0,collapse="" )
    list2 <- apply( Tablexa, 1, paste0,collapse="" )
    Tablexa <- Tablexa[ match(setdiff(list2,list1),list2),]
    Tablexa <- cbind(Tablexa,0)
    colnames(Tablexa) <- colnames(freqTab)
    freqTab <- rbind(freqTab,Tablexa)
    freqTab <- dplyr::arrange_at(freqTab, filtrs)
  }

  freqTab <- freqTab %>% dplyr::group_by_at(paste0("x",Sj)) %>%
                dplyr::mutate(Nx_Sj=sum(Nxa_Sj))
  freqTab <- freqTab %>%
                dplyr::mutate(qax_Sj=dplyr::if_else(Nx_Sj>0,Nxa_Sj/Nx_Sj,1/lenA)) %>%
                   dplyr::ungroup()
  dplyr::as_tibble(freqTab)
}
