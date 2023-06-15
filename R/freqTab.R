#' Table of estimated transition probabilities and frequencies of sample sequences.
#'
#' A table with the sample sequences (including time 0) frequencies , and the estimated transition probabilities.
#'
#' @param S A subset of \eqn{1:d}.
#' @param j An element j in set \eqn{ {1:d} \ {S} }.
#' @param A The states space.
#' @param countsTab A data frame with the sequences of size d+1 in sample and their absolute frequency. The names of the
#' columns in countsTab must be \eqn{ {"x_{-j}" : j in d to 1 }  } then "a" and "Nxa". The use of function [countsTab()]
#' is recommended for generating this data frame.
#' @param complete logical, if TRUE and if there are sequences that didn't appear in the sample,
#'the freqTab will be completed with such sequences appearing 0 times.
#'
#' @return A table with the sequences that appeared in the sample (at times \eqn{ {-S}U{-j}U{0} } for each size d+1 sequence),
#' their frequencies and the MLE estimator of the transition probabilities.
#' @export
#' @importFrom dplyr %>%
#' @examples freqTab(S=c(1,3),j=NULL,A=c(0,1,2),
#' countsTab=countsTab(c(1,0,1,0,2,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0),d=3))
#' freqTab(S=c(1,3),j=NULL,A=c(0,1,2),
#' countsTab=countsTab(c(1,0,1,0,2,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0),d=3),complete=FALSE)
freqTab <- function(S,j=NULL,A,countsTab,complete=TRUE){
  # Cheking inputs
  if(!is.data.frame(countsTab)){stop("countsTab must be a dataframe. Try usings countsTab function.")}
  d <- ncol(countsTab)-2
  if(length(S)!=0){
    if(length(S) < 1  ||
       length(dim(S))!=0 ||
       any( S%%1 != 0 ) ){stop("S must be a vector of at least 1 integer number.")}
  }
  if(length(j)!=0){
    if(length(j)!=1 ||
       j%%1!=0      ||
       j %in% S ){stop("j must be a integer number in the complement of S in (1:d).")}
  }
  Sj <- sort(c(S,j),decreasing = TRUE)
  if(length(Sj)==0){stop("The set {S}U{j} can't be NULL.")}
  if(Sj[1] > d){stop("The set {S}U{j} can't have an element greater than d which is the upper bound for the chains order.")}
  if(!is.logical(complete)){stop("Complete must be a logical argument.")}
  if( !all(unique(unlist(countsTab[,-(d+2)])) %in% A) ){stop("A must contain all elements that appear in the countsTab sequencies.")}
  A <- sort(A)
  lenSj <- length(Sj)
  filtrs <- c(paste0("x",Sj),"a")
  lenA <- length(A)

  freqTab <- countsTab %>%
                dplyr::group_by_at(filtrs) %>%
                dplyr::summarise(Nxa_Sj=sum(Nxa), .groups="drop")

  if ( ( nrow(freqTab) < lenA^(lenSj+1) ) && complete ){
    Tablexa <- expand.grid(rep(list(A),lenSj+1))[,order((lenSj+1):1)]
    list1 <- apply( freqTab[,1:(lenSj+1)],1,paste0,collapse="" ) #list of sequences in freqTab
    list2 <- apply( Tablexa, 1, paste0,collapse="" ) #list of sequences in Tablexa
    Tablexa <- Tablexa[ match(setdiff(list2,list1),list2),] #list of sequences in Tablexa that where not in freqTab
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
  freqTab
}
