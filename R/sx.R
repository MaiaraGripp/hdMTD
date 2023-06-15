#' Calculates sx, and txy=sx+sy
#'
#' @param S A set of relevant lags.
#' @param freqTab A table with frequencies of a Markov chain.
#' @param lenA The size of the states space.
#' @param x_S A sequence indexed by S.
#' @param mu A parameter of CUT.
#' @param alpha A parameter of CUT.
#' @param xi A parameter of CUT.
#'
#' @return Returns sx the for the threshold txy=sx+sy
#' @importFrom dplyr %>%
#'
sx <- function(S,freqTab,lenA,x_S,mu,alpha,xi){
  filtr_S <- paste0("x",S)
  B <- freqTab
  B$test <- apply(B %>%
                    dplyr::select_at(filtr_S),1,is_xS,x_S)
  C <- dplyr::filter(B,test==TRUE)
  Nx <- C$Nx_Sj[seq(1,nrow(C),by=lenA)]
  result <- numeric(lenA)
  for (k in 1:lenA) {
    sums <- 0
    for (i in 1:lenA) {
      sums = sums +
                sqrt(( C$qax_Sj[i+(k-1)*lenA] +
                alpha/Nx[k])*mu/(2*mu-exp(mu)+1))
    }
    result[k]=sqrt(2*alpha*(1+xi)/Nx[k])*sums+alpha*lenA/(3*Nx[k])
  }
  result
}

