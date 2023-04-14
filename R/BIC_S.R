#' Calculates de BIC given a set of relevant lags.
#'
#' Given a sample, this function estimates the BIC of a string whose set of relevant
#' lags is just the first element of S. Then it assumes that the set of relevant lags
#' are the first 2 elements of S and also estimates the BIC, and so on...
#'Lastly, estimate the BIC assuming that the set of relevant lags is exactly the set S.
#'
#'
#' @param X A Markov chain.
#' @param A The states space.
#' @param S A set of relevant lags, if empty S=\eqn{1,2,\dots, d}.
#' @param xi The BIC constant. Defaulted to 1/2. Smaller c `(near 0)` reduces the impact of overparameterization.
#'
#'@details This function assumes that the order of elements in S is determined by the relevance
#'(i.e. the first element of S is the most relevant lag, and probably the
#'one with the greatest oscillation, the last element of S is the least relevant,
#'probably has the smallest oscillation).
#'
#' @return  A vector of estimated BICs, the first entry supposes the first element
#' of S is the only relevant lag. The second entry supposes the first and the second
#' elements of S are the only relevant lags, and so on, until the penultimate entry
#' which supposes S is the relevant lag set. The last entry shows which of the previous
#' entries has the smallest BIC.
#' @export
#'
BIC_S <- function(X,A,S=1:d,xi=1/2){
  A <- sort(A)
  lenS <- length(S)
  base <- shapeSample(X,max(S))

  ML <- numeric(lenS)
  penalty <- numeric(lenS)
  pML <- numeric(lenS)
  for (i in seq_along(S)) {
    b <- base_Sja(S[1:i],j=NULL,A,base,complete = FALSE)
    ML[i] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
    penalty[i] <- n_parameters(1:i,A)*log(length(X))*xi
    pML[i] <- ML[i]+penalty[i]
    names(pML)[i] <- paste0(S[1:i],collapse = ",")
  }
  nam <- paste0("smallest: ",names(pML)[which( min(pML)==pML )])
  vec <- c(pML,min(pML))
  names(vec)[lenS+1] <- nam
  vec
}
