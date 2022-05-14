#' Calculates txy
#'
#' @param R A sx formed matrix
#' @param A_pairs Pairs in A
#'
#' @return The threshold txy=sx+sy
#'
txy <- function(R,A_pairs){
  tn <- matrix(0,nrow=nrow(R),ncol=nrow(A_pairs))
  for (s in 1:nrow(A_pairs)) {
    tn[,s] <- apply(R[,A_pairs[s,]],1,sum)
  }
  colnames(tn)=apply(A_pairs, 1, paste0,collapse="x")
  rownames(tn)=rownames(R)
  tn
  #tn=|s(1111)+s(2111)|s(1111)+s(3111)|s(2111)+s(3111)|
  #   |s(1211)+s(2211)|s(1211)+s(3211)|s(2211)+s(3211)|...
}
