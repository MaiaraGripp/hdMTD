#' Bayesian Information Criterion `(BIC)` of a Markov Chain.
#'
#' @param X Markov chain.
#' @param A The state space.
#' @param d A upper threshold for the chains order.
#' @param l A upper threshold for the number of elements in the relevant lag set.
#' @param c The BIC constant. Defaulted to 1/2. Smaller c `(near 0)` reduces the impact of overparameterization.
#'
#' @return The BIC of a Markov Chain for each possible set of lags with size 1,2...,l.
#' @importFrom utils combn
#' @export
BIC_l <- function(X,A,d,l,c=1/2){ #smaller c gives less weight to param number
  #tests
  if( !is.numeric(l) ||
      length(l)!=1   ||
      l%%1!=0        ||
      l<1 )stop("l must be an integer number greater than 0 and smaller or equal to d.")
  if( !is.numeric(d) ||
      length(d)!=1   ||
      d%%1!=0        ||
      d<l )stop("d must be an integer number greater or equal to l.")
  if( !is.numeric(c) ||
      c<=0           ||
      length(c)!=1 )stop("c must be a numeric positive constant.")
  if( length(A)<=1 ||
      length(dim(A))!=0 )stop("A must be a vector with at least two values.")
  #\

  base <- shapeSample(X,d)
  if(l==1){
    tryCombs <- matrix( c( rep(0,d), rep(n_parameters(1,A)*log(length(X))*c,d) ),
                        byrow = T,nrow = 2 )
    colnames(tryCombs) <- 1:d
    rownames(tryCombs) <- c("log_MV","penalty")
    for (k in 1:d) {
      S <- k
      b <- base_Sja(S,j=NULL,A,base,complete = FALSE)
      tryCombs[1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
    }
    penalizedML <- apply(tryCombs,2,sum)

  }else{#if l>1
    tryCombs <- list()
    for (i in 1:l) {
      ncombs <- choose(d,i)
      tryCombs[[i]] <- matrix( rep(0,2*ncombs) ,byrow = T,nrow = 2 )
      tryCombs[[i]][2,] <- n_parameters(Lambda=(1:i),A)*log(length(X))*c
      aux <- apply(t(combn(1:d,i)), 1, paste0,collapse=",")
      for ( k in 1:ncombs ) {
        S <- as.numeric(unlist(strsplit(aux[k], ",")))
        b <- base_Sja(S,j=NULL,A,base,complete = FALSE)
        tryCombs[[i]][1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
      }
      colnames(tryCombs[[i]]) <- aux
      rownames(tryCombs[[i]]) <- c("log_MV","penalty")
    }
    penalizedML <- sapply(tryCombs, apply, 2,sum)
  }
  penalizedML
}
