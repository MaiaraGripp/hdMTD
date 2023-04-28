#' A function for inference in MTD Markov chains with BIC.
#'
#' Estimates a relevant lag set \eqn{\Lambda} of MTD models using BIC.
#'
#' @param X A Markov chain.
#' @param A The states space.
#' @param S A set of all possible relevant lags, if only an upper threshold (d) for the chains order is known, set \eqn{S=1:d}. If S is a subset
#'  of the set 1:d, the order of the chain will be considered the greatest element of S.
#' @param l If \eqn{lset=FALSE} l is an upper bound for the number of elements to be returned in the estimated relevant lag set.
#'  If \eqn{lset=TRUE} l will be the exact number number of elements to be returned in the estimated relevant lag set.
#' l is a
#' @param xi The BIC constant. Defaulted to 1/2. Smaller xi `(near 0)` reduces the impact of overparameterization.
#' @param nested If TRUE the function will interpret that the elements in S are sorted from most relevant to least relevant and will calculate BIC supposing \eqn{\Lambda=S[1]}, then \eqn{\Lambda=S[1:2], \dots, \Lambda=S[1:l]}.
#' @param smallestBIC If TRUE will return only the set with smallest BIC between all possible relevant lag sets with elements of S and sizes from 1 til l.
#'if \eqn{nested=TRUE} will suppose the elements of S are sorted from most relevant to least relevant and will return only the smallest BIc between the models supposing
#'\eqn{\Lambda=S[1]}, \eqn{\Lambda=S[1:2], \dots \Lambda=S[1:l]}. Can't be TRUE if \eqn{lset=TRUE}.
#' @param lset If TRUE will return the smallest BIC between the models with exactly l elements of S. If \eqn{nested=TRUE} will return only the BIC of the model with \eqn{\Lambda=S[1:l]}.
#' Can't be TRUE if \eqn{smallestBIC=TRUE}
#' @param warning If TRUE may return warnings.
#'
#' @details smallestBIC and lset can't both be TRUE at the same time since, \eqn{lset=TRUE} looks for the smallest BIC between the all combinations of size l sets. While \eqn{smallestBIC=TRUE}
#' looks for the smallest BIC between all combinations of sizes 1 to l sets.
#' @details Note that the order of the chain matters at the counting of sequences step. If we run the function for \eqn{S=c(1,3,6)} the order d will be set to 6, while if we run for \eqn{S=c(1,3,5)}
#' the order will be \eqn{d=5}. That will alter the count of sequences that appear in the sample and might alter BIC values for the estimated lag sets. Meaning that the same output of relevant lags
#' (for example \eqn{1,3}) can have different BIC values if there was an initial order difference.
#'
#' @return Returns estimations, using Bayesian Information Criterion (BIC), of the relevant lag set of size \eqn{1,2, \dots l}.
#'  For example, suppose \eqn{S=c(9,5,1,2)}, \eqn{l=3}. This function
#' will run from \eqn{l=1 to l=3} and calculate the BIC supposing the MTD model has \eqn{\Lambda=c(1)}, then will calculate the BIC for \eqn{\Lambda=c(2)}, then \eqn{\Lambda=c(5)}, and finally \eqn{\Lambda=c(9)}.
#' Then it will suppose \eqn{\Lambda=c(1,2)}, \eqn{\Lambda=c(1,5)}, ... , \eqn{\Lambda=c(5,9)}, and then \eqn{\Lambda=c(1,2,5)} til \eqn{\Lambda=c(2,5,9)}. After calculating all BIC for each model will return the
#' sets with smallest BIC by each l size. And finally the set with smallest BIC between all l sizes. If \eqn{nested=TRUE} the function will suppose the elements of S are sorted from most relevant to least, and will only
#' calculate the relevant lags supposing \eqn{\Lambda=9}, then \eqn{\Lambda=c(9,5)}, \eqn{\Lambda=c(9,5,1)} stopping here since \eqn{l=3}.
#' @export
#'
#' @examples
#' X <- perfectSample(MTDmodel(Lambda=c(1,3,6),A=c(0,1,2)),N=2000)
#' sparseMarkov_BIC(X,S=1:7,l=4)
#' sparseMarkov_BIC(X,S=1:7,l=4,xi=0.8)
#' sparseMarkov_BIC(X,S=c(6,3,1,4,7),l=3)
#' sparseMarkov_BIC(X,S=c(6,3,1,4,7),l=4,nested=TRUE)
#' sparseMarkov_BIC(X,S=c(6,3,1,4,7),l=4,lset=TRUE)
#' sparseMarkov_BIC(X,S=c(6,3,1,4,7),l=4,xi=0.75,smallestBIC=TRUE)
#' sparseMarkov_BIC(X,S=c(6,3,1,4,7),l=2,lset=TRUE, nested = TRUE)
#'
sparseMarkov_BIC <- function(X,A=NULL,S,l=max(S),xi=1/2,nested=FALSE,
                             smallestBIC=FALSE,lset=FALSE, warning=FALSE){
  #Checking inputs
  X <- checkSample(X)
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A not informed. Code will set A <- sort(unique(X)).")
    }
    A <- unique(X)
  }else{
    if( !is.numeric(A) ||
        length(A)<=1   ||
        length(dim(A))!=0 )stop("States space A must be a numeric vector with at least two values.")
  }
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")}
  #if( !is.numeric(d) || d<2 || (d %% 1)!=0  ||  d<max(S)){
  #stop("The order d must be an integer number greater than 2 or the greatest elements in S.")
  #}
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>length(S) ) {
    cat("l value is not valid. l should be a positive integer lower or equal to the number of elements in S.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("BIC constant xi value is not valid. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  if(!is.logical(nested)){stop("nested must me a logical argument")}
  if(!is.logical(warning)){stop("warning must me a logical argument")}
  if(!is.logical(smallestBIC)){stop("smallestBIC must me a logical argument")}
  if(!is.logical(lset)){stop("lset must me a logical argument")}
  if(smallestBIC & lset){stop("Arguments smallestBIC and lset can't both be TRUE at the same time")}

  #gatherin inputs
  A <- sort(A)
  lenS <- length(S)
  d <- max(S)
  base <- shapeSample(X,d)
  #\.

  if(!nested){
    S <- sort(S)

    if(!lset){

      if(l==1){
        tryCombs <- matrix( c( rep(0,lenS), rep(n_parameters(1,A)*log(length(X))*xi,lenS) ),
                            byrow = T,nrow = 2 )
        colnames(tryCombs) <- S
        rownames(tryCombs) <- c("log_ML","penalty")
        for (k in S) {
          b <- base_Sja(k,j=NULL,A,base,complete = FALSE)
          tryCombs[1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
        }
        pML <- apply(tryCombs,2,sum)
      }else{ # if l>1
        tryCombs <- list()
        for (i in 1:l) {
          nCombs <- choose(lenS,i)
          tryCombs[[i]] <- matrix( rep(0,2*nCombs) ,byrow = T,nrow = 2 )
          tryCombs[[i]][2,] <- n_parameters(Lambda=(1:i),A)*log(length(X))*xi
          aux <- apply(t(combn(S,i)), 1, paste0,collapse=",")
          for ( k in 1:nCombs ) {
            G <- as.numeric(unlist(strsplit(aux[k], ",")))
            b <- base_Sja(G,j=NULL,A,base,complete = FALSE)
            tryCombs[[i]][1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
          }
          colnames(tryCombs[[i]]) <- aux
          rownames(tryCombs[[i]]) <- c("log_ML","penalty")
        }
        pML <- sapply(tryCombs, apply, 2,sum)
      }
      smallest <- names(unlist(pML))[order(unlist(pML))][1]
      nam <- sapply(sapply(pML, orderedNames), dplyr::first )
      pML <- sapply(sapply(pML,sort),dplyr::first)
      pML <- c( pML,unlist(pML)[order(unlist(pML))][1] )
      names(pML) <- c(nam,paste0("smallest: ",smallest)) #nested=FAlSE & lset=FALSE

    }else{ #if lset=TRUE

      combs <- combn(S,l)
      nCombs <- choose(lenS,l)
      tryCombs <- matrix(
        c( rep(0,nCombs) , rep(n_parameters((1:l),A)*log(length(X))*xi,nCombs) ),
        byrow = T,nrow = 2
      )
      colnames(tryCombs) <- apply(t(combs), 1, paste0,collapse=",")
      for ( k in 1:nCombs ) {
        G <- t(combs)[k,]
        b <- base_Sja(G,j=NULL,A,base,complete = FALSE)
        tryCombs[1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
      }
      rownames(tryCombs) <- c("log_ML","penalty")
      pML <- apply(tryCombs, 2,sum)
      pML <- sort(pML)[1] #Nested=FALSE & lset=TRUE
      pML <- as.numeric(unlist(strsplit(names(pML), ",")))
    }

  }else{ #if nested=TRUE

    if(!lset){
      ML <- numeric(l)
      penalty <- numeric(l)
      pML <- numeric(l)
      for (i in 1:l) {
        b <- base_Sja(S[1:i],j=NULL,A,base,complete = FALSE)
        ML[i] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
        penalty[i] <- n_parameters(1:i,A)*log(length(X))*xi
        pML[i] <- ML[i]+penalty[i]
        names(pML)[i] <- paste0(S[1:i],collapse = ",")
      }
      nam <- paste0("smallest: ",names(pML)[which( min(pML)==pML )])
      pML <- c(pML,min(pML))
      names(pML)[l+1] <- nam
    }else{
      b <- base_Sja(S[1:l],j=NULL,A,base,complete = FALSE)
      log_ML <- -sum(b$Nxa_Sj*log(b$qax_Sj))
      penalty <- n_parameters((1:l),A)*log(length(X))*xi
      pML <- log_ML+penalty #nested=TRUE & lset=TRUE
      names(pML) <- paste0(S[1:l],collapse = ",")
    }
  }
  if(smallestBIC){
    pML <- pML[l+1]
    pML <- as.numeric(unlist(strsplit(substring(names(pML),11), ",")))
  }
  pML
}

orderedNames <- function(x){
  names(x[order(x)])
}


