#' A function for inference in MTD Markov chains with BIC.
#'
#' Estimates a relevant lag set \eqn{\Lambda} of MTD models using BIC.
#'
#' @param X A Markov chain.
#' @param d An upper threshold for the chains order.
#' @param S A set of all possible relevant lags, if only an upper threshold (d) for the chains order is known, the function will set \eqn{S=1:d}.
#' @param minl Lower value for l. If minl==maxl the function will return the set of size l with elements of S with smallest BIC.
#' If minl<maxl then the function will consider sets of size minl up to size maxl when looking for the set of elements of S with smallest BIC.
#' @param maxl Upper value for l. If minl==maxl the function will return the set of size l with elements of S with smallest BIC.
#' If minl<maxl then the function will consider all sets of sizes minl up to size maxl when looking for the set of elements of S with smallest BIC.
#' @param xi The BIC constant. Defaulted to 1/2. Smaller xi `(near 0)` reduces the impact of overparameterization.
#' @param A The states space. If not informed function will set A=unique(X). You should inform A specially if there is an element in the states space
#' that doesn't appear in the sample X.
#' @param byl If TRUE will look for the set with smallest BIC by each size from size minl to size maxl, and return the set with smallest BIC for each size.
#' If minl==maxl setting byl TRUE or FALSE makes no difference since the function will only calculate the BIC for models with maxl elements in the relevant lag set.
#' @param BICvalue If TRUE the function will also return the calculated values of the BIC for the estimated relevant lag sets.
#' @param warning If TRUE may return warnings.
#'
#' @details Note that the order of the chain matters at the counting of sequences step. If we run the function with some order parameter d, all the sequences of size d will be counted and when
#' calculating the BIC for some subset S of the set 1:d the function will sum over those counts. That means that if the order d changes for say d', some counts will change, especially when d is
#' distant from d'. That might cause the function to calculate a different BIC for the same set S.
#'
#' @return Returns estimations, using Bayesian Information Criterion (BIC), of the relevant lag set of size \eqn{minl,minl+1, \dots maxl}.
#'  For example, suppose \eqn{S=c(9,5,1,2)}, \eqn{minl=1} and \eqn{maxl=3}. This function
#' will run from \eqn{l=1 to l=3} and calculate the BIC supposing the MTD model has \eqn{\Lambda=c(1)}, then will calculate the BIC for \eqn{\Lambda=c(2)}, then \eqn{\Lambda=c(5)}, and finally \eqn{\Lambda=c(9)}.
#' Then it will suppose \eqn{\Lambda=c(1,2)}, \eqn{\Lambda=c(1,5)}, ... , \eqn{\Lambda=c(5,9)}, and then \eqn{\Lambda=c(1,2,5)} til \eqn{\Lambda=c(2,5,9)}. After calculating all BIC for each model will return the
#' sets with smallest BIC between all sizes. Unless option byl=TRUE in which case will return the sets with smallest BIC for each size. Parameters maxl and minl may have the same value, in which case the smallest
#' BIC set will will be chosen among the sets of size maxl.
#' @export
#'
#' @examples
#' X <- perfectSample(MTDmodel(Lambda=c(2,4,6),A=c(0,2,5),w0=0.05),2000)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=1)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=1,BICvalue = TRUE)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=4)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=4, BICvalue = TRUE)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=4,xi=1.2)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=4,byl = TRUE)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=4,byl = TRUE, BICvalue = TRUE)
#'sparseMarkov_BIC(X,d=9,minl=1,maxl=4,xi=0.8,byl = TRUE)
#'sparseMarkov_BIC(X,d=9,S=c(2,8,7,6,4),minl=3,maxl=3)
#'sparseMarkov_BIC(X,d=9,S=c(2,4,7,6,8),minl=1,maxl=3)
#'sparseMarkov_BIC(X,d=9,S=c(2,4,7,6,8),minl=3,maxl=3,byl=TRUE)
#'sparseMarkov_BIC(X,d=9,S=c(2,8,7,6,4),minl=1,maxl=3,xi=0.9,byl=TRUE)
#'sparseMarkov_BIC(X,d=9,S=c(2,8,7,6,4),minl=1,maxl=3,xi=0.9,byl=TRUE, BICvalue = TRUE)
#'
sparseMarkov_BIC <- function(X,d,S=1:d,minl=1,maxl=max(S),
                              xi=1/2,A=NULL,byl=FALSE,BICvalue=FALSE,warning=FALSE){
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
  if( !is.numeric(d) ||
      d<2 ||
      (d %% 1)!=0  ||
      d<max(S)){stop("The order d must be an integer number greater than 2 or the greatest element in S.")
    }
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")
    }
  if ( is.na(minl) ||
       !is.numeric(minl) ||
       minl%%1 != 0 ||
       minl>length(S) ||
       minl <=0 ) {
    stop("minl value is not valid. minl should be a positive integer lower or equal to the number of elements in S.")
  }
  if ( is.na(maxl) ||
       !is.numeric(maxl) ||
       maxl%%1 != 0 ||
       maxl>length(S) ||
       maxl<minl) {
    stop("maxl value is not valid. maxl should be a positive integer lower or equal to the number of elements in S, and greater or equal to minl.")
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("BIC constant xi value is not valid. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  if(!is.logical(byl)){stop("byl must me a logical argument")}
  if(!is.logical(BICvalue)){stop("BICvalue must me a logical argument")}
  if(!is.logical(warning)){stop("warning must me a logical argument")}

  #gatherin inputs
  A <- sort(A)
  lenS <- length(S)
  S <- sort(S)
  base <- shapeSample(X,d)
  #\.

  if(maxl==minl){
        nCombs <- choose(lenS,minl)
        tryCombs <- matrix( c( rep(0,nCombs), rep(n_parameters(1:minl,A)*log(length(X))*xi,nCombs) ),
                            byrow = T,nrow = 2 )
        aux <- apply(t(combn(S,minl)), 1, paste0,collapse=",")
        for ( k in 1:nCombs ) {
          G <- as.numeric(unlist(strsplit(aux[k], ",")))
          b <- base_Sja(G,j=NULL,A,base,complete = FALSE)
          tryCombs[1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
        }
        colnames(tryCombs) <- aux
        rownames(tryCombs) <- c("log_ML","penalty")
        pML <- apply(tryCombs,2,sum)
        pML <- sort(pML)[1]
        if(!BICvalue){
          pML <- as.numeric(unlist(strsplit(names(pML), ",")))
        }
      }else{ #maxl>minl
        tryCombs <- list()
        cont <- 1
        for (i in minl:maxl) {
          nCombs <- choose(lenS,i)
          tryCombs[[cont]] <- matrix( rep(0,2*nCombs) ,byrow = T,nrow = 2 )
          tryCombs[[cont]][2,] <- n_parameters(Lambda=(1:i),A)*log(length(X))*xi
          aux <- apply(t(combn(S,i)), 1, paste0,collapse=",")
          for ( k in 1:nCombs) {
            G <- as.numeric(unlist(strsplit(aux[k], ",")))
            b <- base_Sja(G,j=NULL,A,base,complete = FALSE)
            tryCombs[[cont]][1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
          }
          colnames(tryCombs[[cont]]) <- aux
          rownames(tryCombs[[cont]]) <- c("log_ML","penalty")
          cont <- cont+1
        }
          if(maxl-minl==1){ #this solves a problem that may occur if all vectors in
            #tryCombs list have the same length. If that is the case, the function sapply
            #will format pML differently and remove the colnames, which we cant have...
            if( ncol(tryCombs[[1]])==ncol(tryCombs[[2]]) ){
             tryCombs[[2]] <- cbind(tryCombs[[2]],tryCombs[[2]][,1])
             colnames(tryCombs[[2]])[ncol(tryCombs[[2]])] <- colnames(tryCombs[[2]])[1]
            }
          }
        pML <- sapply(tryCombs, apply, 2,sum)

          if(byl){
            smallest <- names(unlist(pML))[order(unlist(pML))][1]
            nam <- sapply(sapply(pML, orderedNames), dplyr::first )
            pML <- sapply(sapply(pML,sort),dplyr::first)
            pML <- c( pML,unlist(pML)[order(unlist(pML))][1] )
            names(pML) <- c(nam,paste0("smallest: ",smallest))
              if(!BICvalue){
                pML <- names(pML)
              }
          }else{
            pML <- unlist(pML)[order(unlist(pML))][1]
              if(!BICvalue){
                pML <- as.numeric(unlist(strsplit(names(pML), ",")))
              }
          }
    }
   pML
}

orderedNames <- function(x){
  names(x[order(x)])
}

