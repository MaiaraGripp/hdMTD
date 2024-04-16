#' The Bayesian Information Criterion (BIC) method for inference in MTD models.
#'
#' A function for inference in MTD Markov chains with BIC. It estimates a relevant lag set \eqn{\Lambda} of MTD models using the Bayesian Information Criterion.
#' This means that this method selects the set of lags that minimizes the penalized log likelihood for a given sample.
#'
#'
#' @param X A MTD chain sample.
#' @param d An upper bound for the chains order.
#' @param S A subset of 1:d that contains the relevant lags. If only an upper bound (d) for the chains order is known, the function will set \eqn{S=1:d}.
#' @param minl A lower value for l. If minl==maxl, the function will return the set of size l with elements of S with the smallest BIC.
#' If minl<maxl then the function will consider sets of size minl up to size maxl when looking for the set of elements of S with smallest BIC.
#' @param maxl An upper value for l. If minl==maxl, the function will return the set of size l with elements of S that have the smallest BIC.
#' If minl<maxl then the function will consider all sets of sizes minl up to size maxl when looking for the set of elements of S with smallest BIC.
#' @param xi The BIC constant. Defaulted to 1/2. A smaller xi `(near 0)` reduces the impact of overparameterization.
#' @param A The states space. If not informed, function will set A=unique(X). If there is an element in the states space
#' that doesn't appear in the sample X, A should be informed.
#' @param byl If TRUE, the function will look for the set with smallest BIC by each size from size minl to size maxl, and return the set with smallest BIC for each size.
#' If minl==maxl setting byl TRUE or FALSE makes no difference since the function will only calculate the BIC for models with maxl elements in the relevant lag set.
#' @param BICvalue If TRUE, the function will also return the calculated values of the BIC for the estimated relevant lags sets.
#' @param warning If TRUE, may return warnings.
#' @param ... Used to accommodate any extra arguments passed by the [hdMTD()] function.
#'
#' @details Note that the order of the chain matters at the counting of sequences step. If we run the function with a certain order parameter d, all the sequences of size d will be counted, and when
#'calculating the BIC for a subset S (of 1:d), the function will sum over those counts. This means that if the order changes from d to d', some counts will change, especially when d is
#'distant from d'. As a result, the function might calculate a different BIC for the same set S.
#'
#' @return Returns estimations, using Bayesian Information Criterion (BIC), of the relevant lag set of size \eqn{minl,minl+1, \dots maxl}.
#' For example, suppose \eqn{S=c(9,5,1,2)}, \eqn{minl=1}, and \eqn{maxl=3}. This function
#' will iterate from \eqn{l=1 to l=3} and calculate the BIC assuming the MTD model has \eqn{\Lambda=c(1)}, then it will calculate the BIC for \eqn{\Lambda=c(2)}, followed by \eqn{\Lambda=c(5)}, and finally \eqn{\Lambda=c(9)}.
#' Afterwards, it will assume \eqn{\Lambda=c(1,2)}, \eqn{\Lambda=c(1,5)}, ..., \eqn{\Lambda=c(5,9)}, and subsequently \eqn{\Lambda=c(1,2,5)} until \eqn{\Lambda=c(2,5,9)}. After calculating the BIC for each model, it will return the
#' sets with the smallest BIC among all sizes. However, if the option byl=TRUE is selected, it will return the sets with the smallest BIC for each size. Parameters maxl and minl may have the same value, in which case the set with the smallest
#' BIC will be chosen among the sets of size maxl.
#' @export
#'
#' @examples
#' X <- perfectSample(MTDmodel(Lambda=c(2,5),A=c(0,2),lam0=0.05),500)
#'hdMTD_BIC (X,d=6,minl=1,maxl=1)
#'hdMTD_BIC (X,d=5,minl=2,maxl=2,BICvalue = TRUE)
#'
hdMTD_BIC <- function(X,d,S=1:d,minl=1,maxl=max(S),
                          xi=1/2,A=NULL,byl=FALSE,BICvalue=FALSE,warning=FALSE,...){
  #Checking inputs
  X <- checkSample(X)
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A is not informed. Code will set A <- sort(unique(X)).")
    }
    A <- sort(unique(X))
  }
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it must include all states that occur in the sample.")
  }
  if( !is.numeric(d) ||
      d<2 ||
      (d %% 1)!=0  ||
      d<max(S)){stop("The order d must be an integer number greater than 2 or the greatest element in S.")
    }
  if(length(S) < 2  ||
     length(dim(S))!=0 ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")
    }
  if ( is.na(minl) ||
       !is.numeric(minl) ||
       minl%%1 != 0 ||
       minl>length(S) ||
       minl <=0 ) {
      stop("The minl value is not valid. minl should be a positive integer lower or equal to the number of elements in S.")
  }
  if ( is.na(maxl) ||
       !is.numeric(maxl) ||
       maxl%%1 != 0 ||
       maxl>length(S) ||
       maxl<minl) {
    stop("The maxl value is not valid. maxl should be a positive integer lower or equal to the number of elements in S, and greater or equal to minl.")
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("The BIC constant xi value is not valid. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  if(!is.logical(byl)){stop("byl must me a logical argument")}
  if(!is.logical(BICvalue)){stop("BICvalue must me a logical argument")}
  if(!is.logical(warning)){stop("warning must me a logical argument")}

  A <- sort(A)
  lenS <- length(S)
  S <- sort(S)
  base <- countsTab(X,d)

## Calculating penalized loglikelihood
  if(maxl==minl){
        nCombs <- choose(lenS,minl) #number of possible sets with size minl of elements of S
        tryCombs <- matrix( c( rep(0,nCombs), rep(n_parameters(1:minl,A)*log(length(X))*xi,nCombs) ),
                            byrow = T,nrow = 2 )
        #since minl=maxl the number of parameters is the same, and the penalty is the same for any set
        aux <- apply(t(combn(S,minl)), 1, paste0,collapse=",")
        #aux is a vector with all possible sets of size minl with elements of S
        for ( k in 1:nCombs ) {
          G <- as.numeric(unlist(strsplit(aux[k], ",")))
          b <- freqTab(S=G,j=NULL,A=A,countsTab=base,complete = FALSE)
          tryCombs[1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
        }
        colnames(tryCombs) <- aux
        rownames(tryCombs) <- c("log_ML","penalty")
        pML <- apply(tryCombs,2,sum) # -loglikelihood + penalty
        pML <- sort(pML)[1] # max
        if(!BICvalue){
          pML <- as.numeric(unlist(strsplit(names(pML), ",")))
        }
      }else{ #maxl>minl
        tryCombs <- list()
        cont <- 1
        # does the same thing as when minl=maxl but for all minl<=l<=maxl.
        # So tryCombs is now a list, and the penalty changes with l
        for (i in minl:maxl) {
          nCombs <- choose(lenS,i)
          tryCombs[[cont]] <- matrix( rep(0,2*nCombs) ,byrow = T,nrow = 2 )
          tryCombs[[cont]][2,] <- n_parameters(Lambda=(1:i),A)*log(length(X))*xi
          aux <- apply(t(combn(S,i)), 1, paste0,collapse=",")
          for ( k in 1:nCombs) {
            G <- as.numeric(unlist(strsplit(aux[k], ",")))
            b <- freqTab(S=G,j=NULL,A=A,countsTab=base,complete = FALSE)
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


