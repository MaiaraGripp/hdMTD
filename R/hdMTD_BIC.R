#' The Bayesian Information Criterion (BIC) method for inference in MTD models
#'
#' A function for estimating the relevant lag set \eqn{\Lambda} of a Markov chain using
#' Bayesian Information Criterion (BIC). This means that this method selects the set of lags
#' that minimizes a penalized log likelihood for a given sample, see *References* below for
#' details on the method.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param d A positive integer representing an upper bound for the chain order.
#' @param S A numeric vector of positive integers from which this function will select
#' a set of relevant lags. Typically, \code{S} is a subset of \code{1:d}. If \code{S}
#' is not provided, by default \code{S=1:d}.
#' @param minl A positive integer. \code{minl} represents the smallest length of any relevant lag
#'  set this function might return. If \code{minl == maxl}, this function will return exactly
#'  \code{minl} lags in its output, specifically the subset of \code{S} with \code{minl} elements
#'  that has the smallest BIC. If \code{minl < maxl}, the function will consider subsets ranging
#'  from length \code{minl} to length \code{maxl} when searching for the subset of \code{S} with
#'  the smallest BIC.
#' @param maxl A positive integer equal to or greater than \code{minl} but less than the number
#'  of elements in \code{S} (\code{maxl = length(S)} is accepted but in this case the output will
#'  always be \code{S}). \code{maxl} represents the largest length of any relevant lag set this
#'  function might return.
#' @param xi The BIC penalization term constant. Defaulted to 1/2. A smaller \code{xi} `(near 0)`
#' reduces the impact of overparameterization.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A=unique(X)}.
#' @param byl Logical. If \code{TRUE}, the function will look for the set with smallest BIC by each
#' length (from  \code{minl} to \code{maxl}), and return the set with smallest BIC for each length.
#' If \code{minl==maxl} setting \code{byl=TRUE} or \code{FALSE} makes no difference, since the
#' function will only calculate the BIC for sets with \code{maxl} elements in the relevant lag set.
#' @param BICvalue Logical. If \code{TRUE}, the function will also return the calculated values of
#'  the BIC for the estimated relevant lag sets.
#' @param single_matrix Logical. If \code{TRUE}, the chain sample is thought to come from an MTD model
#' where the stochastic matrices \eqn{p_j} are constant across all lags \eqn{j\in \Lambda}. In practice,
#' this means the user believes the stochastic matrices for every lag in \code{S} are the same, which reduces
#' the number of parameters in the penalization term.
#' @param indep_part Logical. If \code{FALSE} there is no independent distribution and \eqn{\lambda_0=0} which
#' reduces the number of parameters in the penalization term.
#' @param zeta A positive integer representing the number of distinct matrices \eqn{p_j}
#' in the MTD. By default, it is equal to \code{maxl}, which is the maximum
#' allowable value and indicates that all matrices \eqn{p_j} are distinct. If
#' \code{zeta=1}, all matrices \eqn{p_j} are identical; if \code{zeta=2}
#' there exists two groups of distinct matrices and so on. When \code{minl<maxl},
#' for each \eqn{minl\leq l\leq maxl}, \code{zeta=min(zeta,l)}. If \code{single_matrix=TRUE}
#' then \code{zeta} is set to 1.
#' @param warning Logical. If \code{TRUE}, informs the user if the state space was set as \code{A=unique(X)}.
#' @param ... Other parameters. This is used to accommodate any additional argument passed
#' to [hdMTD_BIC()] through the [hdMTD()] function.
#'
#' @details Note that the upper bound for the order of the chain (\code{d}) affects the estimation
#' of the transition probabilities. If we run the function with a certain order parameter \code{d},
#' only the sequences of length \code{d} that appeared in the sample will be counted. Therefore,
#' all transition probabilities, and hence all BIC values, will be calculated with respect to
#'that \code{d}. If we use another value for \code{d} to run the function, even if the output
#' agrees with that of the previous run, its BIC value might change a little.
#'
#' @references Imre Csiszár, Paul C. Shields. "The consistency of the BIC Markov order estimator."
#'  The Annals of Statistics, 28(6) 1601-1619 December2000.
#' [DOI 10.1214/aos/1015957472](https://doi.org/10.1214/aos/1015957472).
#'
#' @return Returns a vector with the estimated relevant lag set using BIC. It might return more
#' than one set if \code{minl < maxl} and \code{byl = TRUE}. Additionally, it can return the value
#' of the penalized likelihood for the outputted lag sets if \code{BICvalue = TRUE}.
#' @export
#'
#' @examples
#' X <- testChains[,1]
#'hdMTD_BIC (X,d=6,minl=1,maxl=1)
#'hdMTD_BIC (X,d=3,minl=1,maxl=2,BICvalue = TRUE)
#'
hdMTD_BIC <- function(X,d,S=1:d,minl=1,maxl=length(S),
                      xi=1/2,A=NULL,byl=FALSE,BICvalue=FALSE,
                      single_matrix=FALSE,indep_part=TRUE,
                      zeta=maxl,warning=FALSE,...){
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
  if(!is.logical(indep_part)){stop("indep_part must me a logical argument")}
  if(!is.logical(single_matrix)){stop("single_matrix must me a logical argument")}
  if(single_matrix==FALSE){
    if ( is.na(zeta) ||
         !is.numeric(zeta) ||
         zeta%%1 != 0 ||
         zeta>maxl ||
         zeta<1 ) {
      stop("The zeta value is not valid. zeta should be a positive integer
           representing the number of distinct matrices pj in the MTD.")
    }
  }
  A <- sort(A)
  lenS <- length(S)
  S <- sort(S)
  base <- countsTab(X,d)

## Calculating penalized log-likelihood
  if(maxl==minl){
        nCombs <- choose(lenS,minl) #number of possible sets with length minl of elements of S
        tryCombs <- matrix( c(rep(0,nCombs),
                              rep(n_parameters(Lambda=1:minl,
                                               A=A,
                                               single_matrix=single_matrix,
                                               indep_part=indep_part,
                                               zeta=zeta)*log(length(X))*xi,nCombs)
                             ),byrow = T,nrow = 2 )
        aux <- apply(t(combn(S,minl)), 1, paste0,collapse=",")
        #aux is a vector with all possible sets of length minl with elements of S
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
          tryCombs[[cont]][2,] <- n_parameters(Lambda=(1:i),
                                               A=A,
                                               single_matrix = single_matrix,
                                               indep_part = indep_part,
                                               zeta=min(zeta,i))*log(length(X))*xi

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
            #will format pML differently and remove the colnames, which we can't have...
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


