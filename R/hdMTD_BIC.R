
#' The Bayesian Information Criterion (BIC) method for inference in MTD models
#'
#' A function for estimating the relevant lag set \eqn{\Lambda} of a Markov chain using
#' Bayesian Information Criterion (BIC). This means that this method selects the set of lags
#' that minimizes a penalized log likelihood for a given sample, see *References* below for
#' details on the method.
#'
#' @param X A vector or single-column data frame containing a chain sample (`X[1]` is the most recent).
#' @param d A positive integer representing an upper bound for the chain order.
#' @param S A numeric vector of positive integers from which this function will select
#' a set of relevant lags. Typically, \code{S} is a subset of \code{1:d}. If \code{S}
#' is not provided, by default \code{S=1:d}.
#' @param minl A positive integer. \code{minl} represents the smallest length of any relevant lag
#'  set this function might return. If \code{minl == maxl}, this function will return the
#'  subset of \code{S} of length \code{minl} with the lowest BIC. If \code{minl < maxl}, the function
#'  will consider subsets ranging from length \code{minl} to length \code{maxl} when searching for
#'  the subset of \code{S} with the smallest BIC.
#' @param maxl A positive integer equal to or greater than \code{minl} but less than the number
#'  of elements in \code{S} (\code{maxl = length(S)} is accepted but in this case the output will
#'  always be \code{S}). \code{maxl} represents the largest length of any relevant lag set this
#'  function might return.
#' @param xi The BIC penalization term constant. Defaulted to 1/2. A smaller \code{xi} `(near 0)`
#' reduces the impact of overparameterization.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A=sort(unique(X))}.
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
#' in the MTD, which affects the number of parameters in the penalization term. Defaulted
#' to \code{maxl}. See more in *Details*.
#' @param warning Logical. If \code{TRUE}, the function warns the user when \code{A} is set automatically.
#' @param ... Additional arguments (not used in this function, but maintained for compatibility with [hdMTD()].
#'
#' @details Note that the upper bound for the order of the chain (\code{d}) affects the estimation
#' of the transition probabilities. If we run the function with a certain order parameter \code{d},
#' only the sequences of length \code{d} that appeared in the sample will be counted. Therefore,
#' all transition probabilities, and hence all BIC values, will be calculated with respect to
#'that \code{d}. If we use another value for \code{d} to run the function, even if the output
#' agrees with that of the previous run, its BIC value might change a little.
#'
#' The parameter \code{zeta} indicates the the number of distinct matrices pj in the MTD.
#' If \code{zeta = 1}, all matrices \eqn{p_j} are identical; if \code{zeta = 2} there exists
#' two groups of distinct matrices and so on. The largest value for \code{zeta} is \code{maxl}
#' since this is the largest number of matrices \eqn{p_j}. When \code{minl<maxl},
#' for each \code{minl} \eqn{\leq} \code{l} \eqn{\leq} \code{maxl}, \code{zeta = min(zeta,l)}.
#' If \code{single_matrix = TRUE} then \code{zeta} is set to 1.
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
#' X <- testChains[, 1]
#' hdMTD_BIC (X, d = 6, minl = 1, maxl = 1)
#' hdMTD_BIC (X,d = 3,minl = 1, maxl = 2, BICvalue = TRUE)
#'
hdMTD_BIC <- function(X, d, S = seq_len(d), minl = 1, maxl = length(S),
                      xi = 1/2, A = NULL, byl = FALSE, BICvalue = FALSE,
                      single_matrix = FALSE, indep_part = TRUE,
                      zeta = maxl, warning = FALSE,...){
  # Validate inputs
  X <- checkSample(X)
  check_hdMTD_BIC_inputs(X, d, S, minl, maxl, xi, A, byl, BICvalue, single_matrix,
                         indep_part, zeta, warning)

  # Set the state space if not provided
  if(length(A) == 0) { A <- sort(unique(X)) } else { A <- sort(A) }

  A <- sort(A)
  lenS <- length(S)
  S <- sort(S)
  base <- countsTab(X, d)

  # Compute penalized log-likelihood
  if(maxl == minl){

    nCombs <- choose(lenS, minl) # Number of size minl subsets of S
    Combs <- t(combn(S, minl)) # All size minl subsets of S

    tryCombs <- matrix(rep(0, 2*nCombs), byrow = T, nrow = 2)
    # tryCombs[1, ] will store log likelihoods for each possible size minl subset
    # of S as set of relevant lags. tryCombs[2, ] will store their penalization term
    colnames(tryCombs) <- apply(Combs, 1, paste0, collapse = ",")
    rownames(tryCombs) <- c("log_ML","penalty")

    # Compute the number of parameters in an MTD with minl relevant lags
      n_param <- n_parameters(Lambda = seq_len(minl), A = A,
                              single_matrix = single_matrix,
                              indep_part = indep_part, zeta = zeta)
      # Compute penalizations
      tryCombs[2, ] <- n_param * log(length(X)) * xi

      # Compute log likelihoods for each size minl subset of S as set of relevant lags.
      for (k in seq_len(nCombs)) {
          b <- freqTab(S = Combs[k, ], j = NULL, A = A, countsTab = base, complete = FALSE)
          tryCombs[1, k] <- -sum(b$Nxa_Sj * log(b$qax_Sj))
      }

      pML <- colSums(tryCombs) # -loglikelihood + penalty
      pML <- sort(pML)[1] # min (did not use min(pML) since names(min(pML))=NULL)

      if(!BICvalue){ # only names(pML) is stored
        pML <- as.numeric(unlist(strsplit(names(pML), ",")))
      }

     # If maxl > minl the function repeats the algorithm above for all
     # minl<=i<=maxl as the sizes of subsets of S
  }else{

      tryCombs <- list() # Initiate a list of matrices tryCombs

      for (i in minl:maxl) {
          cont <- i - minl + 1

          nCombs <- choose(lenS, i) # Number of size i subsets of S
          Combs <- t(combn(S, i)) # All size i subsets of S

          tryCombs[[cont]] <- matrix(rep(0, 2*nCombs), byrow = T, nrow = 2)
          colnames(tryCombs[[cont]]) <- apply(Combs, 1, paste0, collapse = ",")
          rownames(tryCombs[[cont]]) <- c("log_ML", "penalty")

          # Compute the number of parameters in an MTD with i relevant lags
          n_param <- n_parameters(Lambda = seq_len(i), A = A, single_matrix = single_matrix,
                                  indep_part = indep_part, zeta = min(zeta, i))
          # Compute penalizations
          tryCombs[[cont]][2, ] <- n_param * log(length(X)) * xi

          # Compute log likelihoods for each size i subset of S as set of relevant lags.
          for (k in seq_len(nCombs)) {
            b <- freqTab(S = Combs[k, ], j = NULL, A = A, countsTab = base, complete = FALSE)
            tryCombs[[cont]][1, k] <- -sum(b$Nxa_Sj * log(b$qax_Sj))
          }
      }

      pML <- sapply(tryCombs, colSums, simplify = FALSE) # -loglikelihood + penalty
          # simplify needs to be FALSE because if lenS is odd, minl=(lenS-1)/2
          # and maxl=minl+1 all matrices in tryCombs will have the same ncol.

      pML <- sapply(pML, function(x) sort(x)[1]) # Lowest BIC value per set size

      if(byl){
        smallest <- pML[order(pML)][1] # Set with the lowest BIC
        pML <- c(pML, smallest)
        names(pML)[length(pML)] <- paste0("smallest: ",names(smallest))
        if(!BICvalue){
          pML <- names(pML)
        }
      }else{
        pML <- sort(pML)[1]
        if(!BICvalue){
          pML <- as.numeric(unlist(strsplit(names(pML), ",")))
        }
      }
  }
  pML
}



