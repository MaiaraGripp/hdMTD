#' EM estimation of MTD parameters
#'
#' Estimation of MTD parameters through the Expectation Maximization (EM) algorithm.
#'
#' @details Regarding the \code{M} parameter: it functions as a stopping
#' criterion within the EM algorithm. When the difference between
#' the log-likelihood computed with the newly estimated parameters
#' and that computed with the previous parameters falls below \code{M},
#' the algorithm halts. Nevertheless, if the value of \code{nIter}
#' (which represents the maximum number of iterations) is smaller
#' than the number of iterations required to meet the \code{M} criterion,
#' the algorithm will conclude its execution when \code{nIter} is reached.
#' To ensure that the \code{M} criterion is effectively utilized, we
#' recommend using a higher value for \code{nIter}, which is set to a
#' default of 100.
#'
#' Concerning the \code{init} parameter, it is expected to be a list
#' comprising either 2 or 3 entries. These entries consist of:
#' an optional vector named \code{p0}, representing an independent
#' distribution (the probability in the first entry of \code{p0} must be
#' that of the smallest element in \code{A} and so on), a required list
#' of matrices \code{pj}, containing a stochastic matrix for each
#' element of \code{S} ( the first matrix must refer to the smallest
#' element of \code{S} and so on), and a vector named \code{lambdas} representing
#' the weights, the first entry must be the weight for \code{p0}, and then one entry
#' for each element in \code{pj} list. If your MTD model does not have an independent
#' distribution \code{p0}, set \code{init$lambda[1]=0}.
#'
#'
#' @references
#' This function was created based on the following article:
#' Lebre, Sophie & Bourguignon, Pierre-Yves. (2008). An EM algorithm for estimation
#' in the Mixture Transition Distribution model. Journal of Statistical Computation
#' and Simulation. 78. \url{https://doi.org/10.1080/00949650701266666}.
#'
#'
#' @param X A vector or single-column data frame containing an MTD chain sample (`X[1]` is the most recent).
#' @param S A numeric vector of positive integers. Typically, \code{S} represents a set of relevant lags.
#' @param M A stopping point for the EM algorithm. If \code{M=NULL} the algorithm will run
#' for a total of \code{nIter} iteractions.
#' @param init A list with initial parameters: \code{p0} (optional), \code{lambdas} (required),
#'  \code{pj} (required). The entries in \code{lambdas} are weights for the distribution \code{p0}
#'  and the distributions present in the list \code{pj}. Therefore, the order in which the elements
#'  appear in the vector \code{lambdas} is important for correct assignment. Please refer to the
#'  *Details* section for more information.
#' @param iter Logical. If \code{TRUE}, returns the number of iterations of the
#' algorithm, that is, the number of times the initial parameters were updated.
#' @param nIter An integer positive number with the maximum number of iterations.
#' @param oscillations Logical. If \code{TRUE}, the function will return the estimated oscillations
#' for the updated model along with the estimated parameters.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A=unique(X)}.
#'
#' @export
#'
#' @return A list with the estimated parameters of the MTD model.
#'
#' @examples
#' # Simulating data.
#' # Model:
#' MTD <- MTDmodel(Lambda=c(1,10),A=c(0,1),lam0=0.01)
#' # Sampling a chain:
#' X <- hdMTD::perfectSample(MTD,N=2000)
#'
#' # Initial Parameters:
#' init <- list('p0'=c(0.4,0.6),'lambdas'=c(0.05,0.45,0.5),
#'   'pj'=list(matrix(c(0.2,0.8,0.45,0.55),byrow = TRUE,ncol=2),
#'    matrix(c(0.25,0.75,0.3,0.7),byrow = TRUE,ncol=2)))
#'
#' # MTDest() ------------------------------------
#' MTDest(X,S=c(1,10),M=1,init)
#' MTDest(X,S=c(1,10),init=init,iter = TRUE)
#' MTDest(X,S=c(1,10),init=init,iter = TRUE,nIter=5)
#' MTDest(X,S=c(1,10),init=init,oscillations = TRUE)
#'
MTDest <- function(X, S, M = 0.01, init, iter = FALSE, nIter = 100, A = NULL, oscillations = FALSE) {

  # Validate inputs
  X <- checkSample(X)
  check_MTDest_inputs(X, S, M, init, iter, nIter, A, oscillations)

  if (is.null(A)) { A <- sort(unique(X)) } else { A <- sort(A) }

  indep <- TRUE
  if (length(init$p0) == 0 || all(init$p0 == 0)) {
      if (init$lambdas[1] == 0 ) {
          indep <- FALSE
          init$p0 <- rep(0, lenA)
      }
  }

  # Creates an MTD object and validates inputs
  initMTD <- MTDmodel(Lambda = S, A = A, lam0 = init$lambdas[1],
                      lamj = init$lambdas[-1], pj = init$pj, p0 = init$p0)
  checkMTD(initMTD)

  rS <- S
  S <- sort(S, decreasing = TRUE)
  lenS <- length(S)
  lenX <- length(X)
  lenA <- length(A)

  S0 <- c(S, 0)
  lenS0 <- lenS + 1
  base <- countsTab(X, S[1])
  baseSja <- freqTab(S, j = NULL, A, base)

  # Check if initial probabilities are compatible with the sample
  Pinit <- t(initMTD$P)
  dim(Pinit) <- NULL
  if(any(round(Pinit, 6) == 0)) {
      if (baseSja$qax_Sj[which(round(Pinit, 5) == 0)] > 0) {
          stop("The initial parameters in init aren't compatible with the sample X.
           There are sequences that appear in the sample for which
           the provided initial probability is zero.
           Run: 'MTDmodel( Lambda = S, A = sort(unique(X)), lam0 = init$lambdas[1],
           lamj = init$lambdas[-1], pj = pj, p0 = init$p0 )$P' and verify if there
           are null entries.")
      }
  }

  # The EM algorithm

  contInt <- 0
  distlogL <- NULL

  # Parameter update loop. It terminates if the distance between the
  # log likelihood of the initial and updated models is less than M, or if the
  # maximum number of iterations (nIter) is reached.
  repeat{

        initMTD <- MTDmodel(rS, A, lam0 = init$lambdas[1], lamj = init$lambdas[-1],
                            pj = init$pj, p0 = init$p0)

        # Compute log likelihood
        if(any(baseSja$Nxa_Sj == 0)){
            # If TRUE is possible that some probabilities in initMTD$P are also 0, and
            # log(P(a|x))*N(xa) = log(0)*0 = -inf*0 = Nan. To prevent this use the auxiliary
            # function prodinf (see utils.R). With it, terms resulting from -inf*0 are set to 0.
            initLogL <- sum(prodinf(log(as.vector(t(initMTD$P))), baseSja$Nxa_Sj))
        }else{
            initLogL <- sum(log(as.vector(t(initMTD$P))) * baseSja$Nxa_Sj)
        }

        # Step E (Expectation)

        # Creates a stochastic matrix with the lambdas
        pSja <- matrix(rep(rev(init$lambdas), lenA^lenS0), byrow = T, ncol = lenS0)
        colnames(pSja) <- S0

        # pSja column referent to lam0 is multiplied by p0
        pSja[, lenS0] <- pSja[, lenS0] * init$p0
        # pSja column referent to lamj is multiplied by pj, j\in S
        indexA <- expand.grid(rep(list(seq_len(lenA)), lenS0))[, order(lenS0:1)]
        cont <- lenS
        for (i in seq_len(lenS)) { #col
            pj <- init$pj[[cont]]
            for (j in seq_len(lenA^lenS0)) { #row
                pSja[j,i] <- pSja[j,i]*pj[indexA[j,i],indexA[j,(lenS0)]]
            }
            cont <- cont-1
        }
        colnames(pSja) <- paste0("lamXp_", S0) # lambdas \times pj()
        # Normalize pSja so each row is a distribution or a vector of 0
        norm <- rowSums(pSja)
        norm[which(norm == 0)] <- 1 # to avoid 0/0
        Pj_xa_S <- pSja/norm
        # Pj_xa_S is a matrix with the probabilities of a hidden state variable Z=j,
        # j \in S\cup{0}, given any possible size lenS0 sequence with elements of A
        colnames(Pj_xa_S) <- paste0("P", S0, "xa_S")
        rownames(Pj_xa_S) <- apply(expand.grid(rep(list(A), lenS0))[, order(lenS0:1)],
                                   1, paste0,collapse = "")


        # End of step E

        # Step M (Maximization)

        # lambdas update
        NxaXPjxa <- Pj_xa_S * baseSja$Nxa_Sj
        colnames(NxaXPjxa) <- paste0("NXP", S0, "xa")
        end_lambdas <- apply(NxaXPjxa, 2, sum)/(lenX-S[1])
        names(end_lambdas) <- paste0("lam-", S0)
        end_lambdas <- rev(end_lambdas)
        # end_lambdas is the updated lambdas vector

        # p0 update
        NPSja <- cbind(baseSja[, seq_len(lenS0)], NxaXPjxa)
        end_p0 <- rep(sum(NPSja$NXP0xa), lenA)
        if (indep) {
            for (i in seq_len(lenA)) {
                end_p0[i] <- sum((NPSja %>%
                  dplyr::filter(a == A[i]))$NXP0xa)/end_p0[i]
            }
        }
        names(end_p0) <- paste0("p_0(", A, ")")
        # end_p0 is the updated independent distribution p0

        # pj update
        end_pj <- list()
        for (j in seq_len(lenS)) {
            aux_pj <- matrix(0, ncol = lenA, nrow = lenA)

          for (i in seq_len(lenA)) {
              aux_NPSja <- NPSja %>%
                dplyr::filter(NPSja[ ,j] == A[i]) # fix past x_j=A[i]

              aux_pj[i, ] <- sum(aux_NPSja %>%
                dplyr::select_at(paste0("NXP", S[j], "xa")))

              for (k in seq_len(lenA)) {
                aux_pj[i, k] <- sum(aux_NPSja %>%
                  dplyr::filter(a == A[k]) %>%
                  dplyr::select_at(paste0("NXP", S[j], "xa")))/aux_pj[i, k]
              }
          }
          colnames(aux_pj) = rownames(aux_pj) <- A
          end_pj[[lenS - j + 1]] <- aux_pj
        }
        names(end_pj) <- paste0("p_-",rev(S))
        # end_pj is the list with updated pj matrices

        # Compute updated MTD model
        endMTD <- MTDmodel(rS, A, lam0 = end_lambdas[1], lamj = end_lambdas[-1],
                           pj = end_pj, p0 = end_p0)
        endlist <- list("lambdas" = end_lambdas, "pj" = end_pj, "p0" = end_p0)

        # Compute updated log likelihood
        if (any(baseSja$Nxa_Sj == 0)) {
            endLogL <- sum( prodinf( log(as.vector(t(endMTD$P))), baseSja$Nxa_Sj ) )
        }else{
            endLogL <- sum( log(as.vector(t(endMTD$P)))*baseSja$Nxa_Sj )
        }

        # Compute difference between updated and former log likelihood
        distlogL[contInt + 1] <- abs(endLogL - initLogL)

        if (length(M) == 1){ # If the log likelihood difference is too small stop EM algorithm
            if(distlogL[contInt + 1] < M){ break }
        }

        init <- endlist # initial parameters are replaced with updated parameters
        contInt <- contInt + 1 # updates number of iterations count
        if (contInt == nIter) { break } # if the number of iterations reached nIter stop EM algorithm
    }

    if (oscillations) { # Calculates and returns the oscillations of the updated model
        init$oscillations <- oscillation(initMTD)
    }
    if (iter) { # Returns the number of iterations and their log likelihood distances
        init$iterations <- contInt
        init$distlogL <- distlogL
    }
    init
}


