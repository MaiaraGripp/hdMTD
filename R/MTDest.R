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
#' Concerning the \code{init} parameter, it is expected to be a list with 2 or
#' 3 entries. These entries consist of:
#' an optional vector named \code{p0}, representing an independent
#' distribution (the probability in the first entry of \code{p0} must be
#' that of the smallest element in \code{A} and so on), a required list
#' of matrices \code{pj}, containing a stochastic matrix for each
#' element of \code{S} ( the first matrix corresponds to the smallest lag in
#' \code{S} and so on), and a vector named \code{lambdas} representing
#' the weights, the first entry must be the weight for \code{p0}, and then one entry
#' for each element in \code{pj} list. If your MTD model does not have an independent
#' distribution \code{p0}, set \code{init$lambdas[1]=0}.
#'
#'
#' @references
#' Lebre, Sophie and Bourguignon, Pierre-Yves. (2008).
#' An EM algorithm for estimation in the Mixture Transition Distribution model.
#' *Journal of Statistical Computation and Simulation*, *78*(1), 1-15.
#' \doi{10.1080/00949650701266666}
#'
#'
#' @param X A vector or single-column data frame containing an MTD chain sample
#'  (`X[1]` is the most recent).
#' @param S A numeric vector of distinct positive integers, sorted increasingly.
#' Typically, \code{S} represents a set of relevant lags.
#' @param M A stopping point for the EM algorithm. If \code{M=NULL} the algorithm
#'  will run for a total of \code{nIter} iterations (i.e., no convergence check).
#' @param init A list with initial parameters: \code{p0} (optional), \code{lambdas}
#' (required), \code{pj} (required). The entries in \code{lambdas} are weights
#' for the distribution \code{p0} and the distributions present in the list
#' \code{pj}. Therefore, the order in which the elements appear in the vector
#' \code{lambdas} is important for correct assignment. See  *Details*.
#' @param iter Logical. If \code{TRUE} include iteration diagnostics in the output
#' (number of updates \code{iterations}; vector of absolute log-likelihood changes
#' \code{deltaLogLik}; and \code{lastComputedDelta}, the last delta log-likelihood
#' compared against \code{M}).
#' @param nIter An integer positive number with the maximum number of EM iterations.
#' @param A Optional integer vector giving the state space. If omitted, defaults
#' to \code{sort(unique(X))}.
#' @param oscillations Logical. If \code{TRUE}, also compute oscillations for the
#' fitted model (see \code{\link{oscillation}}).
#'
#' @export
#'
#' @return
#' An S3 object of class \code{"MTDest"} (a list) with at least the following elements:
#' \itemize{
#'   \item \code{lambdas}: estimated mixture weights (independent part first, if any).
#'   \item \code{pj}: list of estimated transition matrices \eqn{p_j}.
#'   \item \code{p0}: estimated independent distribution (if applicable).
#'   \item \code{S}: the lag set used for estimation.
#'   \item \code{A}: the state space used for estimation.
#'   \item \code{logLik}: log-likelihood of the final fitted model.
#'   \item \code{call}: the matched call.
#'   \item \code{oscillations} (optional): returned if \code{oscillations = TRUE}.
#'   \item \code{iterations}, \code{deltaLogLik}, \code{lastComputedDelta} (optional):
#'         returned if \code{iter = TRUE}. Here, \code{iterations} is the number of
#'         parameter updates performed; \code{deltaLogLik} stores the successive
#'         absolute log-likelihood changes for the accepted updates; and
#'         \code{lastComputedDelta} is the last computed change (which may be below
#'         \code{M} when the loop stops by convergence).
#' }
#'
#' @examples
#' # Simulating data.
#' set.seed(1)
#' MTD <- MTDmodel(Lambda = c(1, 10), A = c(0, 1), lam0 = 0.01)
#' X <- perfectSample(MTD, N = 2000)
#'
#' # Initial Parameters:
#' init <- list(
#'   'p0' = c(0.4, 0.6),
#'   'lambdas' = c(0.05, 0.45, 0.5),
#'   'pj' = list(
#'       matrix(c(0.2, 0.8, 0.45, 0.55), byrow = TRUE, ncol = 2),
#'       matrix(c(0.25, 0.75, 0.3, 0.7), byrow = TRUE, ncol = 2)
#'     )
#'  )
#'
#' fit <- MTDest(X, S = c(1, 10), init = init, iter = TRUE)
#' str(fit, max.level = 1)
#' fit$logLik
#'
#' fit2 <- MTDest(X, S = c(1, 10), init = init, oscillations = TRUE)
#' fit2$oscillations
#'
MTDest <- function(X, S, M = 0.01, init, iter = FALSE, nIter = 100, A = NULL, oscillations = FALSE) {

  # Validate inputs
  X <- checkSample(X)
  check_MTDest_inputs(X, S, M, init, iter, nIter, A, oscillations)

  if (is.null(A)) { A <- sort(unique(X)) } else { A <- sort(A) }
  lenA <- length(A)

  indep <- TRUE # Used for p0 update
  if (length(init$p0) == 0 || all(init$p0 == 0)) {
    if (init$lambdas[1] == 0 ) {
      indep <- FALSE
      init$p0 <- rep(0, lenA)
    }
  }
  # Creates an MTD object to validate inputs with checkMTD
  initMTD <- MTDmodel(Lambda = S, A = A, lam0 = init$lambdas[1],
                      lamj = init$lambdas[-1], pj = init$pj,
                      p0 = init$p0, indep_part = indep)
  checkMTD(initMTD)

  rS <- S
  S <- sort(S, decreasing = TRUE)
  lenS <- length(S)
  lenX <- length(X)

  S0 <- c(S, 0)
  lenS0 <- lenS + 1
  base <- countsTab(X, S[1])
  baseSja <- freqTab(S, j = NULL, A, base)
  indexA <- expand.grid(rep(list(seq_len(lenA)), lenS0))[, order(lenS0:1)]

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

  ## --- EM algorithm

  contInt <- 0
  deltaLogLik <- numeric(nIter)

  # Parameter update loop. It terminates if the distance between the
  # log-likelihood of the initial and updated models is less than M, or if the
  # maximum number of iterations (nIter) is reached.
  repeat{

    initMTD <- suppressWarnings(
      MTDmodel(rS, A, lam0 = init$lambdas[1], lamj = init$lambdas[-1],
               pj = init$pj, p0 = init$p0))

    # Compute log-likelihood
    if(any(baseSja$Nxa_Sj == 0)){
      # If TRUE is possible that some probabilities in initMTD$P are also 0, and
      # log(P(a|x))*N(xa) = log(0)*0 = -inf*0 = Nan. To prevent this use the auxiliary
      # function prodinf (see utils.R). With it, terms resulting from -inf*0 are set to 0.
      initLogLik <- sum(prodinf(log(as.vector(t(initMTD$P))), baseSja$Nxa_Sj))
    }else{
      initLogLik <- sum(log(as.vector(t(initMTD$P))) * baseSja$Nxa_Sj)
    }

    # - Step E (Expectation)

    # Creates a stochastic matrix with the lambdas
    pSja <- matrix(rep(rev(init$lambdas), lenA^lenS0), byrow = T, ncol = lenS0)

    # pSja column referent to lam0 is multiplied by p0
    pSja[, lenS0] <- pSja[, lenS0] * init$p0
    # pSja column referent to lamj is multiplied by pj, j\in S
    cont <- lenS
    for (i in seq_len(lenS)) { #col
      pj <- init$pj[[cont]]
      for (j in seq_len(lenA ^ lenS0)) { #row
        pSja[j, i] <- pSja[j, i] * pj[indexA[j, i], indexA[j, (lenS0)]]
      }
      cont <- cont-1
    }
    # Normalize pSja so each row is either a distribution or a vector of 0s.
    norm <- rowSums(pSja)
    norm[which(norm == 0)] <- 1 # to avoid 0/0
    Pj_xa_S <- pSja/norm
    # Pj_xa_S is a matrix with the conditional probabilities of a hidden variable Z
    # (a state variable, i.e., Z\in S\cup{0}), given any possible size lenS0
    # sequence with elements from A.

    # - End of step E

    # - Step M (Maximization)

    # Updates
    # lambdas update
    NxaXPjxa <- Pj_xa_S * baseSja$Nxa_Sj
    colnames(NxaXPjxa) <- paste0("NXP", S0, "xa")
    end_lambdas <- apply(NxaXPjxa, 2, sum)/(lenX-S[1])
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
      end_pj[[lenS - j + 1]] <- aux_pj
    }
    # end_pj is the list with updated pj matrices

    # Computations
    # Updated MTD model
    endMTD <- suppressWarnings(
      MTDmodel(rS, A, lam0 = end_lambdas[1], lamj = end_lambdas[-1],
               pj = end_pj, p0 = end_p0))
    endlist <- list("lambdas" = endMTD$lambdas, "pj" = endMTD$pj, "p0" = endMTD$p0)

    # Updated log-likelihood
    if (any(baseSja$Nxa_Sj == 0)) {
      endLogLik <- sum( prodinf( log(as.vector(t(endMTD$P))), baseSja$Nxa_Sj ) )
    }else{
      endLogLik <- sum( log(as.vector(t(endMTD$P)))*baseSja$Nxa_Sj )
    }

    # Compute difference between updated and former log-likelihood
    lastComputedDelta <- abs(endLogLik - initLogLik)

    # Stop criteria 1 (convergence)
    if (!is.null(M)){
      if(lastComputedDelta < M){ break }
    }

    # If lastComputedDelta >= M:
    init <- endlist # initial parameters are replaced with updated parameters
    contInt <- contInt + 1 # updates number of iterations count
    deltaLogLik[contInt] <- lastComputedDelta # Stores deltaLogLik

    # Stop criteria 2 (max number of iteration)
    if (contInt == nIter) { break }
  }
  # Ensuring pj row and colnames are ok
  for (i in seq_len(lenS)) {
    colnames(init$pj[[i]]) = rownames(init$pj[[i]]) <- A
  }
  # Recomputes MTD with final parameters
  finalMTD <- MTDmodel(rS, A, lam0 = init$lambdas[1], lamj = init$lambdas[-1],
                       pj = init$pj, p0 = init$p0)
  # Outputs
  out <- list("lambdas" = finalMTD$lambdas, "pj" = finalMTD$pj, "p0" = finalMTD$p0)

  # Recomputes LogLik
  if (any(baseSja$Nxa_Sj == 0)) {
    out$logLik <- sum(prodinf(log(as.vector(t(finalMTD$P))), baseSja$Nxa_Sj))
  } else {
    out$logLik <- sum(log(as.vector(t(finalMTD$P))) * baseSja$Nxa_Sj)
  }
  # Optional outputs
  if (iter) {
    out$iterations <- contInt
    out$deltaLogLik <- if (contInt > 0) deltaLogLik[seq_len(contInt)] else numeric(0)
    out$lastComputedDelta <- lastComputedDelta
  }
  if (oscillations) {
    out$oscillations <- oscillation(finalMTD)
  }
  out$call <- match.call()
  out$S    <- rS
  out$A    <- A
  class(out) <- c("MTDest", "list")
  out
}

