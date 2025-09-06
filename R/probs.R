#' Estimated transition probabilities
#'
#' Computes the Maximum Likelihood estimators (MLE) for an MTD Markov chain with
#' relevant lag set \code{S}.
#'
#' @param X A vector or single-column data frame containing a sample of a Markov chain (`X[1]` is the most recent).
#' @param S A numeric vector of unique positive integers. Typically, \code{S} represents
#' a set of relevant lags.
#' @param matrixform Logical. If \code{TRUE}, the output is formatted as a stochastic
#' transition matrix.
#' @param A A numeric vector of distinct integers representing the state space.
#' If not provided, this function will set \code{A <- sort(unique(X))}.
#' @param warning Logical. If \code{TRUE}, the function warns the user when the state
#' space is automatically set as \code{A <- sort(unique(X))}.
#'
#' @return A data frame or a matrix containing estimated transition probabilities:
#'
#' - If \code{matrixform = FALSE}, the function returns a data frame with three columns:
#'   - The past sequence \eqn{x_S} (a concatenation of past states).
#'   - The current state \eqn{a}.
#'   - The estimated probability \eqn{\hat{p}(a | x_S)}.
#'
#' - If \code{matrixform = TRUE}, the function returns a stochastic transition matrix,
#'   where rows correspond to past sequences \eqn{x_S} and columns correspond to states in \code{A}.
#'
#' @details
#' The probabilities are estimated as:
#' \deqn{\hat{p}(a | x_S) = \frac{N(x_S a)}{N(x_S)}}
#' where \eqn{N(x_S a)} is the number of times the sequence \eqn{x_S} appeared in the sample
#' followed by \eqn{a}, and \eqn{N(x_S)} is the number of times \eqn{x_S} appeared
#' (followed by any state). If \eqn{N(x_S) = 0}, the probability is set to \eqn{1 / |A|}
#' (assuming a uniform distribution over \code{A}).
#'
#' @export
#'
#' @examples
#' # Simulate a chain from an MTD model
#' set.seed(1)
#' M <- MTDmodel(Lambda = c(1, 15, 30), A = c(1, 2, 3), lam0 = 0.05)
#' X <- perfectSample(M, N = 400)
#'
#' # Estimate probabilities for different subsets S
#' probs(X, S = c(1, 30))
#' probs(X, S = c(1, 15, 30), matrixform = TRUE)
#'
probs <- function(X, S, matrixform = FALSE, A = NULL, warning = FALSE){

    X <- checkSample(X) # Validate and preprocess the input sample
    check_probs_inputs(X, S, matrixform, A, warning) # Validate input parameters

    # Set the state space if not provided
    if (length(A) == 0) {A <- sort(unique(X))} else {A <- sort(A)}

    S <- sort(S, decreasing = TRUE) # Ensure S is sorted in decreasing order
    lenS <- length(S)

    # Compute frequency tables
    base <- countsTab(X, max(S))
    base <- freqTab(S = S, A = A, countsTab = base, complete = TRUE)

    # Construct output data frame
    probs <- data.frame(apply(base[, seq_len(lenS)], 1, paste0, collapse = ""),
                        base[, lenS + 1], base$qax_Sj)
    # probs = | Concatenated past states | Current state | Estimated probability |
    names(probs) <- c(paste("past_{", paste0(-S, collapse = ","),"}"), "a", "p(a|past)")

    # Convert output to stochastic matrix if requested
    if (matrixform) {
        Pest <- probs$`p(a|past)`
        dim(Pest) <- c(length(A), length(A)^lenS)
        Pest <- t(Pest)
        colnames(Pest) <- A
        rownames(Pest) <- unique(probs[, 1])
        probs <- Pest
    }
    return(probs)
}
