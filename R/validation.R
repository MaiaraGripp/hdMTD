# validation.R - Input validation functions for the package
#
# This file contains helper functions for validating inputs in various functions
# within the package. These functions are not exported and are used internally
# to ensure that user inputs are correctly formatted and consistent.

# 1 - checkMTD()
# 2 - check_MTDmodel_inputs()
# 3 - check_freqTab_inputs()
# 4 - check_dTVsample_inputs()
# 5 - check_probs_inputs()
# 6 - check_hdMTD_FS_inputs()
# 7 - check_hdMTD_CUT_inputs()
# 8 - check_hdMTD_BIC_inputs()
# 9 - check_MTDest_inputs()

#########################################################################

checkMTD <- function(MTD){
  # Verifies if an object of class MTD is correctly structured, containing
  # all the necessary parameters, and if these parameters satisfy their
  # respective constraints.

  #Verifies if the object is a list with class MTD.
  if (!is.list(MTD)) {
    stop("MTD must be a list. Try using MTDmodel() to create an MTD.")
  }
  if (!inherits(MTD, "MTD")) {
    stop("MTD must be an object of class MTD. Try using MTDmodel() to create an MTD.")
  }

  # Checks if Lambda is a numeric vector of unique positive integers in ascending order
  if ( any(MTD$Lambda <= 0) || !all(MTD$Lambda%%1 == 0) || !is.vector(MTD$Lambda) ||
       !is.numeric(MTD$Lambda) || length(MTD$Lambda) != length(unique(MTD$Lambda)) ) {
    stop("Lambda must be a vector of unique positive integers. Try using MTDmodel() to create an MTD.")
  }
  if ( any(sort(MTD$Lambda) != MTD$Lambda) ) {
    stop("Lambda must be sorted in ascending order. Try using MTDmodel() to create an MTD.")
  }

  lenL <- length(MTD$Lambda)

  # Checks if A is a numeric vector of integers (length ≥ 2), sorted in ascending order
  if (length(MTD$A) <= 1 || !is.vector(MTD$A) || any(MTD$A%%1 != 0) || length(MTD$A) != length(unique(MTD$A))) {
    stop("State space A must be a vector (length >= 2) of unique integers. Try using MTDmodel() to create an MTD.")
  }
  if (any(sort(MTD$A) != MTD$A)) {
    stop("State space A must be sorted in ascending order. Try using MTDmodel() to create an MTD.")
  }

  lenA <- length(MTD$A)

  # Checks if p0 is a numeric nonnegative vector of length 1 or length(A), summing to 1
  if (!is.numeric(MTD$p0) || !is.vector(MTD$p0) || !all(MTD$p0 >= 0)) {
    stop("p0 must be a nonnegative numeric vector. Try using MTDmodel() to create an MTD.")
  }
  if (!length(MTD$p0) %in% c(1, lenA)) {
    stop(paste0("p0 must be either a scalar 0 or a numeric vector of length ",
        lenA, ". Try using MTDmodel() to create an MTD."))
  }
  if (round(sum(MTD$p0), 3) != 1 & sum(MTD$p0) != 0) {
    stop("The elements in p0 must either sum to 1 or all be 0. Try using MTDmodel() to create an MTD.")
  }

  # Checks if lambdas is a numeric nonnegative vector of length
  # (length(Lambda) + 1) that sums to 1
  if (!is.numeric(MTD$lambdas) || round(sum(MTD$lambdas), 3) != 1 ||
      !all(MTD$lambdas >= 0) || length(MTD$lambdas) != (lenL + 1)) {
    stop(paste0("lambdas must be a vector of length ", lenL + 1, " (the number of
  relevant lags in Lambda plus 1), consisting of nonnegative numbers that sum to 1.
  The first element of the lambdas vector is the weight for the independent
  distribution p0, if your MTD model does not include an independent distribution,
  set lambdas[1] to 0. Try using MTDmodel() to create an MTD."
    ))
  }

  # Checks if pj is a list with length(Lambda) elements, each containing a
  # stochastic matrix of size length(A) x length(A)
  if(!is.list(MTD$pj) || length(MTD$pj) != lenL ||
     !all(sapply(MTD$pj, is.matrix)) || !all(sapply(MTD$pj,dim) == c(lenA,lenA))) {
    stop(paste0("pj must be a list with ", lenL, " stochastic matrices ", lenA,
         "x",lenA,". Try using MTDmodel() to create an MTD."))
  }
  aux <- do.call(rbind, MTD$pj)
  if(!is.numeric(aux) || !all(round(apply(aux, 1, sum),3) == 1) || !all(aux>=0)) {
    stop(paste0("pj must be a list with ", lenL, " stochastic matrices ", lenA,
    "x",lenA,". In other words, each matrix row must sum up to 1. Try using MTDmodel() to create an MTD."))
  }
}
# Note: This package includes a function called MTDmodel(), which outputs
# a properly structured MTD object that does not require additional checks.
# However, since the user can create the MTD object manually, this checkMTD()
# is used, within the functions that use MTD objects as inputs, to prevent errors.


#########################################################################
#########################################################################
#########################################################################

check_MTDmodel_inputs <- function(Lambda, A, lam0, lamj, pj, p0, single_matrix, indep_part){
  # Validates the inputs in MTDmodel() function.

  if(!is.numeric(Lambda) || any(Lambda <= 0) || !all(Lambda %% 1 == 0) ||
     !is.vector(Lambda) || length(Lambda) != length(unique(Lambda)) ) {
    stop("Lambda must be a vector of unique positive integers.")
  }

  if( length(A) <= 1 || any(A %% 1 != 0) || length(A) != length(unique(A)) || any(A < 0) ) {
    stop("A must be a vector of nonnegative integers with length(A)>=2.")
  }

  if( !is.logical(indep_part) ) stop("Argument indep_part must be TRUE or FALSE.")

  if ( !is.null(p0) ) {
    if ( all(p0 >= 0) ) {
      if ( sum(p0) != 0 && round(sum(p0),5) != 1 ) {
        stop("p0 must add to 1 or be NULL or a vector of 0.")
        if ( round(sum(p0),5) == 1 && length(p0) != length(A) ) {
          stop("If a distribution p0 is provided it must add to 1 and have a probability
          for each element in A.")
        }
      }
    } else { stop( "p0 must be either NULL or a numeric vector of nonnegative numbers.") }
  }

  if ( !is.null(lam0) ) {
    if ( length(lam0) != 1 || !is.numeric(lam0) || lam0 < 0 || lam0 >= 1 ){
      stop("lam0 is either NULL or a number in the range [0,1).")
    }
  }

  if ( !is.null(lamj) ) {
    if( !is.numeric(lamj) || length(lamj) != length(Lambda) || !all(lamj > 0) || sum(lamj)>1 ) {
      stop(paste0("lamj must be NULL or a vector of positive numbers. The length of lamj must be ",length(Lambda),
                  " and sum(lamj) cannot be greater than 1."))
    }
  }

  if( !is.logical(single_matrix) ) stop("single_matrix must be TRUE or FALSE.")

  if ( !is.null(pj) ) {
    if ( !is.list(pj) ) stop("pj must be either NULL or a list of matrices.")
    if ( any( lapply(pj, is.matrix) == FALSE ) ) stop("pj must be either NULL or a list of matrices.")
    if ( single_matrix && length(pj) != 1 ) stop("Since single_matrix=TRUE, pj must be NULL or be a list with a single stochastic matrix.")
    if ( !single_matrix && length(pj) != length(Lambda) ) stop("pj must be NULL or be a list with ", length(Lambda), "matrices.")
    aux <- do.call(rbind,pj)
    if( !all(round(apply(aux, 1, sum), 3) == 1)  || !all(aux >= 0) || !is.numeric(aux) ||
        ncol(aux) != length(A) || any(sapply(pj, dim) != length(A)) ) {
      stop(paste0("pj must be a list of stochastic matrices ", length(A), "x",length(A)))
    }
  }
}

#########################################################################
#########################################################################
#########################################################################

check_freqTab_inputs <- function(S, j, A, countsTab, complete) {
  # Validates the inputs in freqTab() function.

  if(!is.data.frame(countsTab)) {
    stop("countsTab must be a tibble or a dataframe. Try using countsTab() function.")
  }
  if(!is.null(S) && length(S) > 0) {
    if(!is.vector(S) || any(S %% 1 != 0) || any(S < 1))
      stop("S must be a vector of positive integers or NULL.")
  }

  if(!is.null(j) && length(j) > 0) {
    if(length(j) != 1 || j %% 1 != 0 || j %in% S)
      stop("j must be a single integer not present in S.")
  }

  Sj <- c(S, j)
  if(length(Sj) == 0) {
    stop("The set {S}U{j} can't be NULL.")
  }
  d <- ncol(countsTab) - 2
  if(max(Sj) > d) {
    stop("The set {S}U{j} cannot have an element greater than d.")
  }
  if(!is.logical(complete)) {
    stop("Complete must be a logical argument.")
  }
  if(!all(unique(unlist(countsTab[, -(d+2)])) %in% A)) {
    stop("A must contain all elements that appear in the countsTab sequences.")
  }
  if(length(A) <= 1 || any(A %% 1 != 0) || length(A)!=length(unique(A)) || any(A < 0)  ) {
    stop("A must be a vector of length greater than 1 composed of unique nonnegative integers.")
  }
}

#########################################################################
#########################################################################
#########################################################################

check_dTVsample_inputs <- function(S, j, A, base, lenA, A_pairs, x_S) {
  # Validates the inputs in dTV_sample() function.

  if( length(S) > 0 ){
      if( !is.vector(S) || any(S%%1 != 0) || any(S<1) ) {
        stop("S must be a positive integer vector, a number or NULL.")
      }
      if( length(x_S)!=length(S) ) {
        stop("x_S must be a sequence of length(S) elements.")
      }
  } else {
      if( !is.null(S) ) {
        stop("S must be a positive integer vector, a number or NULL.")
      }
  }

  if( length(j)!=1 || j%%1!=0 || j %in% S || j < 1 ) {
    stop("j must be a integer number in the complement of S.")
    }

  if( is.null(A) ){
      if( length(lenA)==0 || length(A_pairs)==0 ) {
        stop("Either the state space A must be provided, or both lenA (the number of elements in A) and A_pairs (all possible pairs of elements from A) must be provided.")
      }
      if( !is.numeric(lenA) || lenA < 2 || lenA %% 1 != 0 ) {
        stop("lenA must be an integer number >= 2.")
      }
      if( !is.matrix(A_pairs) || ncol(A_pairs) != 2 || any(A_pairs%%1 != 0) ) {
        stop("A_pairs must be a matrix with two columns containing unique integer pairs.")
      }
      if( length(S)>0 && !all(x_S %in% A_pairs) ) {
        stop("x_S must be a sequence of elements from the state space A.")
      }
  } else {
      if( length(A) <= 1 || !is.vector(A) || any( A %% 1 != 0 ) || length(A) != length(unique(A)) ) {
        stop("A must be a vector of length greater than 1 composed of unique integers.")
      }
      if( length(lenA) != 0 || length(A_pairs) != 0 ) {
        warning("Since the state space A was provided, this function will set lenA <- length(A) and A_pairs <-  t(utils::combn(A, 2))}, even if you have provided at least one of them.")
      }
      if( length(S)>0 && !all(x_S %in% A) ) {
        stop("x_S must be a sequence of elements from A.")
      }
  }
}

#########################################################################
#########################################################################
#########################################################################


 check_probs_inputs <- function(X, S, matrixform, A, warning) {
   # Validates the inputs in probs() function.

   if( length(S) < 1 || !is.numeric(S) || any(S <= 0) || any( S%%1 != 0) || length(S) != length(unique(S)) ){
     stop("S must be a numeric vector of unique positive integers with length >= 1.")
   }

   if( length(A) > 0 ){
     if( length(A) <= 1 || any( A%%1 != 0 ) || length(A) != length(unique(A)) ) {
        stop("A must be a vector of distinct integers and length >=2.")
     }
     if ( !all( A %in% unique(X) ) ) {
        warning("Some elements in A do not appear in the sample.")
     }
     if ( !all( unique(X) %in% A ) ) {
        stop("The sample contains elements that do not appear in A.")
     }
   } else if (warning) {
     warning("States space A not provided. The function will set A <- sort(unique(X)).")
   }

   if(!is.logical(matrixform)){
     stop("matrixform should be either TRUE or FALSE.")
   }
 }

 #########################################################################
 #########################################################################
 #########################################################################

 check_hdMTD_FS_inputs <- function(X, d, l, A, elbowTest, warning) {
   # Validates the inputs in hdMTD_FS() function.

   if( length(d) != 1 || !is.numeric(d) || d < 2 || (d %% 1) != 0 ){
     stop("The order d must be an integer equal to or greater than 2.")
   }

   if( length(l) != 1 || !is.numeric(l) || l%%1 != 0 || l > d ) {
     stop("The l value is not valid for FS method. l should be a positive integer smaller or equal to d.")
   }

   if( length(A) > 0 ){
       if( length(A) <= 1 || any( A%%1 != 0 ) || length(A) != length(unique(A)) ) {
         stop("A must be a vector of distinct integers and length >=2.")
       }
       if ( !all( A %in% unique(X) ) ) {
         warning("Some elements in A do not appear in the sample.")
       }
       if ( !all( unique(X) %in% A ) ) {
         stop("The sample contains elements that do not appear in A.")
       }
   } else if (warning) {
      warning("States space A not provided. The function will set A <- sort(unique(X)).")
   }

   if(!is.logical(elbowTest)){
     stop("elbowTest should be either TRUE or FALSE.")
   }
 }


 #########################################################################
 #########################################################################
 #########################################################################

 check_hdMTD_CUT_inputs <- function(X, d, S, alpha, mu, xi, A, warning) {
   # Validates the inputs in hdMTD_CUT() function.

   if( length(d) != 1 || !is.numeric(d) || d < 2 || (d %% 1) != 0 ){
     stop("The order d must be an integer equal to or greater than 2.")
   }

   if(length(S) < 2 || !is.numeric(S) || any( S%%1 != 0) || max(S) > d ||
      length(S) != length(unique(S)) || any( S <= 0)) {
     stop("S must be a vector of distinct positive integers containing at least 2 values.")
   }

   if( length(A) > 0 ){
     if( length(A) <= 1 || any( A%%1 != 0 ) || length(A) != length(unique(A)) ) {
       stop("A must be a vector of distinct integers and length >=2.")
     }
     if ( !all( A %in% unique(X) ) ) {
       warning("Some elements in A do not appear in the sample.")
     }
     if ( !all( unique(X) %in% A ) ) {
       stop("The sample contains elements that do not appear in A.")
     }
   } else if (warning) {
     warning("States space A not provided. The function will set A <- sort(unique(X)).")
   }

   if ( is.na(alpha) || !is.numeric(alpha) || alpha <= 0 || length(alpha) != 1 ) {
     stop("alpha must be a positive real number.")
   }
   if ( is.na(mu) || !is.numeric(mu) || mu <= 0 || mu <= (exp(mu) - 1)/2 || length(mu) != 1 ) {
     stop("mu must be a positive real number smaller than (exp(mu)-1)/2.")
   }
   if ( is.na(xi) || !is.numeric(xi) || xi <= 0 || length(xi) != 1 ) {
     stop("xi must be a positive real number.")
   }
 }

 #########################################################################
 #########################################################################
 #########################################################################

 check_MTDest_inputs <- function(X, S, M, init, iter, nIter, A, oscillations) {
   # Validates the inputs of MTDest() function

   if(!is.numeric(S) || any(S <= 0) || !all(S %% 1 == 0) ||
      !is.vector(S) || length(S) != length(unique(S)) || any(sort(S) != S) ) {
     stop("S must be a vector of unique positive integers sorted from smallest to largest.")
   }

   if (!is.null(A)) {
     if( length(A) <= 1 || any( A%%1 != 0 ) || length(A) != length(unique(A)) ) {
       stop("A must be a vector of distinct integers and length >=2.")
     }
     if ( !all( A %in% unique(X) ) ) {
       warning("Some elements in A do not appear in the sample.")
     }
     if ( !all( unique(X) %in% A ) ) {
       stop("The sample contains elements that do not appear in A.")
     }
   }

   if(!is.list(init)){
     stop("init must be a list with the initial parameters for the EM algorithm.")
   }
   if(!all(names(init) %in% c("p0", "pj", "lambdas"))){
     stop("The init list entrances must be labeled 'p0', 'pj', and 'lambdas', and at least 'pj' and 'lambdas' must be provided.")
   }

   if( is.null(init$lambdas) || length(init$lambdas) != (length(S) + 1) || !all(init$lambdas >= 0)) {
     stop("The parameter init$lambdas must be a numeric, non-negative vector of length ", length(S)+1)
   }

   if (is.null(init$pj) || length(init$pj) != length(S) || !is.list(init$pj)) {
     stop("The parameter init$pj must be a list with ", length(S), " matrices.")
   }

   if( length(init$p0) == 0 && init$lambdas[1] > 0) {
     stop("You did not provide a distribution init$p0, but provided a positive weight for it (init$lambdas[1]>0).")
   }
   if( length(init$p0) != 0 && sum(init$p0) > 0 && init$lambdas[1] == 0) {
     stop("You have provided values for init$p0, but provided weight 0 for them (init$lambdas[1]=0).")
   }

   if (!is.logical(oscillations)) stop("oscillations must be TRUE or FALSE.")
   if (!is.logical(iter)) stop("Iter must be TRUE or FALSE.")

   if( length(nIter) != 1 || !is.numeric(nIter) || nIter %% 1 != 0 || nIter <= 0 ){
     stop("nIter must be a positive integer.")
     }
   if(!is.null(M)){
     if(length(M) != 1 || !is.numeric(M) || M <= 0 ){
       stop("M is either NULL or a positive real number.")
       }
   }
 }
