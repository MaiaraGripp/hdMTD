#' EM estimation of MTD parameters
#'
#' Estimation of MTD parameters through EM algorithm
#'
#' @param X A MTD chain sample
#' @param S The relevant lag set
#' @param M A stopping point for the EM algorithm. If NULL the algorithm will run
#' for a total of nIter iteractions.
#' @param init A list with initial parameters.
#' @param iter If True, returns the number of iterations of the
#' algorithm, that is, the number of times the initial parameters
#' were updated.
#' @param nIter An integer positive number with the maximum number of iterations.
#' @param oscillations If TRUE the function will return the estimated oscillations
#' for the updated model along with the estimated parameters.
#' @param A The states space
#'
#' @details Regarding the 'M' parameter: it functions as a stopping
#'criterion within the EM algorithm. When the difference between
#' the log-likelihood computed with the newly estimated parameters
#'and that computed with the previous parameters falls below M,
#'the algorithm halts. Nevertheless, if the value of nIter
#' (which represents the maximum number of iterations) is smaller
#'than the number of iterations required to meet the M criterion,
#'the algorithm will conclude its execution when nIter is reached.
#'To ensure that the M criterion is effectively utilized, we
#'recommend using a higher value for nIter, which is set to a
#' default of 100.
#'
#'Concerning the 'init' parameter, it is expected to be a list
#'comprising either 2 or 3 entries. These entries consist of:
#'an optional vector named 'p0' (representing an independent
#'distribution, ordered from the smallest to the greatest
#'element of A), a required list of matrices 'p_j' (a list containing
#'a stochastic matrix for each element of S, also ordered from the
#'smallest to the greatest element of S), and a vector named 'lambdas'
#'(representing the weights, first the weight for p0, and
#' then for each element in p_j).
#' @export
#' @return A list with the estimated parameters of the MTD model
#'
MTDest <- function(X,S,M=0.01,init,iter=FALSE,nIter=100,A=NULL,oscillations=FALSE){
  if(length(S) < 1  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){
    stop("S must be informed. S should be a number or a vector of positive integer numbers representing the relevant lags.")
  }
  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)

  X <- checkSample(X)
  lenX <- length(X)

  if(length(A)==0){
    A <- unique(X)
  }
  if( length(A)<=1   ||
      any(A%%1!=0)    )stop("The informed states space A must be a numeric vector with at least two integers.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it must include all states that occur in the sample.")
  }
  A <- sort(A)
  lenA <- length(A)

  if(!is.list(init)){stop("init must be a list with the initial parameters for the EM algorithm.")}
  if(!all(names(init) %in% c("p0","p_j","lambdas"))){
    stop("The init parameters must be names p0, p_j and lambdas, and at least p_j and lambdas must be informed.")
  }
  indep <- TRUE
  if(length(init$p0)==0){
    if( length(init$lambdas)==(lenS+1) & init$lambdas[1]!=0 ){
      stop("You didn't provide a distribution p0, but you entered a positive weight for this distribution. If your MTD model doesn't have an independent distribution, either set lambdas[1]=0 or provide a vector of lambdas with the same number of elements as S.")}
      indep <- FALSE
      init$p0 <- rep(0,lenA)
      if(length(init$lambdas)==lenS){
        init$lambdas <- c(0,init$lambdas)
      }
  }
  if( sum(init$p0)==0 ){indep <- FALSE}
  checkMTD(MTDmodel(S,A,w0=init$lambdas[1],
                    w_j = init$lambdas[-1],
                    p_j=init$p_j,
                    p0=init$p0)
           )
  if (!is.logical(iter)) {
    stop("Iter must be logical.")
  }
  if (!is.logical(oscillations)) {
    stop("oscillations must be logical.")
  }
  if(length(nIter)!=1 ||
     nIter %% 1 !=0 ||
     nIter <= 0 ){stop("nIter must be a positive integer number.")}
  if(length(M)!=0){
    if(length(M)!=1 ||
      !is.numeric(M) ||
       M<=0 ){stop("M is either NULL or a positive real number.")}
  }

  S0 <- c(S,0)
  lenS0 <- length(S0)
  base <- countsTab(X,S[1])
  baseSja <- freqTab(S,j=NULL,A,base)
  contInt <- 0
  distlogL <- NULL

  repeat{
    initMTD <- MTDmodel(S,A,w0=init$lambdas[1],
                        w_j = init$lambdas[-1],
                        p_j=init$p_j,
                        p0=init$p0)
    initLogL <- sum( log(as.vector(t(initMTD$P)))*baseSja$Nxa_Sj )

    #preciso separar o caso de não ter distribuição indep p0!!!

    #PASSO E
    pSja <- matrix(rep(rev(init$lambdas),lenA^(lenS0)),byrow = T,ncol=(lenS+1))
    colnames(pSja) <- S0
    pSja[,lenS0] <- pSja[,lenS0]*init$p0 #init$p0 must be ordered according to A

    indexA <- expand.grid(rep(list(1:lenA),lenS0))[,order((lenS0):1)]
    cont <- lenS
    for (i in 1:lenS) { #col
      pj <- init$p_j[[cont]]
      for (j in 1:lenA^(lenS0)) { #row
        pSja[j,i] <- pSja[j,i]*pj[indexA[j,i],indexA[j,(lenS0)]]
      }
      cont <- cont-1
    }
    colnames(pSja) <- paste0("lamXp_",S0)
    norm <- apply(pSja, 1, sum)
    Pj_xa_S <- pSja/norm
    colnames(Pj_xa_S) <- paste0("P",S0,"xa_S")
    # Fim do passo E

    # Passo M
    NxaXPjxa <- Pj_xa_S*baseSja$Nxa_Sj
    colnames(NxaXPjxa) <- paste0("NXP",S0,"xa")
    #eq (14)
    end_lambdas <- apply(NxaXPjxa,2,sum)/(lenX-S[1])
    names(end_lambdas) <- paste0("lam-",S0)
    end_lambdas <- rev(end_lambdas)
    #eq (15) para p0
    NPSja <- cbind(baseSja[,1:lenS0],NxaXPjxa)
    end_p0 <- rep(sum(NPSja$NXP0xa),lenA)
    if(indep){
      for (i in 1:lenA) {
        end_p0[i] <- sum((NPSja %>%
                            dplyr::filter(a==A[i]))$NXP0xa)/end_p0[i]
      }
    }
    names(end_p0) <- paste0("p_0(",A,")")
    #eq (15) para pj
    end_pj <- list()
    for (j in 1:lenS) {
      aux_pj <- matrix(0,ncol = lenA,nrow = lenA)
      for (i in 1:lenA) {
        aux_NPSja <- NPSja %>%
          dplyr::filter(NPSja[,j]==A[i]) #fix past x_j=A[i]
        aux_pj[i,] <- sum( aux_NPSja %>%
                             dplyr::select_at(paste0("NXP",S[j],"xa"))  )
        for (k in 1:lenA) {
          aux_pj[i,k] <- sum(aux_NPSja %>%
                               dplyr::filter(a==A[k]) %>%
                               dplyr::select_at(paste0("NXP",S[j],"xa"))
          )/aux_pj[i,k]
        }
      }
      colnames(aux_pj)=rownames(aux_pj) <- A
      end_pj[[lenS-j+1]] <- aux_pj
    }
    names(end_pj) <- paste0("p_-",rev(S))
    #Estimar a verossimilhança:
    endMTD <- MTDmodel(S,A,w0=end_lambdas[1],
                       w_j = end_lambdas[-1],
                       p_j=end_pj,
                       p0=end_p0)
    endlist <- list("lambdas"=end_lambdas,"p_j"=end_pj,"p0"=end_p0)

    endLogL <- sum( log( as.vector(t(endMTD$P)) )*baseSja$Nxa_Sj )

    distlogL[contInt+1] <- abs(endLogL-initLogL)
    if(length(M)==1){
      if(distlogL[contInt+1]<M){break}
    }
    contInt <- contInt+1
    init <- endlist
    if (contInt==nIter) {break}
  }
  if(oscillations){
    init$oscillations <- oscillation(initMTD)
  }
  if(iter){
    init$iterations <- contInt
    init$distlogL <- distlogL
    }
  init
}

