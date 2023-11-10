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
#'
#' @export
#'
#' @return A list with the estimated parameters of the MTD model
#' @examples
#' MTD <- MTDmodel(Lambda=c(1,10),A=c(0,1),w0=0.01)
#'X <- perfectSample(MTD,N=2000)
#'init <- list('p0'=c(0.4,0.6),'lambdas'=c(0.05,0.45,0.5),
#''p_j'=list(matrix(c(0.2,0.8,0.45,0.55),byrow = TRUE,ncol=2),
#'matrix(c(0.25,0.75,0.3,0.7),byrow = TRUE,ncol=2)))
#'MTDest(X,S=c(1,10),M=1,init)
#'MTDest(X,S=c(1,10),init=init,iter = TRUE)
#'MTDest(X,S=c(1,10),init=init,iter = TRUE,nIter=5)
#'MTDest(X,S=c(1,10),init=init,oscillations = TRUE)
MTDest <- function(X,S,M=0.01,init,iter=FALSE,nIter=100,A=NULL,oscillations=FALSE){
  if(length(S) < 1  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){
    stop("S must be informed. S should be a number or a vector of positive integer numbers representing the relevant lags.")
  }
  rS <- sort(S)
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
    stop("The init list entrances must be labeled 'p0', 'p_j', and 'lambdas', and at least 'p_j' and 'lambdas' must be provided.")
  }

  if(length(init$lambdas)!=(lenS+1)||
     !is.vector(init$lambdas)      ||
     !is.numeric(init$lambdas)     ||
     !all(init$lambdas>=0)         ||
     round(sum(init$lambdas),3)!=1)stop(paste0("The parameter init$lambdas must be a numeric, non-negative vector of length ", lenS+1," that must sum up to 1."))

  indep <- TRUE
  if(length(init$p0)==0 & init$lambdas[1]!=0){
    stop("You didn't provide a distribution p0, but you entered a positive weight for this distribution. If your MTD model doesn't have an independent distribution, set lambdas[1]=0.")
  }
  if(length(init$p0)==0 & init$lambdas[1]==0){
    indep <- FALSE
    init$p0 <- rep(0,lenA)
  }
  if( !length(init$p0) %in% c(1,lenA) ||
      !sum(init$p0) %in% c(0,1)      ||
      !is.numeric(init$p0)           ||
      !is.vector(init$p0)            ||
      !all(init$p0>=0)               ){
    stop("The parameter init$p0 must either 0 or a non-negative vector that sums up to 1.")
  }
  if( sum(init$p0)==0 ){indep <- FALSE}

  initMTD <- MTDmodel(rS,A,w0=init$lambdas[1],
                      w_j = init$lambdas[-1],
                      p_j=init$p_j,
                      p0=init$p0)
  checkMTD(initMTD)

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

  Pinit <- t(initMTD$P) #check if initial probabilities are compatible with initial sample
  dim(Pinit) <- NULL
  if(any(round(Pinit,5)==0)){
    if( baseSja$qax_Sj[which(round(Pinit,5)==0)]>0 ){
      stop("The initial parameters aren't compatible with the sample. There are sequences that appear in the sample for which the provided initial probability is zero. Run: 'MTDmodel(Lambda=S,A=sort(unique(X)),w0=init$lambdas[1],w_j=init$lambdas[-1],p_j = p_j, p0=init$p0)$P' and verify if there are null entries.")
    }
  }

  contInt <- 0
  distlogL <- NULL


  repeat{
    initMTD <- MTDmodel(rS,A,w0=init$lambdas[1],
                        w_j = init$lambdas[-1],
                        p_j=init$p_j,
                        p0=init$p0)

    if(any(baseSja$Nxa_Sj==0)){ #if TRUE is possible that some probabilities in
      #the initial P matrix are 0. So the loglikelihood will be -inf*0=Nan. The condition
      #in Nxa=0 is better since Nxa=0 is necessary for P=0. And P=0 with Nxa>0 should trigger error.
      initLogL <- sum( prodinf( log(as.vector(t(initMTD$P))) , baseSja$Nxa_Sj ) )
    }else{
      initLogL <- sum( log(as.vector(t(initMTD$P)))*baseSja$Nxa_Sj )
    }



    #step E
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
    #each element in pSja is a numerator of eq (13) for a different
    #past i^0_m (rows) and g={0,Lambda} columns
    colnames(pSja) <- paste0("lamXp_",S0)
    norm <- apply(pSja, 1, sum)
    norm[which(norm==0)] <- 1 #avoid 0/0

    Pj_xa_S <- pSja/norm #eq(13) for past i^0_m (rows) and g={0,Lambda} columns
    colnames(Pj_xa_S) <- paste0("P",S0,"xa_S")
    rownames(Pj_xa_S) <- apply(expand.grid(rep(list(A),lenS0))[,order((lenS0):1)],
                               1, paste0,collapse = "")
    #end step E

    #Step M
    NxaXPjxa <- Pj_xa_S*baseSja$Nxa_Sj
    colnames(NxaXPjxa) <- paste0("NXP",S0,"xa")
    #eq (14)
    end_lambdas <- apply(NxaXPjxa,2,sum)/(lenX-S[1])
    names(end_lambdas) <- paste0("lam-",S0)
    end_lambdas <- rev(end_lambdas)
    #eq (15) for p0
    NPSja <- cbind(baseSja[,1:lenS0],NxaXPjxa)
    end_p0 <- rep(sum(NPSja$NXP0xa),lenA)
    if(indep){
      for (i in 1:lenA) {
        end_p0[i] <- sum((NPSja %>%
                            dplyr::filter(a==A[i]))$NXP0xa)/end_p0[i]
      }
    }
    names(end_p0) <- paste0("p_0(",A,")")
    #eq (15) for pj
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
    #likelihoods
    endMTD <- MTDmodel(rS,A,w0=end_lambdas[1],
                       w_j = end_lambdas[-1],
                       p_j=end_pj,
                       p0=end_p0)
    endlist <- list("lambdas"=end_lambdas,"p_j"=end_pj,"p0"=end_p0)


    if(any(baseSja$Nxa_Sj==0)){ #if TRUE is possible that some probabilities in
      #the initial P matrix are 0. So the loglikelihood will be -inf*0=Nan. The condition
      #in Nxa=0 is better since Nxa=0 is necessary for P=0. And P=0 with Nxa>0 should trigger error.
      endLogL <- sum( prodinf( log(as.vector(t(endMTD$P))) , baseSja$Nxa_Sj ) )
    }else{
      endLogL <- sum( log(as.vector(t(endMTD$P)))*baseSja$Nxa_Sj )
    }

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

prodinf <- function(x,y){
  prinf <- numeric(length(x))
  for (i in 1:length(x)) {
    if (is.infinite(x[i]) && y[i] == 0) {
      prinf[i] <- 0
    } else {
      prinf[i] <- x[i] * y[i]
    }
  }
  prinf
}
