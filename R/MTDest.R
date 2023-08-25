#' EM estimation
#'
#' Estimation of MTD parameters through EM algorithm
#'
#' @param X A MTD chain sample
#' @param S The relevant lag set
#' @param M A stopping point for the EM algorithm.
#' @param init A list with initial parameters.
#' @param A The states space
#' @param indep_part If False the MTD model does not have
#' parameters independent from the past and p0=0.
#'
#' @details About the parameter M: it serves as a stopping
#'criterion for the EM algorithm. If the distance between
#'the log-likelihood calculated with the newly estimated
#'parameters and the log-likelihood calculated with the
#'previous parameters is less than M, the algorithm stops.
#'
#'About the parameter init: It's a list that must have 3
#'entries: a vector named p0 (representing an independent
#'distribution, ordered from the smallest to the greatest
#'element of A), p_j (a list containing a stochastic matrix
#'for each element of S, also ordered from the smallest to
#'the greatest element of S), and a vector named lambdas
#'(representing the weights, first the weight for p0, and
#' then for each element in p_j).
#'
#' @return A list with the estimated parameters of the MTD model
#' @export
#'
MTDest <- function(X,S,M=2,init,A=NULL,indep_part=TRUE){
  if(length(S) < 1  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){
    stop("S must be informed. S should be a number or a vector of
         positive integer numbers representing the relevant lags.")
  }
  X <- checkSample(X)
  if(length(A)==0){
    A <- unique(X)
  }
  if( length(A)<=1   ||
      any(A%%1!=0)    )stop("The informed states space A must be a numeric
                              vector with at least two integers.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it must include all states that occur in the sample.")
  }
  A <- sort(A)
  lenA <- length(A)
  lenX <- length(X)
  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)
  S0 <- c(S,0)
  lenS0 <- length(S0)
  base <- countsTab(X,S[1])
  baseSja <- freqTab(S,j=NULL,A,base)

# init is a list with the initial parameters
# init must have 3 entries, p0, p_j and lambdas
  #must check restriction for init

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
      for (i in 1:lenA) {
        end_p0[i] <- sum((NPSja %>%
                            dplyr::filter(a==A[i]))$NXP0xa)/end_p0[i]
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
      endlist <- list("lambdas"=end_lambdas,"p_j"=end_pj,"p0"=end_p0)
      #Estimar a verossimilhança:
      endMTD <- MTDmodel(S,A,w0=end_lambdas[1],
                         w_j = end_lambdas[-1],
                         p_j=end_pj,
                         p0=end_p0)

      endLogL <- sum( log( as.vector(t(endMTD$P)) )*baseSja$Nxa_Sj )

      if(abs(endLogL-initLogL)<=M){break}
      init <- endlist
  }
  endlist
}
require(devtools)
use_git()
