#' Makes a table with certain sequences in sample and their frequencies.
#'
#' @param S A set of lags in \eqn{1:d}.
#' @param j A lag \eqn{j \in (1:d) \setminus S}.
#' @param A The states space.
#' @param base A table with sequences of size d+1 in sample and their absolute frequency.
#' @param complete logical, if TRUE and there are sequences that didn't appear in the sample,
#'the base will be completed with such sequences appearing 0 times.
#'
#' @return A table with sequences x indexed by \eqn{S\cup j} followed by any element a, their frequencies and the MLE estimator of the transition probabilities.
#' @export
#' @importFrom dplyr %>%
#'
base_Sja <- function(S,j=NULL,A,base,complete=TRUE){

  #faz um teste para verificar se é possível contruir a maior base (l) acho melhor
  #fazer isso dentro do loop do FS
  #xa=try(expand.grid(rep(list(A),lenSj+1)),silent = TRUE) #change lenSj to lS+1
  #if(class(xa)=="try-error"){stop(paste0("For dim_S=",lenSj-1," the data set with all x_Sja sequences is too large, try a lower lS."))}
  A <- sort(A)
  Nx_Sj <- Nxa <- Nxa_Sj <- NULL #to solve problem of visible binding.
  Sj <- sort(c(S,j),decreasing = TRUE)
  if ( length(Sj) == 0 ) { stop("Set S U{j} can't be empty.") }
  lenSj <- length(Sj)
  filtrs <- c(paste0("x",Sj),"a")
  lenA <- length(A)

  base <- base %>%
            dplyr::group_by_at(filtrs) %>%
            dplyr::summarise(Nxa_Sj=sum(Nxa), .groups="drop")

  if ( ( nrow(base) < lenA^(lenSj+1) ) && complete ){ # if TRUE: this means that there are some
#sequences that didn't appear in the sample, so base needs to be completed.
    Tablexa <- expand.grid(rep(list(A),lenSj+1))[,order((lenSj+1):1)]
    list1 <- apply( base[,1:(lenSj+1)],1,paste0,collapse="" ) #list of sequencies in base
    list2 <- apply( Tablexa, 1, paste0,collapse="" ) #list of sequencies in Tablexa
    Tablexa <- Tablexa[ match(setdiff(list2,list1),list2),] #list of sequencies in Tablexa that where not in base
    Tablexa <- cbind(Tablexa,0) #attaching number of occurrences of these missing sequencies
    colnames(Tablexa) <- colnames(base)
    base <- rbind(base,Tablexa)
    base <- dplyr::arrange_at(base, filtrs)
  }

  base <- base %>% dplyr::group_by_at(paste0("x",Sj)) %>%
                   dplyr::mutate(Nx_Sj=sum(Nxa_Sj))
  base <- base %>% dplyr::mutate(qax_Sj=dplyr::if_else(Nx_Sj>0,Nxa_Sj/Nx_Sj,1/lenA)) %>%
                   dplyr::ungroup()
  base
}
