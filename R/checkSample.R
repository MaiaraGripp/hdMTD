#' Checks if a sample is suitable for use
#'
#' @param X A vector with a Markov Chain.
#'
#' @return Returns possible sample errors.
#'
checkSample <- function(X){
  if ( !is.numeric(X) ) { stop("X must be a numeric vector.") }
  if ( length(ncol(X)) != 0 ) { stop("X must have only 1 dimension.") }
  if ( any(is.na(X)) ) { stop("NA values are not allowed in the sample.") }
  if ( length(unique(X)) == 1 ) { stop("All elements in the sample are the same.") }
}

#colocar mais restrições? Tamanho minimo de amostra...

#acrescentar input message=TRUE nos parâmetros, não esquecer de colocar @param message A logic argument to exhibit a message
# embaixo de @param X na documentação.
#if( length(levels(X)) > 6){warning("States space size is ",length(levels(X)), ", might be too big.")}
#if( message==TRUE & any( levels(X)!= 1:length(levels(X)) ) ){
#message("The elements ",paste0(levels(X)," ")," in sample where replaced by ",paste0(seq(1,length(levels(X)))," ")," respectively.")
#}
