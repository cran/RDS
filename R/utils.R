
#' converts to character with minimal loss of precision for numeric variables
#' @param x the value
#' @param ... passed to either format or as.character.
as.char <- function(x,...){
	if(is.numeric(x))
		if(is.integer(x))
			ifelse(is.na(x),NA,format(x,trim=TRUE,scientific=FALSE,...))
		else
			ifelse(is.na(x),NA,format(x,trim=TRUE,scientific=FALSE,digits=15,...))
	else
		as.character(x,...)
}



