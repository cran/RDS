#' Extract or replace the *ult*imate (last) element of a vector or a list, or an element counting from the end.
#'
#' @param x a vector or a list.
#' @param i index from the end of the list to extract or replace (where 1 is the last element, 2 is the penultimate element, etc.).
#'
#' @return An element of `x`.
#'
#' @examples
#' x <- 1:5
#' (last <- ult(x))
#' (penultimate <- ult(x, 2)) # 2nd last.
#'
#' \dontshow{
#' stopifnot(last==5)
#' stopifnot(penultimate==4)
#' }
#'
#' @export
ult <- function(x, i=1L){
  x[[length(x)-i+1L]]
}

#' Set the class of the control list
#' 
#' This function sets the class of the control list, with the default being the
#' name of the calling function.
#' 
#' 
#' @param myname Name of the class to set.
#' @param control Control list. Defaults to the \code{control} variable in the
#' calling function.
#' @return The control list with class set.
#' @seealso check.control.class, print.control.list
#' @keywords utilities
#' @export
set.control.class <- function(myname=as.character(RDS::ult(sys.calls(),2)[[1L]]), control=get("control",pos=parent.frame())){
  class(control) <- c(myname, "control.list", "list")
  control
}

# Disable partial matching in control lists.
`$.control.list` <- getElement
