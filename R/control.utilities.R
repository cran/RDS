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

#' Named element accessor for ergm control lists
#' 
#' Utility method that overrides the standard `$' list accessor to disable
#' partial matching for ergm \code{control.list} objects
#' 
#' Executes \code{\link[base]{getElement}} instead of \code{\link[base]{$}} so
#' that element names must match exactly to be returned and partially matching
#' names will not return the wrong object.
#' 
#' @param object list-coearceable object with elements to be searched
#' @param name literal character name of list element to search for and return
#' @return Returns the named list element exactly matching \code{name}, or
#' \code{NULL} if no matching elements found
#' @author Pavel N. Krivitsky
#' @seealso see \code{\link{getElement}}
#' @name control.list.accessor
#' @export
`$.control.list` <- function(object, name) object[[name, exact = TRUE]]
