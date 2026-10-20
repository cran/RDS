

#' Homophily Configuration Graph Estimates
#' 
#' This function computes the Homophily Configuration Graph type estimates for a categorical variable.
#' 
#' 
#' @aliases RDS.HCG.estimates
#' @param rds.data An \code{rds.data.frame} with recruitment time set.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical variable to be analyzed.
#' @param subset An optional criterion to subset \code{rds.data} by. It is
#' an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{subset = seed > 0} to
#' exclude seeds. It can also be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param small.fraction Should a small sample fraction be assumed
#' @param N Population size. If NULL, this value is taken to be the population.size.mid 
#' attribute of the data and if that is not set, no finite population
#' correction is used.
#' @param to.factor force variable to be a factor
#' @param cont.breaks If variable is numeric, how many discretization points should be used in the calculation of the weights.
#' @return 
#' An object of class \code{rds.weighted.estimate}.
#' 
#' @author Ian E. Fellows
#' @seealso \code{\link{RDS.I.estimates}}, \code{\link{RDS.II.estimates}}, \code{\link{RDS.SS.estimates}}
#' @examples
#' 
#' data(fauxtime)
#' RDS.HCG.estimates(rds.data=fauxtime,outcome.variable='var1')
#' @export 
RDS.HCG.estimates <- function(rds.data, outcome.variable, N=NULL,subset=NULL, small.fraction=FALSE,
                            to.factor=FALSE, cont.breaks = 3){
  if(!has.recruitment.time(rds.data) && !small.fraction)
    stop("rds.data does not have recruitment time set")
  se <- substitute(subset)
  subset <- eval(se,rds.data,parent.frame())
  weight.type <- "HCG"
  if(length(outcome.variable) == 1){
    result <- RDS.estimates.local(rds.data,outcome.variable,
                                  subset=subset, weight.type=weight.type, N=N, small.fraction=small.fraction,to.factor=to.factor,cont.breaks=cont.breaks)
  }
  else {
    result <- lapply(outcome.variable,function(g){
      RDS.estimates.local(rds.data,g,subset=subset, weight.type=weight.type, N=N, 
                          small.fraction=small.fraction,to.factor=to.factor,cont.breaks=cont.breaks)
    })
    names(result) <- outcome.variable
    
  }
  return(result)
}
