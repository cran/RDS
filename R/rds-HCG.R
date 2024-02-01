

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
#' @param empir.lik Should confidence intervals be estimated using 
#' empirical likelihood.
#' @param N Population size to be used to calculate the empirical likelihood interval. If NULL, this value is
#' taken to be the population.size.mid attribute of the data and if that is not set, no finite population
#' correction is used.
#' @param to.factor force variable to be a factor
#' @param cont.breaks If variable is numeric, how many discretization points should be used in the calculation of the weights.
#' @return If the \code{empir.lik} is true, an object of class
#' \code{rds.interval.estimate} is returned. This is a list with components
#' \itemize{ \item\code{estimate}: The numerical point estimate of proportion
#' of the \code{trait.variable}.  \item\code{interval}: A matrix with six
#' columns and one row per category of \code{trait.variable}: \itemize{
#' \item\code{point estimate}: The HT estimate of the population mean.
#' \item\code{95\% Lower Bound}: Lower 95\% confidence bound.  \item\code{95\%
#' Upper Bound}: Upper 95\% confidence bound.  \item\code{Design Effect}: The
#' design effect of the RDS.  \item\code{s.e.}: Standard error.  \item\code{n}:
#' Count of the number of sample values with that value of the trait.  } }
#' 
#' Otherwise an object of class \code{rds.HCG.estimate} object is returned.
#' 
#' @author Ian E. Fellows
#' @seealso \code{\link{RDS.I.estimates}}, \code{\link{RDS.II.estimates}}, \code{\link{RDS.SS.estimates}}
#' @examples
#' 
#' data(fauxtime)
#' RDS.HCG.estimates(rds.data=fauxtime,outcome.variable='var1')
#' @export 
RDS.HCG.estimates <- function(rds.data, outcome.variable, N=NULL,subset=NULL, small.fraction=FALSE,
                            empir.lik=TRUE, to.factor=FALSE, cont.breaks = 3){
  if(!has.recruitment.time(rds.data) && !small.fraction)
    stop("rds.data does not have recruitment time set")
  se <- substitute(subset)
  subset <- eval(se,rds.data,parent.frame())
  weight.type <- "HCG"
  if(length(outcome.variable) == 1){
    result <- RDS.estimates.local(rds.data,outcome.variable,
                                  subset=subset, empir.lik=empir.lik,weight.type=weight.type, N=N, small.fraction=small.fraction,to.factor=to.factor,cont.breaks=cont.breaks)
  }
  else {
    result <- lapply(outcome.variable,function(g){
      RDS.estimates.local(rds.data,g,subset=subset,
                          empir.lik=empir.lik,weight.type=weight.type, N=N, small.fraction=small.fraction,to.factor=to.factor,cont.breaks=cont.breaks)
    })
    names(result) <- outcome.variable
    
  }
  return(result)
}
