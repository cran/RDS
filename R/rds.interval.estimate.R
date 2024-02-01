#' An object of class rds.interval.estimate
#' 
#' This function creates an object of class \code{rds.interval.estimate}.
#' 
#' 
#' @param estimate The numerical point estimate of proportion of the
#' \code{trait.variable}.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical variable to be analyzed.
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I
#' (DS)"}, and \code{"Arithemic Mean"}. If \code{NULL} it defaults to
#' \code{"Gile's SS"}.
#' @param uncertainty A string giving the type of uncertainty estimator to use.
#' The options are \code{"SRS"}, \code{"Gile"} and \code{"Salganik"}. This is
#' usually determined by \code{weight.type} to be consistent with the
#' estimator's origins. The estimators \code{"RDS-I"}, \code{"RDS-I (DS)"}, \code{"RDS-II"} default to
#' \code{"Salganik"}, "Arithmetic Mean" defaults to \code{"SRS"} and "Gile's
#' SS" defaults to the \code{"Gile"} bootstrap.
#' @param weights A numerical vector of sampling weights for the sample, in
#' order of the sample.  They should be inversely proportional to the
#' first-order inclusion probabilites, although this is not assessed or
#' inforced.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{pop.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing it defaults to 1000.
#' @param conf.level The confidence level for the confidence intervals. The
#' default is 0.95 for 95\%.
#' @param csubset A character string representing text to add to the output label. Typically
#' this will be the expression used it define the subset of the data used  for the estimate.
#' @return An object of class \code{rds.interval.estimate} is returned. This is
#' a list with components \itemize{ \item\code{estimate}: The numerical point
#' estimate of proportion of the \code{trait.variable}.  \item\code{interval}:
#' A matrix with six columns and one row per category of \code{trait.variable}:
#' \itemize{ \item\code{point estimate}: The HT estimate of the population
#' mean.  \item\code{95\% Lower Bound}: Lower 95\% confidence bound.
#' \item\code{95\% Upper Bound}: Upper 95\% confidence bound.  \item\code{Design
#' Effect}: The design effect of the RDS.  \item\code{s.e.}: Standard error.
#' \item\code{n}: Count of the number of sample values with that value of the
#' trait.  } }
#' @export
#' @author Mark S. Handcock
#' @references Gile, Krista J., Handcock, Mark S., 2010. Respondent-driven Sampling:
#' An Assessment of Current Methodology, Sociological Methodology, 40,
#' 285-327. <doi:10.1111/j.1467-9531.2010.01223.x>
#'
#' Gile, Krista J., Beaudry, Isabelle S. and Handcock, Mark S., 2018 
#' Methods for Inference from Respondent-Driven Sampling Data,
#' Annual Review of Statistics and Its Application
#' <doi:10.1146/annurev-statistics-031017-100704>.
#' 
#' Salganik, M., Heckathorn, D. D., 2004. \emph{Sampling and estimation in
#' hidden populations using respondent-driven sampling}. Sociological
#' Methodology 34, 193-239.
#' 
#' Volz, E., Heckathorn, D., 2008. \emph{Probability based estimation theory
#' for Respondent Driven Sampling}. The Journal of Official Statistics 24 (1),
#' 79-97.
#' @keywords manip survey
#' @examples
#' 
#' data(faux)
#' RDS.I.estimates(rds.data=faux,outcome.variable='X',smoothed=TRUE)
#' 

rds.interval.estimate <- function(estimate, outcome.variable, 
                                  weight.type, uncertainty, weights, N=NULL,conf.level=.95,csubset="") {
  
  .Object <- list()
  
  .Object$estimate <- estimate[1:(length(estimate)/6)]
  .Object$interval <- estimate
  .Object$outcome.variable <- outcome.variable
  .Object$weight.type <- weight.type
  .Object$uncertainty <- uncertainty
  .Object$weights <- weights
  .Object$N <- N
  .Object$conf.level <- conf.level
  .Object$csubset <- csubset
  class(.Object) <- "rds.interval.estimate"
  return(.Object)
}

#' Prints an \code{rds.interval.estimate} object
#' @param x an \code{rds.interval.estimate} object
#' @param as.percentage logical. Print the interval estimates
#' as percentages (as distinct from proportions). 
#' The default, NULL, means that it will determine if
#' the variable is discrete or continuous and only print them as percentages if they are discrete.
#' @param ... unused
#' @export
#' @method print rds.interval.estimate
print.rds.interval.estimate <- function(x, as.percentage=NULL, ...) {
  fmt <- function(x,...){
    format(x,...,scientific=FALSE)
  }
  matest <- matrix(x$interval, ncol = 6, byrow = FALSE)
  if (nrow(matest) > 1) {
    rownames(matest) <- names(x$interval)[1:nrow(matest)]
  }
  else {
    rownames(matest) <- x$outcome.variable
    names(x$interval) <- rep(x$outcome.variable, length(x$interval))
  }
  rownames(matest)[is.na(rownames(matest))] <- "NA"
  mnames <- max(nchar(names(x$interval)[1:nrow(matest)]))
  colnames(matest) <- c("point", "lower", "upper", "Design Effect", 
                        "s.e.", "n")
  
  if(is.null(as.percentage)){
    as.percentage <- attr(x, "is.cts")
    as.percentage <- if(is.logical(as.percentage) && (as.percentage==FALSE)){
           as.percentage <- TRUE
     }else{as.percentage <- FALSE}
  }

  if(as.percentage){
    matest[,c(1,2,3,5)] <- 100*matest[,c(1,2,3,5)]
  }

  nsamples <- sum(matest[,ncol(matest)])
  fmatest <- rbind(matest,"")
  fmatest[nrow(fmatest),ncol(fmatest)-1] <- "Total"
  if(as.percentage){
    fmatest[,ncol(fmatest)] <- fmt(c(matest[,ncol(matest)],nsamples), width = 5, digits = 5)
    fmatest[-nrow(fmatest),  1] <- fmt(matest[,  1], width = 8, digits=3)
    fmatest[-nrow(fmatest),2:4] <- fmt(matest[,2:4], width = 8, digits=3)
    fmatest[-nrow(fmatest),  4] <- fmt(matest[,  4], width = 6, digits=3)
    fmatest[-nrow(fmatest),  5] <- fmt(matest[,  5], width = 8, digits=3)
  }else{
    fmatest[,ncol(fmatest)] <- fmt(c(matest[,ncol(matest)],nsamples), width = 5, digits = 5)
    fmatest[-nrow(fmatest),  1] <- fmt(matest[,  1], width = 8, digits=4)
    fmatest[-nrow(fmatest),2:4] <- fmt(matest[,2:4], width = 8, digits=4)
    fmatest[-nrow(fmatest),  4] <- fmt(matest[,  4], width = 6, digits=3)
    fmatest[-nrow(fmatest),  5] <- fmt(matest[,  5], width = 8, digits=3)
  }
  
  clp <- x$conf.level*100
  tmp <- as.data.frame(fmatest,stringsAsFactors=FALSE)
  tmp[-nrow(tmp),2] <- paste("(",tmp[,2],", ",tmp[,3],")",sep="")[-nrow(tmp)]
  tmp <- tmp[,-3]
  colnames(tmp) <- c("Estimate", paste0(clp,"% Interval"), "Design Effect", "Std. Error", "N")
  if (is.element("DeducerRichOutput", .packages())) {
    colnames(fmatest) <- c("Point Estimate", paste0(" ",clp,"% Lower\n Bound"), 
                           paste0(" ",clp,"% Upper\n Bound"), "Estimated\n Design Effect", 
                           "Standard Error", "Sample Size")
    if(as.percentage){
      colnames(fmatest)[1] <- "Point Estimate (%)"
    }
    
    get("print_to_html")(fmatest, caption.placement = "top", digits = c(8, 
                                                                        5, 5, 5, 1, 5, 0), caption = paste(x$weight.type, 
                                                                                                           "Estimate for", x$outcome.variable,
                                                                                                           switch(((x$csubset=="")|(x$csubset=="NULL"))+1,paste("[",x$csubset,"]",sep=""),NULL)))
  }
  else {
    if(as.percentage){
      colnames(tmp)[1] <- "Estimate (%)"
    }
    cat(paste(c(x$weight.type, "Estimate for", 
                x$outcome.variable,
                switch(((x$csubset=="")|(x$csubset=="NULL"))+1,paste("[",x$csubset,"]",sep=""),NULL),"\n")))
    print(tmp)
    if(!is.null(x$N))
      cat("* Using population size estimate:",x$N,"\n")
  }
  return(invisible(tmp))
}

#' Is an instance of rds.interval.estimate
#' @param x An object to be tested.
#' @export
is.rds.interval.estimate <- function(x) inherits(x,"rds.interval.estimate")

#' Is an instance of rds.interval.estimate.list
#' This is a (typically time ordered) sequence of RDS estimates of a comparable quantity
#' @param x An object to be tested.
#' @export
is.rds.interval.estimate.list <- function(x) inherits(x,"rds.interval.estimate.list")

#' Convert the output of print.rds.interval.estimate from a character data.frame to a numeric matrix
#' @param x An object, typically the result of print.rds.interval.estimate.
#' @param proportion logical, Should the outcome be treated as a proportion and converted to a percentage.
#' @export
export.rds.interval.estimate <- function(x,proportion=TRUE){
  a <- as.matrix(x)
  a <- a[-nrow(a),,drop=FALSE]
  b <- matrix(0,ncol=6,nrow=nrow(a))
  for(i in 1:nrow(a)){
    b[i,2:3] <- eval(parse(text=paste("c",a[i,2],sep="")))
  }
  b[,c(1,4:6)] <- as.numeric(a[,c(1,3:5)])
  colnames(b) <- c(colnames(a)[1],"95% lower","95% upper",colnames(a)[3:ncol(a)])
  rownames(b) <- rownames(a)
  if(proportion){
    b[,c(1,2,3,5)] <- 100*b[,c(1,2,3,5)]
  }
  as.data.frame(b)
}
