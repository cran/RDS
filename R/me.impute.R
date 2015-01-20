#
# Calculate the Measurement error model MLE
#
#
#' Estimates each person's personal network size (degree) based on their self-reported degree and the 
#'  number of their (direct) recruits. It uses the time the person was recruited as a factor in 
#'  determining the number of recruits they produce.
#' @param rds.data An rds.data.frame
#' @param max.coupons The number of recruitment coupons distributed to each 
#' 		enrolled subject (i.e. the maximum number of recruitees for any subject).
#'              By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param type.impute The type of imputation based on the conditional distribution. 
#'     It can be of type \code{distribution},\code{mode},\code{median}, or \code{mean} 
#'     with the first , the default, being a random draw from the conditional distribution.
#' @param recruit.time vector; A numerical value for the data/time that the person was interviewed.
#' @param reflect.time logical; If \code{FALSE} then the \code{recruit.time} is the time before the 
#' end of the study (instead of the time since the survey started or chronological time).
#' @param unit.scale numeric; If not \code{NULL} it sets the numeric value of the scale parameter of the
#' distribution of the unit sizes.
#' For the negative bionimial, it is the multiplier on the variance of the negative binomial 
#' compared to a Poisson (via the Poisson-Gamma mixture representation). Sometimes the scale is 
#' unnaturally large (e.g. 40) so this give the option of fixing it (rather than use the MLE of it). The
#' model is fit with the parameter fixed at this passed value.
#' @param unit.model The type of distribution for the unit sizes.
#'  It can be of \code{nbinom}, meaning a negative binomial. In this case \code{unit.scale} is the multiplier 
#' on the variance of the negative binomial compared to a Poisson of the same mean.
#' The alernative is \code{cmp}, meaning a Conway-Maxwell-Poisson distribution. In this case \code{ unit.scale}
#' is the scale parameter compared to a Poisson of the same mean (values less than one mean 
#' under-dispersed and values over one mean over-dispersed).
#' @param guess vector; if not \code{NULL}, the initial parameter values for the MLE fitting. 
#' @param verbose logical; if this is \code{TRUE}, the program will print out additional
#   information about the fitting process.
#' @export
#' @examples
#' \dontrun{
#' data(fauxmadrona)
#' # The next line fits the model for the self-reported personal
#' # network sizes and imputes the personal network sizes 
#' # It may take up to 60 seconds.
#' idegree <- impute.degree(fauxmadrona)
#' # frequency of estimated personal network sizes
#' table(idegree)
#' }
impute.degree <-function(rds.data,max.coupons=NULL,
	type.impute = c("distribution","mode","median","mean"),
        recruit.time=NULL,reflect.time=TRUE, unit.scale=NULL, 
	unit.model = c("nbinom","cmp"),
	guess=NULL,
        verbose=FALSE){
 	if(!is(rds.data,"rds.data.frame"))
 		stop("rds.data must be of type rds.data.frame")   
	
	if(missing(unit.model)){
		unit.model <- "nbinom"
	}
	unit.model <- match.arg(unit.model, c("nbinom","cmp"))
	n <- nrow(rds.data)
 	if(is.null(attr(rds.data,"network.size.variable")))
 		stop("rds.data must have a network.size attribute.")
        nr <- get.number.of.recruits(rds.data)
 	if(is.null(max.coupons)){
 	  max.coupons <- attr(rds.data,"max.coupons")
 	  if(is.null(max.coupons)){
 		max.coupons <- max(nr,na.rm=TRUE)
          }
        }
        if(is.null(recruit.time)){
	 recruit.time <- 1:n
        }
        if(is.character(recruit.time)){
	 recruit.time <- as.numeric(rds.data[[recruit.time]])
        }
        if(reflect.time){
	 recruit.time <- max(recruit.time)-recruit.time
        }
	network.size <- as.numeric(rds.data[[attr(rds.data,"network.size.variable")]])
	remvalues <- is.na(network.size)
	if(any(remvalues)){
          warning(paste(sum(remvalues),"of",nrow(rds.data),
                   "network sizes were missing. These will be imputed from the marginal distribution"), call. = FALSE)
        }
	if(missing(type.impute)){type.impute <- "distribution"}
        type.impute <- match.arg(type.impute,
	 c("distribution","mode","median","mean"))
        if(is.na(type.impute)) { # User typed an unrecognizable name
          stop(paste('You must specify a valid type.impute. The valid types are "distribution","mode","median", and "mean"'), call.=FALSE)
        }

        nw <- get.wave(rds.data)
        ns <- get.seed.id(rds.data)
        is.seed <- (get.rid(rds.data)=="seed")
        nsize <- pmax(network.size,nr+!is.seed)
        if(is.null(guess)){
         if(is.null(unit.scale)){
          if(unit.model=="cmp"){
            guess <- c(38,-5,0,0.14,1.5,2)
          }else{
            guess <- c(38,-5,0,0.14,1.5,2)
          }
         }else{
          if(unit.model=="cmp"){
            guess <- c(38,-5,0,0.14,2)
          }else{
            guess <- c(38,-5,0,0.14,2)
          }
         }
        }
        fit <- memle(guess=guess,network.size=nsize[!remvalues],num.recruits=nr[!remvalues],recruit.time=recruit.time[!remvalues],max.coupons=max.coupons,unit.scale=unit.scale,unit.model=unit.model)
        if(verbose){
         print(summary(fit))
        }
        a=dmepdf(fit$coef,nsize,nr,recruit.time,unit.scale=unit.scale,unit.model=unit.model)

	is <- switch(type.impute, 
		`distribution` = {
			ff <- function(x){sample.int(n=length(x),size=1,prob=x)};
			apply(a,2,ff)
				 },
		`mode` = {
			apply(a,2,which.max)
				 },
		`median` = {
			ff <- function(x){seq_along(x)[match(TRUE,cumsum(x) >= 0.5)]};
			apply(a,2,ff)
				 },
		`mean` = {
			apply(a,2,function(x){sum((1:length(x))*x)})
				 }
			 )

	return(is)
}

"is.psd" <- function(V, tol = 1e-12){
  if(is.null(V)){return(FALSE)}
  if(sum(is.na(V))>0){
   ev <- FALSE
  }else{
   ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
   ev <- all(ev/max(ev) > tol) & all(ev >=  - tol * abs(ev[1]))
   if(is.na(ev)){ev <- FALSE}
   if(ev != TRUE){ev <- FALSE}
  }
  ev
}
#
# Complete data log-likelihoods
#
llmeall <- function(v,x,network.size,num.recruits,recruit.time,max.coupons=3,cutoff=0,cutabove=1000,np=6,unit.model="nbinom"){
 llme <- switch(unit.model,
  "cmp"=llcmpme, "nbinom"=llnbme)
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llme(v=v,network.size=network.size,num.recruits=num.recruits,recruit.time=recruit.time,max.coupons=max.coupons,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
dmepdfR <- function(v,network.size,num.recruits,recruit.time,max.coupons=3,K=100,nb.scale=NULL){
 g <- 1:K
 mean = v[1]
 beta0 = v[2]
 beta1 = v[3]
 ln.sd=v[4]
 if(is.null(nb.scale)){
   scale=v[5]
   opt=v[6]
 }else{
   scale=nb.scale
   opt=v[5]
 }
#
 n <- length(network.size)
 g <- 1:K
 rt <- sort(unique(recruit.time))
 rtprob <- beta0 + beta1*rt
 rtprob <- exp(rtprob)/(1+exp(rtprob))
#
 rx <- sort(unique(num.recruits))
 tr <- sort(unique(network.size))
 m.recruit.time <- match(recruit.time,rt)
 m.num.recruits <- match(num.recruits,rx)
 m.network.size <- match(network.size,tr)
 pdf <- dnbinom(x=g,mu=mean,size=scale)
 b1=array(0,dim=c(length(rx),length(g),length(rt)))
 for(x1 in seq_along(rx)){for(x2 in seq_along(g)){for(x3 in seq_along(rt)){
   if((g[x2]<=max.coupons)|(rx[x1]<max.coupons)){
    b1[x1,x2,x3] <- dbinom(x=rx[x1],size=g[x2],prob=rtprob[x3])
   }else{
    b1[x1,x2,x3] <- 1-pbinom(q=max.coupons-1,size=g[x2],prob=rtprob[x3])
   }
 }}}
 b2=outer(log(tr),log(g),function(x1,x2){dnorm(x=x1-x2+log(opt),sd=ln.sd,log=FALSE)})
 b=matrix(0,ncol=length(num.recruits),nrow=length(g))
 for(x1 in seq_along(network.size)){
  if(is.na(m.network.size[x1])){
   a = b1[m.num.recruits[x1],,m.recruit.time[x1]]*pdf
  }else{
   a = (b2[m.network.size[x1],]*b1[m.num.recruits[x1],,m.recruit.time[x1]])*pdf
  }
  b[,x1] = a/sum(a)
 }
 b
}
#  File R/summary.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
###############################################################################
# The <summary.ergm> function prints a 'summary of model fit' table and returns
# the components of this table and several others listed below
#
# --PARAMETERS--
#   object     : an ergm object
#
#
# --IGNORED PARAMETERS--
#   ...        : used for flexibility
#   digits     : significant digits for the coefficients;default=
#                max(3,getOption("digits")-3), but the hard-coded value is 5
#   correlation: whether the correlation matrix of the estimated parameters
#                should be printed (T or F); default=FALSE
#   covariance : whether the covariance matrix of the estimated parameters
#                should be printed (T or F); default=FALSE
#   eps        : the numerical tolerance inputted to the native R function
#                <format.pval>; the code using 'eps' is now commented out;
#                default=.0001
#
# --RETURNED--
#   ans: a "summary.ergm" object as a list containing the following
#      formula         : object$formula
#      randomeffects   : object$re
#      digits          : the 'digits' inputted to <summary.ergm> or the default
#                        value (despite the fact the digits will be 5)
#      correlation     : the 'correlation' passed to <summary.ergm>
#      degeneracy.value: object$degenarcy.value
#      covariance      : the 'covariance' passed to <summary.ergm>
#      iterations      : object$iterations
#      samplesize      : NA if 'pseudolikelihood'=TRUE, object$samplesize otherwise
#      message         : a message regarding the validity of the standard error
#                        estimates
#      aic             : the AIC goodness of fit measure
#      bic             : the BIC goodness of fit measure
#      coefs           : the dataframe of parameter coefficients and their
#                        standard erros and p-values
#      asycov          : the asymptotic covariance matrix
#      asyse           : the asymptotic standard error matrix
#
################################################################################

summary.me <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE,
                          eps=0.0001)
{
  control <- object$control

  if(is.null(object$hessian) && is.null(object$covar)){
    object$covar <- diag(NA, nrow=length(object$coef))
  }
  
  if(is.null(object$covar)){
    asycov <- try(robust.inverse(-object$hessian), silent=TRUE)
    if(inherits(asycov,"try-error")){
      asycov <- diag(1/diag(-object$hessian))
    }
  }else{
    asycov <- object$covar
  }
  colnames(asycov) <- rownames(asycov) <- names(object$coef)
  
  asyse <- diag(asycov)
  asyse[asyse<0|is.infinite(object$coef)] <- NA
  asyse <- sqrt(asyse)
  asyse <- matrix(asyse, ncol=length(asyse))
  colnames(asyse) <- colnames(asycov)

  ans <- list(formula=object$formula,
              digits=digits, correlation=correlation,
              covariance=covariance,
              iterations=object$iterations,
              control=object$control)
  
# nodes<- network.size(object$network)
# dyads<- network.dyadcount(object$network,FALSE)-network.edgecount(NVL(get.miss.dyads(object$constrained, object$constrained.obs),network.initialize(1)))
  dyads<- object$df
  df <- length(object$coef)

  rdf <- dyads - df
  tval <- object$coef / asyse
  pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)

  count <- 1
  templist <- NULL
  while (count <= length(names(object$coef)))
    {
     templist <- append(templist,c(object$coef[count],
          asyse[count],pval[count]))
     count <- count+1
    }

  tempmatrix <- matrix(templist, ncol=3,byrow=TRUE)
  colnames(tempmatrix) <- c("Estimate", "Std. Error", "p-value")
  rownames(tempmatrix) <- names(object$coef)

  devtext <- "Deviance:"
  mle.lik<-object$loglik
  null.lik<-object$loglik.null

  ans$null.lik.0 <- is.na(null.lik)

  ans$devtable <- c("",apply(cbind(paste(format(c("    Null", "Residual"), width = 8), devtext), 
                                   format(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik), digits = digits), " on",
                                   format(c(dyads, rdf), digits = digits)," degrees of freedom\n"), 
                             1, paste, collapse = " "),"\n")
    
  ans$aic <- -2*mle.lik +2*(dyads) #AIC(mle.lik)
  ans$bic <- -2*mle.lik +log(df)*(dyads) #BIC(mle.lik)
  
  ans$coefs <- as.data.frame(tempmatrix)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.me"
  ans
}

#  File R/print.summary.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
###############################################################################
# The <print.summary.ergm> function prints a subset of the information given
# by <summary.ergm>
#
# --PARAMETERS--
#   x           : a "summary.ergm" object, as returned by <summary.ergm>
#   digits      : the number of significant digits for the coefficients;
#                 default=max(3, getOption("digits")-3)
#   correlation : whether the correlation matrix of the estimated parameters
#                 should be printed (T or F); default=FALSE
#   covariance  : whether the covariance matrix of the estimated parameters
#                 should be printed (T or F); default=FALSE
#   signif.stars: whether stars are to be printed on summary tables of
#                 coefficients (T or F); default=getOption("show.signif.stars")
#   eps.Pvalue  : the tolerance to be passed to the R <printCoefmat> function;
#                 default=.0001
#   ...         : additional parameters to be passed to <printCoefmat> 
#
# --RETURNED--
#   x
#
###############################################################################

print.summary.me <- function (x, 
              digits = max(3, getOption("digits") - 3),
              correlation=FALSE, covariance=FALSE,
              signif.stars= getOption("show.signif.stars"),
              eps.Pvalue=0.0001, print.header=TRUE, print.formula=TRUE, print.fitinfo=TRUE, print.coefmat=TRUE, print.message=TRUE, print.deviances=TRUE, print.drop=TRUE, ...){
  if(missing(digits)) digits <- x$digits
  
  control <- x$control
  if(print.header){
    cat("\n==========================\n")
    cat("Summary of model fit\n")
    cat("==========================\n\n")
  }
  
  if(print.formula){
    cat("Formula:   ")
    print(x$formula)
    cat("\n")
  }

  if(print.fitinfo){
    if (!is.null(x$iterations)) {
      cat("Iterations: ", x$iterations, "\n")
    }
  }

  if(print.coefmat){
    printCoefmat(x$coefs, digits=digits, signif.stars=signif.stars,
                 P.values=TRUE, has.Pvalue=TRUE, na.print="NA",
                 eps.Pvalue=eps.Pvalue, ...)
  }

  if(print.message){
    if(!is.null(x$message)){ 
      cat(x$message)
    }
    cat("\n")
  }

  if(print.deviances){
    if(!is.null(x$devtable)){
      cat(x$devtable)

      if(x$null.lik.0) cat("Note that the null model likelihood and deviance are defined to be 0.\n\n")
      
      cat(paste("AIC:", format(x$aic, digits = digits), "  ", 
                "BIC:", format(x$bic, digits = digits), "  ",
                "(Smaller is better.)", "\n", sep=" "))
    } 
  }

  if((missing(covariance)&x$covariance)|covariance == TRUE){
    cat("Asymptotic covariance matrix:\n")
    print(x$asycov)
  }
  
  if((missing(correlation)&x$correlation)|correlation == TRUE){
    cat("\nAsymptotic correlation matrix:\n")
    asycor <- x$asycov / crossprod(x$asyse)
    dimnames(asycor) <- dimnames(x$asycov)
    print(asycor)
  }
  
  invisible(x)
}
