#' Compute a test of trend in prevalences based on a likelihood-ratio statistic
#' 
#' This function takes a series of point estimates and their associated standard errors and
#' computes the p-value for the test of a monotone decrease in the 
#' population prevalences (in sequence order). 
#' The p-value for a monotone increase is
#' also reported. An optional plot of the estimates and the null distribution of the test statistics is provided.
#' More formally, let the \eqn{K} population prevalences in sequence order be \eqn{p_1, \ldots, p_K}.
#' We test the null hypothesis:\cr
#' \deqn{H_0 : p_1 = \ldots = p_K}
#' vs
#' \deqn{H_1 : p_1 \ge p_2 \ldots \ge p_K}
#' with at least one equality strict. A likelihood ratio statistic for this test has 
#' been derived (Bartholomew 1959).
#' The null distribution of the likelihood ratio statistic is very complex 
#' but can be determined by a simple Monte Carlo.\cr
#' The function requires the isotone library.
#' 
#' @aliases LRT.trend
#' @aliases LRT.trend.test
#' @param data A two row matrix or data.frame of prevalence estimates and
#' their standard errors. The first row is the prevalence estimates and the
#' second are the standard errors. The column are the comparison groups in the
#' order (e.g., time) there are to be assessed. The row names of \code{data}
#' should be "estimate" and "sigma". This is 
#' @param variables A character vector of column names it select from \code{data}.
#' @param confidence.level The confidence level for the confidence intervals. The default is 0.95 for 95\%.
#' @param number.of.bootstrap.samples The number of Monte Carlo draws to
#' determine the null distribution of the likelihood ratio statistic.
#' @param plot A character vector of choices, a subset of \code{estimates}, \code{distributions}.
#' If \code{estimates} is given then a plot of the estimates and nominal 95\% confidence bands (as error bars) is produced.
#' If \code{distributions} is given then a plot is produced of the null distributions of 
#' the likelihood
#' ratio statistic with the observed statistics plotted as a vertical dashed line.
#' @param seed The value of the random number seed. Preset by default to allow reproducibility.
#' @return 
#' A list with components
#' \itemize{ \item\code{pvalue.increasing}: The p-value for the test of a monotone increase in population prevalence.
#' \item\code{pvalue.decreasing}: The p-value for the test of a monotone decrease in population prevalence.
#' \item\code{L}: The value of the likelihood-ratio statistic.
#' \item\code{x}: The passed vector of prevalence estimates in the order (e.g., time).
#' \item\code{sigma} The passed vector of standard error estimates corresponding to \code{x}.
#' }
#' 
#' @author Mark S. Handcock
#' @references Bartholomew, D. J. (1959). A test of homogeneity for ordered alternatives. Biometrika 46 36-48. 
#' @keywords survey manip
#' @examples
#' 
#' d <- t(data.frame(estimate=c(0.16,0.15,0.3), sigma=c(0.04,0.04,0.1)))
#' colnames(d) <- c("time_1","time_2","time_3") 
#' LRT.trend.test(d,number.of.bootstrap.samples=1000)
#' @export 
LRT.trend.test <- function(data,variables=colnames(data),confidence.level=0.95,
   number.of.bootstrap.samples=5000,plot=NULL,seed=1){
   out <- LRT.trend(x=as.numeric(data["estimate",variables]),
                    sigma=as.numeric(data["sigma",variables]),
                    number.of.bootstrap.samples=number.of.bootstrap.samples,
                    confidence.level=confidence.level,
                    plot=plot,seed=seed)
   invisible(out)
}
# Distribution
#' @export 
LRT.trend <- function(x,sigma,number.of.bootstrap.samples=5000,confidence.level=0.95,plot=NULL,seed=1,smooth=0.35){
# require(isotone, quietly=TRUE, warn.conflicts=FALSE)
# if(plot){
#   suppressMessages(require(locfit, quietly = TRUE))
# }
  a <- 1/sigma^2
  xbar <- sum(a*x)/sum(a)
  k <- length(x)
  set.seed(seed)
  r <- matrix(stats::rnorm(k*number.of.bootstrap.samples),ncol=k)
  r <- sweep(r,2,sigma,"*") + xbar
  L <- t(apply(r,1,LRT.value.trend,sigma=sigma))
  obsL <- LRT.value.trend(x,sigma)
  options("scipen"=100)
  cat(sprintf("In a test of the null hypothesis of equal proportions against the alternative that at least two are unequal, the p-value is %s.\n",format.pval(mean(obsL[3] <= L[,3]),eps=0.0001)))
  if(mean(obsL[3] <= L[,3]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of equal proportions in favor of the alternative that at least two are unequal (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of equal proportions (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
  cat(sprintf("In a test of the null hypothesis of equal proportions against the alternative of a monotone trend (either up or down), the p-value is %s.\n",format.pval(mean(obsL[4] <= L[,4]),eps=0.0001)))
  if(mean(obsL[4] <= L[,4]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of equal proportions in favor of the alternative of a monotone trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of equal proportions (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
  cat(sprintf("In a test of the null hypothesis of equal proportions against the alternative of an increasing trend, the p-value is %s.\n",format.pval(mean(obsL[1] <= L[,1]),eps=0.0001)))
  if(mean(obsL[1] <= L[,1]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of equal proportions in favor of the alternative of an increasing trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of equal proportions (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
  cat(sprintf("In a test of the null hypothesis of equal proportions against the alternative of an decreasing trend, the p-value is %s.\n",format.pval(mean(obsL[2] <= L[,2]),eps=0.0001)))
  if(mean(obsL[2] <= L[,2]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of equal proportions in favor of the alternative of a decreasing trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of equal proportions (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
#
  if("estimates" %in% plot){
    yminus <- x - 1.96*sigma
    yplus  <- x + 1.96*sigma
    errorbar(x=seq_along(x),y=x,yminus=yminus,
     yplus=yplus,xlab="sequence",ylab="estimate",
     main="Trend of estimates")
  }
  if("distributions" %in% plot){
    binn <- 500
    maxl <- max(1.1*obsL[3],stats::quantile(L[,3],0.95))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,3],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,3],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,3]),xlim=c(0,max(obsL[3],stats::quantile(L[,3],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Overall Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[3],lty=2)
    maxl <- max(1.1*obsL[1],stats::quantile(L[,1],0.95))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,1],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,1],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,1]),xlim=c(0,max(obsL[1],stats::quantile(L[,1],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Increasing Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[1],lty=2)
    maxl <- max(1.1*obsL[2],stats::quantile(L[,2],0.95))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,2],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,2],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,2]),xlim=c(0,max(obsL[2],stats::quantile(L[,2],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Decreasing Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[2],lty=2)
    maxl <- max(1.1*obsL[4],stats::quantile(L[,4],0.95))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,4],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,4],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,1]),xlim=c(0,max(obsL[1],stats::quantile(L[,1],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Monotone Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[4],lty=2)
  }
  invisible(list(pvalue.increasing=mean(obsL[1] <= L[,1]),
       pvalue.decreasing=mean(obsL[2] <= L[,2]),
       pvalue.overall=mean(obsL[3] <= L[,3]),
       L=L,x=x,sigma=sigma))
}
#' Compute a test of trend in prevalences based on a likelihood-ratio statistic
#' 
#' This function takes a series of point estimates and their associated standard errors and
#' computes the p-value for the test of a monotone decrease in the 
#' population prevalences (in sequence order). 
#' The p-value for a monotone increase is
#' also reported.
#' More formally, let the \eqn{K} population prevalences in sequence order be \eqn{p_1, \ldots, p_K}.
#' We test the null hypothesis:\cr
#' \deqn{H_0 : p_1 = \ldots = p_K}
#' vs
#' \deqn{H_1 : p_1 \ge p_2 \ldots \ge p_K}
#' with at least one equality strict. A likelihood ratio statistic for this test has 
#' been derived (Bartholomew 1959).
#' The null distribution of the likelihood ratio statistic is very complex 
#' but can be determined by a simple Monte Carlo.\cr
#' The function requires the isotone library.
#' 
#' @aliases LRT.value.trend
#' @param x A vector of prevalence estimates in the order (e.g., time).
#' @param sigma A vector of standard error estimates corresponding to \code{x}.
#' @return 
#' A list with components
#' \itemize{ \item\code{pvalue.increasing}: The p-value for the test of a monotone increase in population prevalence.
#' \item\code{pvalue.decreasing}: The p-value for the test of a monotone decrease in population prevalence.
#' \item\code{L}: The value of the likelihood-ratio statistic.
#' \item\code{x}: The passed vector of prevalence estimates in the order (e.g., time).
#' \item\code{sigma} The passed vector of standard error estimates corresponding to \code{x}.
#' }
#' 
#' @author Mark S. Handcock
#' @references Bartholomew, D. J. (1959). A test of homogeneity for ordered alternatives. Biometrika 46 36-48. 
#' @keywords survey manip
#' @examples
#' 
#' \dontrun{
#' x <- c(0.16,0.15,0.3)
#' sigma <- c(0.04,0.04,0.1)
#' LRT.value.trend(x,sigma)
#' }
#' @export 
LRT.value.trend <- function(x,sigma){
  a <- 1/sigma^2
# xbar <- mean(x)
  xbar <- sum(a*x)/sum(a)
  s <- seq_along(x)
  Atot <- cbind(s[-length(s)], s[-1])
  fit.ls1 <- isotone::activeSet(Atot, "LS", y = x, weights = a)
# m=fit.ls1$x
# LRT.increasing <- sum(a*(x-xbar)^2)-sum(a*(x-m)^2)
  LRT.increasing <- sum(a*(x-xbar)^2)-fit.ls1$fval
  Atot <- cbind(s[-1],s[-length(s)])
  fit.ls2 <- isotone::activeSet(Atot, "LS", y = x, weights = a)
  #m=fit.ls2$x
  #LRT.decreasing <- sum(a*(x-xbar)^2)-sum(a*(x-m)^2)
  LRT.decreasing <- sum(a*(x-xbar)^2)-fit.ls2$fval
  LRT.overall <- sum(a*(x-xbar)^2)
  LRT.monotone <- sum(a*(x-xbar)^2)-min(fit.ls1$fval,fit.ls2$fval,na.rm=TRUE)
  H0.monotone <- min(fit.ls1$fval,fit.ls2$fval,na.rm=TRUE)
  H0.increasing <- fit.ls1$fval
  H0.decreasing <- fit.ls2$fval
  return(c(LRT.increasing,LRT.decreasing,LRT.overall,LRT.monotone,H0.increasing,H0.decreasing,H0.monotone))
}
#
LRT.trend.null <- function(x,sigma,number.of.bootstrap.samples=5000,confidence.level=0.95,plot=NULL,seed=1,smooth=0.35){
# require(isotone, quietly=TRUE, warn.conflicts=FALSE)
# if(plot){
#   suppressMessages(require(locfit, quietly = TRUE))
# }
  a <- 1/sigma^2
  xbar <- sum(a*x)/sum(a)
  k <- length(x)
  set.seed(seed)
  s <- seq_along(x)
  Atot <- cbind(s[-length(s)], s[-1])
  fit.ls1 <- isotone::activeSet(Atot, "LS", y = x, weights = a)
  m.inc=fit.ls1$x
  r <- matrix(stats::rnorm(k*number.of.bootstrap.samples),ncol=k)
  r <- sweep(r,2,sigma,"*") + m.inc
  L.inc <- t(apply(r,1,LRT.value.trend,sigma=sigma))
  Atot <- cbind(s[-1],s[-length(s)])
  fit.ls2 <- isotone::activeSet(Atot, "LS", y = x, weights = a)
  m.dec=fit.ls2$x
  r <- matrix(stats::rnorm(k*number.of.bootstrap.samples),ncol=k)
  r <- sweep(r,2,sigma,"*") + m.dec
  L.dec <- t(apply(r,1,LRT.value.trend,sigma=sigma))
  L <- L.inc
  L[,5] <- L.dec[,5]
  L[,6] <- pmin(L.inc[,6],L.dec[,6])
  obsL <- LRT.value.trend(x,sigma)
#
  options("scipen"=100)
  cat(sprintf("In a test of the null hypothesis of an increasing trend in proportions against the complementary hypothesis, the p-value is %s.\n",format.pval(mean(obsL[4] <= L[,4]),eps=0.0001)))
  if(mean(obsL[4] <= L[,4]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of an increasing trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of an increasing trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
  cat(sprintf("In a test of the null hypothesis of an decreasing trend in proportions against the complementary hypothesis, the p-value is %s.\n",format.pval(mean(obsL[5] <= L[,5]),eps=0.0001)))
  if(mean(obsL[5] <= L[,5]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of an decreasing trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of an decreasing trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
  cat(sprintf("In a test of the null hypothesis of a monotone trend in proportions against the complementary hypothesis, the p-value is %s.\n",format.pval(mean(obsL[6] <= L[,6]),eps=0.0001)))
  if(mean(obsL[6] <= L[,6]) < (1-confidence.level)){
    cat(sprintf("We reject the null hypothesis of a monotone trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }else{
    cat(sprintf("We do not reject the null hypothesis of a monotone trend (at the %s%% level).\n\n",format(100*(1-confidence.level))))
  }
#
#
#
  if("estimates" %in% plot){
    yminus <- x - 1.96*sigma
    yplus  <- x + 1.96*sigma
    errorbar(x=seq_along(x),y=x,yminus=yminus,
     yplus=yplus,xlab="sequence",ylab="estimate",
     main="Trend of estimates")
  }
  if("distributions" %in% plot){
    binn <- 500
    maxl <- max(1.1*obsL[6],stats::quantile(L[,3],0.99))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,6],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,6],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,6]),xlim=c(0,max(obsL[6],stats::quantile(L[,6],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Overall Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[6],lty=2)
    maxl <- max(1.1*obsL[4],stats::quantile(L[,4],0.99))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,4],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,4],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,4]),xlim=c(0,max(obsL[4],stats::quantile(L[,4],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Increasing Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[4],lty=2)
    maxl <- max(1.1*obsL[5],stats::quantile(L[,5],0.99))
    r <- seq(0, maxl, length = binn + 1)[-1] - 0.5*maxl/binn
#   yl <- locfit::locfit.raw(locfit::lp(L[,5],nn=0.1,h=0.8), xlim=c(0,maxl))
#   gpdf <- predict(yl, newdata=r)
    gpdf=bgk_kde(L[,5],n=2^(ceiling(log(maxl)/log(2))),MIN=0)
    gpdf <- stats::spline(x=gpdf[1,],y=gpdf[2,],xout=r)$y
    scalef <- binn/sum(gpdf)
    gpdf <- gpdf * scalef
#   maxl <- r[which.max(cumsum(gpdf)>binn*0.99)]
#  plot(density(L[,5]),xlim=c(0,max(obsL[5],stats::quantile(L[,5],0.99))),
   plot(x=r,y=gpdf,xlim=c(0,maxl), type="l",
     ylab="Density",
     xlab="Likelihood Ratio Statistic",main="Null Distribution of the Decreasing Test Statistic",sub="The vertical line is the observed statistic")
   abline(v=obsL[5],lty=2)
  }
  invisible(list(pvalue.increasing=mean(obsL[1] <= L[,1]),
       pvalue.decreasing=mean(obsL[2] <= L[,2]),
       pvalue.overall=mean(obsL[3] <= L[,3]),
       L=L,x=x,sigma=sigma))
}
