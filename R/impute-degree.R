
#' Imputes missing degree values
#' @param rds.data an rds.data.frame
#' @param trait.variable the name of the variable in rds.data to stratify the imputation by
#' @param N population size
#' @param method If mean, the weighted mean value is imputed, otherwize a quantile is used.
#' @param quantile If method is "quantile", this is the quantile that is used. Defaults to median
#' @param recruitment.lower.bound If TRUE, then for each individual, the degree is taken to be the minimum of
#' the number of recruits plus one, and the reported degree
#' @param round.degree Should degrees be integer rounded.  
#' @details 
#' This function imputes degree values using the weighted mean or quantile values of the non-missing degrees.
#' Weights are calcualted using Gile's SS if N is not NULL, or RDS-II if it is. If a trait variable is specified,
#' means and quantile are calculated within the levels of the trait variable
#' @examples 
#' data(faux)
#' rds.data <- faux
#' rds.data$network.size[c(1,2,30,52,81,101,108,111)] <- NA
#' impute.degree(rds.data)
#' impute.degree(rds.data,trait.variable="X")
#' impute.degree(rds.data,trait.variable="X",method="quantile")
#' @export impute.degree
impute.degree <- function(rds.data, trait.variable=NULL, N=NULL, method=c("mean", "quantile"), 
                          quantile=0.50, recruitment.lower.bound=TRUE, round.degree=TRUE){
  assert.valid.rds.data.frame(rds.data)
  
  method <- match.arg(method)
  
  n <- nrow(rds.data)
  if(is.null(trait.variable)){
    by.var <- rep(1,n)
  }else{
    by.var <- as.numeric(as.factor(rds.data[[trait.variable]]))
  }
  nlev <- length(unique(by.var))
  
  nrecruit <- get.number.of.recruits(rds.data)
  degree <- get.net.size(rds.data)
  if(recruitment.lower.bound)
    degree <- ifelse(degree < nrecruit + 1, nrecruit + 1, degree)
  d <- na.omit(degree)
  if(!is.null(N))
    wts <- gile.ss.weights(d, N)
  else
    wts <- vh.weights(d, N)
  nadeg <- is.na(degree)
  imp.value <- rep(NA,nlev)
  for(i in 1:nlev){
    if(method=="mean")
      imp.value[i] <- wtd.mean(d[by.var == i], wts[by.var == i])
    else 
      imp.value[i] <- wtd.quantile(d[by.var == i], wts[by.var == i], probs=quantile)    
  }
  imputed <- imp.value[by.var]
  newdeg <- degree
  newdeg[nadeg] <- imputed[nadeg]
  if(round.degree) 
    newdeg <- round(newdeg)
  newdeg
}



