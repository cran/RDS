#' Get Horvitz-Thompson estimator assuming inclusion probability proportional
#' to the inverse of network.var (i.e. degree).
#' @param rds.data An rds.data.from
#' @param group.variable The grouping variable.
#' @param network.var The network.size variable.
#' @export
get.h.hat <- function(rds.data,
                      group.variable,
                      network.var=attr(rds.data,"network.size")){
  
  if(is(rds.data,"rds.data.frame")){
    stopifnot(group.variable %in% names(rds.data))
    stopifnot(network.var %in% names(rds.data))
  }
  else{
    stop("rds.data must be of type rds.data.frame")
  }
  
  harmonic.mean <- function(x){
    x <- x[!is.na(x)]
    length(x)/sum(1/x) 	
  }
  
  
  degrees <- rds.data[[network.var]]
  temp <- function(x){  
    deg <- degrees[as.character(rds.data[[group.variable]]) == x]
    return(harmonic.mean(deg))
  }	
  
  group.names <- sort(unique(as.character(rds.data[[group.variable]])))
  group.names <- group.names[!is.na(group.names)]
  
  h.hat <- sapply(group.names,temp)
  names(h.hat) <- group.names
  return(h.hat)    	
}


# Compute the horvitz thompson estimate
HT.estimate<-function(weights,outcome){
  nas <- is.na(weights)|is.na(outcome)
  if(is.factor(outcome)){
    num <- rep(0,length=length(levels(outcome[!nas])))
    a <- tapply(weights[!nas],as.numeric(outcome[!nas]),sum)
    num[as.numeric(names(a))] <- a
  }else{
    num<-sum(outcome[!nas]*weights[!nas])
  }
  den<-sum(weights[!nas])
  num/den
}




#' Counts the number or recruiter->recruitee transitions between different levels
#' of the grouping variable.
#' @param rds.data An rds.data.frame
#' @param group.variable The name of a categorical variable in rds.data
#' @export
#' @examples 
#' data(faux)
#' count.transitions(faux,"X")
count.transitions <- function(rds.data, group.variable)
{
  stopifnot(is.rds.data.frame(rds.data))
  
  #extract data
  id <- get.id(rds.data)
  recruiter.id <- get.rid(rds.data)
  grp <- as.factor(rds.data[[group.variable]])
  #       Make sure the factor labels are alphabetic!
  grp=factor(grp,levels=levels(grp)[order(levels(grp))])
  
  #count transitions
  ri <- match(recruiter.id,id)
  rgrp <- grp[ri]
  res <- table(rgrp,grp)
  for(i in 1:nrow(res)){
    if(sum(res[i,])==0){res[i,] <- 0.1}
  }
  class(res) <- "matrix"
  res
}




#' Markov chain statistionary distribution
#' @param mle The transition probabilities
#' @export
#' @return A vector of proportions representing the proportion in each group at the
#' stationary distribution of the Markov chain.
get.stationary.distribution <- function(mle){
  stat.dist <- NULL	
  eigen.analysis <- eigen(t(mle))	
  
  # If -1 is an eigenvalue then we have a periodic chain ... which is bad.
  values.close.to.minus.one <- abs(eigen(t(mle))$values + 1.0) < 
    sqrt(.Machine$double.eps)
  if(sum(values.close.to.minus.one) == 1){
    cat(" Markov chain appears to be periodic.\n")
    stat.dist <- rep(NA,nrow(mle))
    attr(stat.dist,'status') <- "Periodic Markov MLE"
    return(stat.dist)
  }
  
  # The stationary distribution will be given by
  # the eigenvector with eigenvalue equal to one.    
  values.close.to.one <- abs(eigen(t(mle))$values - 1.0) <
    sqrt(.Machine$double.eps)
  
  # If there is exactly one of these, then we can take it as eigenvalue corresponding to the
  # stationary distribution.
  if(sum(values.close.to.one) == 1){
    stationary.dist.eigen.vector <- eigen.analysis$vectors[,which(values.close.to.one)]
    stationary.dist.eigen.vector <- Re(stationary.dist.eigen.vector)
    stat.dist <- stationary.dist.eigen.vector/sum(stationary.dist.eigen.vector)
    attr(stat.dist,'status') <- "MLE gives unique stationary distribution."
    names(stat.dist) <- colnames(mle)	
  }
  return(stat.dist)	
}





