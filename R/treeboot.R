

treeboot <- function(seed.rows, recruits, network.size, outcomes, population.size, fixed.size=TRUE){
  resample <- function(x, size=length(x), replace=TRUE){
    x[sample.int(length(x), size=size, replace=replace)]
  }
  
  n <- length(recruits)
  nseed <- length(seed.rows)
  
  samp.rows <- as.list(resample(seed.rows))
  samp.ids <- as.list(1:length(seed.rows))
  
  sample.id <- function(node.row, id, recr.id, wave){
    if(is.null(node.row)){
      recr <- resample(seed.rows)
      result <- list()
    }else{
      recr <- resample(recruits[[node.row]])
      result <- list(c(node.row, id, recr.id, wave))
    }
    for(i in seq_along(recr)){
      result <- c(result, sample.id(recr[i], paste0(id,i), id, wave + 1))
    }
    result
  }
  if(fixed.size){
    result <- list()
    boot.seeds <- 0
    boot.n <- 0
    while(TRUE){
      seed <- resample(seed.rows, 1)
      tree <- sample.id(seed, as.character(boot.seeds), "_", 0)
      if(boot.n + length(tree) > n){
        wave <- sapply(tree,function(x) as.numeric(x[4]))
        for(i in max(wave):1){
          ntrim <- boot.n + length(tree) - n
          nw <- sum(wave==i)
          if(nw < ntrim){
            tree[wave==i] <- NULL
            wave <- sapply(tree,function(x) as.numeric(x[4]))
          }else{
            ind <- resample(which(wave==i), ntrim, replace=FALSE)
            tree[ind] <- NULL
            break
          }
        }
      }
      result <- c(result, tree)
      boot.n <- length(result)
      boot.seeds <- boot.seeds + 1
      if(length(result) >= n)
        break
    }
  }else{
    result <- sample.id(NULL,"", NA, -1)
  }
  df <- data.frame(row=sapply(result,function(x) as.numeric(x[1])),
                   id=sapply(result,function(x) x[2]),
                   recruiter.id=sapply(result,function(x) x[3]))
  df <- cbind(df,outcomes[df$row, ,drop=FALSE])
  df$network.size.variable <-  network.size[df$row]
  df$row <- NULL
  bootstrapped.data <- as.rds.data.frame(df,
                                         population.size=population.size,
                                         check.valid=FALSE) 
  bootstrapped.data
}

treeboot.bootstrap.estimates <- function(rds.data, group.variable, number.of.bootstrap.samples,
                                         estimator.name, N=NULL, fixed.size=TRUE, ...){
  
  if(is(rds.data,"rds.data.frame")){
    stopifnot(group.variable %in% names(rds.data))
    network.size <- attr(rds.data, "network.size.variable")
  }else{
    stop("rds.data must be of type rds.data.frame")
  }
  
  # First extract the needed information from the rds data frame object. 
  group.memberships <- as.vector(rds.data[[group.variable]])
  group.names <- unique(group.memberships)
  
  # NB: if one wants to treat missing values as a separate group,
  # then the following will need to be changed. 
  group.names <- group.names[!is.na(group.names)]
  
  id <- get.id(rds.data)
  recruiter.id <- get.rid(rds.data)   
  
  # Now rationalize the recruiter.id information so that it corresponds to recruiter
  # row.
  recruiter.row <- match(recruiter.id,id)
  # The zeros in the recruiter id will be mapped to NA's
  recruiter.row[is.na(recruiter.row)] <- 0
  n <- length(id)
  id.row <- 1:n
  seed.rows <- which(get.wave(rds.data) == 0)
  
  recruits <- lapply(id.row, function(i) which(recruiter.row==i))
  outcomes <- data.frame(rds.data)[group.variable]
  ###############################################################################
  # This internal function creates a single bootstrap sample.  The return value 
  # is a data frame with id, recruiter.id, group variable and network size.  
  #bootstrapper <- function(){
  #  treeboot(seed.rows, recruits, rds.data[[network.size]], 
  #           outcomes, get.population.size(rds.data), fixed.size)
  #}
  
  f <- function(){
    boot <- treeboot(seed.rows, recruits, rds.data[[network.size]], 
                     outcomes, get.population.size(rds.data), fixed.size)
    RDS.estimates.local(
      rds.data=boot,
      outcome.variable=group.variable,
      weight.type=estimator.name,
      empir.lik=FALSE,
      N=N,
      ...)@estimate
  }
  bs.results <- replicate(number.of.bootstrap.samples, f())
  
  value <- matrix(0,ncol=length(group.names),nrow=number.of.bootstrap.samples)
  colnames(value) <- group.names
  colnames(value)[colnames(value)=="NA.NA"] <- "NA"
  if(is.matrix(bs.results)){
    for(i in 1:nrow(bs.results)){
      value[,i] <-  unlist(bs.results[i,])
    }
  }else{
    if(is.list(bs.results)){
      for(i in 1:number.of.bootstrap.samples){
        value[i,names(bs.results[[i]])] <-  unlist(bs.results[[i]])
      }
    }else{
      value[,1] <- bs.results
      value[,2] <- bs.results
    }
  }
  
  return(value) 
  
}




treeboot.bootstrap.se <- function(rds.data,group.variable,
                                  number.of.bootstrap.samples,estimator.name,N=NULL,...){
  result <- treeboot.bootstrap.estimates(
    rds.data=rds.data,
    group.variable=group.variable,
    number.of.bootstrap.samples=number.of.bootstrap.samples,
    estimator.name=estimator.name,
    N=N,
    ...)
  
  result <- result[apply(!is.nan(result),1,all),]
  if(nrow(result)>1){
    a=(sqrt(diag(stats::var(result))))	
  }else{
    a=(sqrt(stats::var(as.numeric(result))))	
  }
  attr(a,"bsresult") <- result
  return(a)
}
