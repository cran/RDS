HCG.bootstrap <- function(rds.data, group.variable, number.of.bootstrap.samples, N, fun=function(x) x, 
                         small.fraction=is.null(N) || nrow(rds.data) / N < 0.05, 
                         cont.breaks=3, hcg.est=NULL, control=control.rds.estimates(), verbose=TRUE, ...){
  
  if(is(rds.data,"rds.data.frame")){
    stopifnot(group.variable %in% names(rds.data))
    network.size <- attr(rds.data, "network.size.variable")
  }else{
    stop("rds.data must be of type rds.data.frame")
  }
  is.cts <- is.numeric(rds.data[[group.variable]])
  n <- nrow(rds.data)
  

  #set N very large in case of infinite population
  if(is.null(N)){
    N <- n * 1000000
    small.fraction <- TRUE
  }

  wave <- get.wave(rds.data)

  
  if(!has.recruitment.time(rds.data)){
    #if(!small.fraction)
    #  warning("Performing HCG bootstrap on data without recruitment time where the smaple fraction is not insignificant. Recruitment wave ordering will be used.")
    time <- rep(1, n)
  }else{
    time <- get.recruitment.time(rds.data)  
  }
  rds <- rds.data[order(time, wave),]
  if(!has.recruitment.time(rds)){
    time <- rep(1, n)
  }else{
    time <- get.recruitment.time(rds)  
  }
  
  rid <- get.rid(rds)
  id <- get.id(rds)
  degree <- get.net.size(rds)
  if(is.cts){
    out <- .cut2(rds[, group.variable], g=cont.breaks)
    
    #create a map in each outcome group between degrees and continuous outcome values
    deg.out.map <- list()
    for(i in 1:cont.breaks){
      d <- stats::na.omit(degree[as.numeric(out) == i])
      o <- stats::na.omit(rds[as.numeric(out) == i, group.variable])
      lis <- list()
      for(j in seq_along(d)){
        if(length(lis) <= d[j] || is.null(lis[[d[j]]])){
          tr <- try(lis[[d[j]]] <- o[j])
		  if(inherits(tr,"try-error"))
			browser()
	    }else
          lis[[d[j]]] <- c(lis[[d[j]]], o[j])
      }
      deg.out.map[[i]] <- lis
    }
  }else
    out <- as.factor(rds[, group.variable])
  outNum <- as.numeric(out)
  recruiter.ind <- match(rid, id)
  seed <- which(!(rid %in% id))
  if(is.null(hcg.est)){
    hcg <- hcg.estimate(id, rid, time, degree, out, N, small.fraction=small.fraction)
    theta <- hcg$theta
    weights <- hcg$weights
    yhat <- hcg$yhat
  }else{
    theta <- hcg.est$theta
    weights <- hcg.est$weights
    yhat <- hcg.est$yhat
  }
  lev <- levels(out)
  nlev<- length(lev)
  
  deg.by.group <- list()
  edge.end.by.group <- rep(NA, nlev)
  for(i in 1:nlev){
    wts <- weights[out==lev[i]]
    deg <- degree[out==lev[i]]
    tab <- wtd.table(deg, wts, normwt=TRUE)
    deg.count <- tab$sum.of.weights * yhat[i] * N / sum(tab$sum.of.weights)
    deg.value <- tab$x
    edge.end.by.group[i] <- round(sum(deg.count*deg.value))
    deg.by.group[[i]] <- data.frame(deg.value, deg.count)
  }
  
  out.na <- is.na(out)
  
  results <- list()
  doboot <- function(){
    running.deg.by.group <- deg.by.group
    samp.edge.end.by.group <- rep(0, nlev)
    out.boot <- rep(NA, n)
    deg.boot <- rep(NA, n)
    
    out.boot[seed] <- sample.int(nlev, size=length(seed), prob=yhat, replace = TRUE)#outNum[seed]#
    for(s in seed){
      d <- na.omit(degree[out == levels(out)[out.boot[s]]])
      deg.boot[s] <- sample(d,size = 1,prob = d)
    }
    #if(any(is.na(out.boot[seed]))){
    #  isna <- is.na(out.boot[seed])
    #  out.boot[seed][isna] <- sample.int(nlev, size=sum(isna),prob=yhat, replace = TRUE)
    #}
    #deg.boot[seed] <- degree[seed]
    for(i in (1:n)[-seed]){
      from <- out.boot[recruiter.ind[i]]
      if(small.fraction)
        trans.prob <- theta[from,]
      else
        trans.prob <- pmax(0, theta[from,] * (edge.end.by.group - samp.edge.end.by.group) / edge.end.by.group)
      trans.prob[!is.finite(trans.prob)] <- 0 # handle factor levels with no observations
      to <- out.boot[i] <- sample.int(nlev, size = 1, prob = trans.prob)
      deg.ind <- try(sample.int(nrow(deg.by.group[[to]]), size = 1, 
                                prob=running.deg.by.group[[to]]$deg.count * running.deg.by.group[[to]]$deg.value), silent=TRUE)
      if(inherits(deg.ind, "try-error")){
        d <- na.omit(degree[out == to])
        deg.boot[i] <- sample(d,size = 1,prob = d)
      }else
        deg.boot[i] <- running.deg.by.group[[to]]$deg.value[deg.ind]
      
      if(!small.fraction){
        # remove observed node
        running.deg.by.group[[to]]$deg.count[deg.ind] <- max(0, running.deg.by.group[[to]]$deg.count[deg.ind] - 1)
        samp.edge.end.by.group[to] <- samp.edge.end.by.group[to] + deg.boot[i]
      }
    }
    if(is.cts){
      #Within group/degree, individuals are equally likely to be recruited
      out.boot.cont <- rep(NA, length(out.boot))
      for(r in 1:length(out)){
        outs <- deg.out.map[[out.boot[r]]][[deg.boot[r]]]
        out.boot.cont[r] <- outs[sample.int(length(outs), size=1)]
      }
      out.boot <- out.boot.cont
    }else{
      out.boot <- as.factor(out.boot)
      levels(out.boot) <- lev
    }
    boot.rds <- rds
    
    out.boot[out.na] <- NA
    
    boot.rds[[group.variable]] <- out.boot
    boot.rds[[attr(boot.rds,"network.size.variable")]] <- deg.boot
    #results[[j]] <- fun(boot.rds)
    fun(boot.rds)
  }
  for(j in 1:number.of.bootstrap.samples){
    for(k in 1:10){
      tr <- try(results[[j]] <- doboot(), silent = TRUE)
      if(!inherits(tr,"try-error"))
        break
    }
    if(inherits(tr,"try-error"))
      stop(tr)

    if(verbose) {
      if(j == trunc(j/(number.of.bootstrap.samples/10))*(number.of.bootstrap.samples/10)){
        cat(paste((100/10)*trunc(j/(number.of.bootstrap.samples/10)),'% completed ...\n',sep=""))
      }
    }
  }
  results
}


HCG.bootstrap.se <- function(rds.data, group.variable,
                                  number.of.bootstrap.samples,estimator.name,N=NULL, 
                             small.fraction=is.null(N) || nrow(rds.data) / N < 0.05,
                             to.factor=FALSE, control=control.rds.estimates(), cont.breaks=3, ...){
  estimate <- function(boot,theta){
    #if(estimator.name == "HCG"){
    #  hcg.estimate(get.id(boot), get.rid(boot), get.recruitment.time(boot, wave.fallback=TRUE), get.net.size(boot), 
    #               boot[[group.variable]], N, small.fraction=small.fraction)$yhat
    #}else{
      if(!missing(theta))
        control$hcg.theta.start <- theta
      RDS.estimates.local(
        rds.data=boot,
        outcome.variable=group.variable,
        weight.type=estimator.name,
        empir.lik=FALSE,
        N=N,
        to.factor=to.factor,
        cont.breaks=cont.breaks,
        control=control,
        ...)@estimate
    #}
  }
  if(to.factor){
    rds.data[[group.variable]] <- as.factor(rds.data[[group.variable]])
  }
  
  control$hcg.reltol <- control$hcg.BS.reltol
  result <- HCG.bootstrap(rds.data, group.variable, number.of.bootstrap.samples, N, fun=estimate,
                         small.fraction=small.fraction, cont.breaks=cont.breaks, control=control, ...)
  result <- do.call(rbind, result)
  
  result <- result[apply(!is.nan(result),1,all), , drop=FALSE]
  if(nrow(result)>1){
    a=(sqrt(diag(stats::var(result))))	
  }else{
    a=(sqrt(stats::var(as.numeric(result))))	
  }
  attr(a,"bsresult") <- result
  return(a)
}


