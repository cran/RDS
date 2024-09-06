
hcg.estimate <- function(subject, recruiter, time, degree, out, N, small.fraction=FALSE, tol=.00005, max.iter=50, 
                         reltol=sqrt(.Machine$double.eps), BS.reltol=sqrt(.Machine$double.eps), max.optim=500, 
                         theta.start=NULL, weights.include.seeds = TRUE){
  rds <- data.frame(subject, recruiter, time, degree, out)
  rds <- rds[order(rds$time),]
  rInd <- match(rds$recruiter, rds$subject)
  seed <- !(recruiter %in% subject)
  rrInd <- match(rds$recruiter[rInd],rds$subject)
  out <- as.factor(rds$out)
  outNum <- as.numeric(out)
  outNotMiss <- !is.na(out)
  dg <- rds$degree
  rOut <- out[rInd]
  rOutNum <- as.numeric(rOut)
  rOutNotMiss <- !is.na(rOut)
  rrOut <- out[rrInd]
  rd <- rds$degree[rInd] 
  tots <- table(out)
  totSampDeg <- tapply(dg,out,sum)
  totSampDeg[is.na(totSampDeg)] <- 0
  #takes the transition probabilities (theta) in matrix form and makes a vector
  flattenTheta <- function(theta){
    ttmp <- rep(0,nlev*(nlev-1))
    for(i in 1:nlev)
      ttmp[((i-1)*(nlev-1)+1):(i*(nlev-1))] <- theta[i,(1:nlev)[-i]]
    ttmp
  }
  
  #Takes theta in vector form and makes a matrix
  unflattenTheta <- function(tflat){
    theta <- matrix(0,nrow=nlev,ncol=nlev)
    for(i in 1:nlev)
      theta[i,(1:nlev)[-i]] <- tflat[((i-1)*(nlev-1)+1):(i*(nlev-1))]
    diag(theta) <- 1 - rowSums(theta)
    theta
  }
  
  #calculate ybar with theta and dbar
  calcEst <- function(theta, dbar){
    ybar <- 1 / sapply(1:nlev,function(j) sum(theta[j,]*dbar[j] / (theta[,j]*dbar),na.rm=TRUE))
    ybar <- ybar / sum(ybar)
    
    #upper and lower bounds for ybar
    ub <- rep(1, length(ybar))
    lb <- pmax(as.vector(tots) / N, totSampDeg / (N*dbar))
    #If estimate violates contraints, find a ybar that satisfies them
    if(sum(lb) > 1)
      warning("Algorithm: Current dbar not consistant with population size")
    if(anyNA(lb))
      stop("Algorithm: Error in calculating lower bound")
    if(any(ybar < lb)){
      opt <- function(par) sum((par - ybar)^2) + 10*(sum(par) - 1)^2
      grad <- function(par) 2 * (par - ybar) + 20*(sum(par) - 1)
      eqFun <- function(par) sum(par)
      ui <- matrix(0,length(ybar),length(ybar))
      diag(ui) <- 1
      tr <- try(par <- stats::constrOptim((ub+lb)/2, f=opt,grad=grad, ui=ui, ci=lb)$par, silent = TRUE)
      if(inherits(tr,"try-error")){
        #warning("Unable to find par value satisfying constraints.")
        #print(tr)
        return(ybar)
      }
      #cat("\nsum:", sum(par))
      return(par)
    }else{
      return(ybar)
    }
  }
  
  # The log likelihood
  llik <- function(tflat){
    theta <- unflattenTheta(tflat)
    
    if(any(theta<=0) | any(theta>=1)){
      #warning("theta out of range")
      return(1000000000*(sum(abs(theta[theta<0]-1)) + sum(theta[theta>1]) + 1))
    }
    #print(theta)
    ll <- 0
    ybar <- calcEst(theta,dbar)
    Ns <- ybar*N
    ds <- dbar
    unsmp <- Ns*ds - totSampDeg
    if(any(unsmp<0)){
      warning("likelihood evaluated at inconsistant ybar")
      return(1000000000*(sum(-unsmp[unsmp<0]) + 1))
    }
    for(i in 1:length(out)){
      if(i==1){
        s <- rep(0,nlev)
      }else{
        s[outNum[i-1]] <- s[outNum[i-1]] + dg[i-1]
      }
      if(outNotMiss[i] && rOutNotMiss[i]){
        grp <- rOutNum[i]
        oGrp <- outNum[i]
        u <- Ns*ds - s
        num <- theta[grp,oGrp] * u[oGrp] / (Ns[oGrp]*ds[oGrp])
        denom <- sum(theta[grp,] *  u / (Ns*ds))
        t <- num / denom
        if(is.na(t) || !is.finite(t) || t<=0 ||t>1)
          next #browser()
        ll <- ll + log(t)
      }
    }
    -ll
  }
  t <- table(rOut,out)
  nlev <- ncol(t)
  if(!small.fraction && nlev > 10){
    warning("HCG can not be calculated for factors with more than 10 levels. Falling back to small sample fraction approximation")
    small.fraction <- TRUE
    stop()
  }
  if(small.fraction){
    if(any(t==0))
      t <- t + 1
    theta <- prop.table(t,1)
    dbar <- tapply(dg[!seed],out[!seed],function(x) sum(length(x)) / sum(1/x))
    dbar[is.na(dbar)] <- sum(length(dg)) / sum(1/dg)
    yhat <- 1 / sapply(1:nlev,function(j) sum(theta[j,]*dbar[j] / (theta[,j]*dbar),na.rm=TRUE))
    yhat <- yhat / sum(yhat)
    weights <- rep(0,length(out))
    for(i in 1:nlev){
      w <- (1 / dg[!seed & outNum==i & outNotMiss])
      w <- w / sum(w)
      weights[outNum==i & !seed & outNotMiss] <- (w * yhat[i])
    }
    orgOrder <- match(subject, rds$subject)
    names(yhat) <- levels(out)
    return(list(yhat=yhat,
         weights=weights[orgOrder] / sum(weights), 
         dbar=dbar, 
         theta=theta))
  }
  t <- t + 1
  if(is.null(theta.start)){
    theta <- prop.table(t,1)
  }else{
    theta <- theta.start
  }
  dbar <- tapply(dg[!seed & outNotMiss],out[!seed & outNotMiss],function(x) length(x) / sum(1/x))
  if(anyNA(dbar)){
    dbar[is.na(dbar)] <- sum(length(dg)) / sum(1/dg, na.rm=TRUE)
  }
  yhat <- 1 / sapply(1:nlev,function(j) sum(theta[j,]*dbar[j] / (theta[,j]*dbar),na.rm=TRUE))
  yhat <- yhat / sum(yhat)
  yhatLast <- Inf
  for(i in 1:max.iter){
    #print(yhat)
    wts <- lapply(1:nlev,function(i) {
      dd <- dg[outNotMiss & outNum==i & !seed]
      if(length(unique(dd)) == 1)
        return(rep(1, sum(outNotMiss & outNum==i & !seed)))
      w <- gile.ss.weights(dd, max(tots[i],ceiling(yhat[i]*N)))
      if(is.null(w))
        w <- rep(1, sum(outNotMiss & outNum==i & !seed))
      w
    })
    dbar <- sapply(1:nlev,function(i) sum(dg[outNotMiss & outNum==i & !seed]*wts[[i]], na.rm=TRUE)/sum(wts[[i]], na.rm=TRUE))
    if(anyNA(dbar)){
      dbar[is.na(dbar)] <- sum(length(dg)) / sum(1/dg, na.rm=TRUE)
    }
    dbar <- pmax(dbar, totSampDeg / (yhat*N))
#   browser()
    tflat <- flattenTheta(theta)
    tr <- try(opt <- optim(par=tflat,fn=llik,control=list(reltol=reltol,maxit=max.optim)), silent = TRUE)
    if(inherits(tr,"try-error") || opt$convergence != 0){
      tr <- try(opt <- optim(par=tflat,fn=llik,method="BFGS",control=list(reltol=reltol,maxit=max.optim)), silent = TRUE)
      if(inherits(tr,"try-error") || opt$convergence != 0){
        warning("Optimization failed to converge, trying coordinate descent.")
        #browser()
        for(j in 1:nlev){
          inds <- ((i-1)*(nlev-1)+1):(i*(nlev-1))
          opt <- optim(par=tflat[inds],fn=function(x){
            tflat[inds] <- x
            llik(tflat)
          }, control=list(reltol=reltol,maxit=max.optim))
          tflat[inds] <- opt$par
        }
      }
    }
    theta <- unflattenTheta(opt$par)
    yhat <- calcEst(theta,dbar)
    #print(sum((yhat - yhatLast)^2))
    if(sum((yhat - yhatLast)^2) < tol)
      break
    yhatLast <- yhat
  }
  
  if(weights.include.seeds){
    wts <- lapply(1:nlev,function(i) {
      dd <- dg[outNotMiss & outNum==i]
      if(length(unique(dd)) == 1)
        return(rep(1, sum(outNotMiss & outNum==i)))
      w <- gile.ss.weights(dd, max(tots[i],ceiling(yhat[i]*N)))
      if(is.null(w))
        w <- rep(1, sum(outNotMiss & outNum==i))
      w
    })
    weights <- rep(0,length(out))
    for(i in 1:nlev){
      w <- wts[[i]] / sum(wts[[i]])
      weights[ outNotMiss & outNum==i] <- (w * yhat[i])
    }
  }else{
    weights <- rep(0,length(out))
    for(i in 1:nlev){
      w <- wts[[i]] / sum(wts[[i]])
      weights[ outNotMiss & outNum==i & !seed] <- (w * yhat[i])
    }
  }
  weights <- weights / sum(weights)
  orgOrder <- match(subject, rds$subject)
  names(yhat) <- levels(out)
  list(yhat=yhat,
       weights=weights[orgOrder] / sum(weights), 
       dbar=dbar, 
       theta=theta)
}
