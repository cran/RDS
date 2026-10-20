

#' Compute estimates of the sampling weights of the respondent's observations
#' based on various estimators
#' 
#' 
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"},
#' and \code{"Arithemic Mean"}. It defaults to \code{"Gile's
#' SS"}.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing, the weights will sum to 1. Note that
#' this parameter is required for Gile's SS.
#' @param subset An optional criterion to subset \code{rds.data} by. It is
#' an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{subset = seed > 0} to
#' exclude seeds. It can also be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param control A list of control parameters for algorithm
#' tuning. Constructed using\cr
#' \code{\link{control.rds.estimates}}.
#' @param outcome.variable The name of the variable defining the groups. Required
#' for the group-model weight types ("RDS-I", "RDS-I (DS)" and "HCG").
#' @param hcg.wave.fallback Logical. Internal use: allows the HCG weights to fall
#' back to recruitment-wave ordering when recruitment times are absent. Used by
#' the bootstrap machinery; direct HCG weight computations require recruitment
#' time unless \code{small.fraction=TRUE}.
#' @param cont.breaks The number of groups a continuous outcome is discretized
#' into for the group-model weight types.
#' @param ... Additional parameters passed to the individual weighting algorithms.
#' @return A vector of weights for each of the respondents. It is of the same
#' size as the number of rows in \code{rds.data}. Subjects with missing outcome
#' values receive weight 0 for the group-model weight types; they are retained
#' internally (as their own group) while the weights are fit, so that
#' recruitment chains passing through them are preserved.
#' @seealso \code{\link{rds.I.weights}}, \code{\link{gile.ss.weights}}, \code{\link{vh.weights}}
#' @export
compute.weights <- function(rds.data,
                            weight.type = c("Gile's SS","RDS-I","RDS-I (DS)","RDS-II","Arithmetic Mean","HCG"),
                            N = NULL, subset=NULL, control=control.rds.estimates(),
                            outcome.variable=NULL, hcg.wave.fallback=FALSE,
                            cont.breaks=3, ...){
  args <- list(...)
  if(!is.null(N) && N < nrow(rds.data)){
    stop(sprintf("The population size, %d, is less than the sample size, %d. The population size must be at least as large as the sample size for the estimate to make sense.",N,nrow(rds.data)))
  }

  n <- nrow(rds.data)
  if(is.null(attr(rds.data,"network.size.variable")))
    stop("rds.data must have a network.size attribute.") #wait-is network.size a vector?
  init.deg <- as.numeric(rds.data[[attr(rds.data,"network.size.variable")]])
  deg <- init.deg
  remvalues <- deg==0 | is.na(deg)
  if(any(remvalues)){
    warning(paste(sum(remvalues),"of",nrow(rds.data),
                  "network sizes were missing or zero. The estimator will presume these are",max(deg,na.rm=TRUE)), call. = FALSE)
    deg[remvalues] <- max(deg,na.rm=TRUE)
  }
  if(missing(weight.type)){
    weight.type <- "Gile's SS"
  }
  weight.type <- match.arg(weight.type,
                           c("Gile's SS","RDS-I", "RDS-II", "RDS-I (DS)","Arithmetic Mean","HCG"))
  if(is.na(weight.type)) { # User typed an unrecognizable name
    stop(paste('You must specify a valid weight.type. The valid types are "Gile\'s SS","RDS-I", "RDS-II", "RDS-I (DS)", "HCG", and "Arithmetic Mean"'), call.=FALSE)
  }
  if(weight.type %in% c("Gile's SS") && is.null(N)){
    stop(paste(weight.type,"estimator requires an estimated population size (N)"))
  }
  small.fraction <- "small.fraction" %in% names(args) && isTRUE(args$small.fraction)
  if(weight.type %in% c("HCG")){
    if(is.null(N) && !small.fraction){
      stop(paste(weight.type,"estimator requires an estimated population size (N)"))
    }
    if(!small.fraction && !hcg.wave.fallback && !has.recruitment.time(rds.data)){
      stop(paste(weight.type,"weights require subject recruitment time unless",
                 "small.fraction=TRUE. Set the recruitment time attribute of the data."))
    }
  }

  # Evaluate any subset expression against the original (un-recoded) data.
  se <- substitute(subset)
  if(!is.null(se)){
    subset <- eval(se, rds.data, parent.frame())
    subset[is.na(subset)] <- FALSE
  }

  # Boundary: for the group-model weight types, subjects with missing
  # outcomes are recoded into an internal group so the recruitment tree
  # stays intact while the weights are fit; their weights are zeroed below.
  uses.outcome <- weight.type %in% c("RDS-I","RDS-I (DS)","HCG")
  excl <- NULL
  if(uses.outcome){
    if(is.null(outcome.variable))
      stop(paste(weight.type, "weights require outcome.variable"))
    if(is.null(rds.data[[outcome.variable]]))
      stop(paste("outcome.variable", outcome.variable, "not present in data"))
    out0 <- rds.data[[outcome.variable]]
    if(is.numeric(out0)){
      # Match RDS.estimates.local: continuous outcomes are discretized into
      # cont.breaks groups for the group-model weights; excluded subjects
      # form the internal group of the cut factor.
      excl <- is.na(out0)
      if(!is.null(se))
        excl <- excl | !subset
      rds.data[[outcome.variable]] <- suppressWarnings(.cut2.excl(out0, g=cont.breaks, mask=excl))
    }else{
      if(!is.factor(out0))
        out0 <- as.factor(out0)
      excl <- .excluded.mask(out0, if(is.null(se)) NULL else subset)
      if(any(excl))
        rds.data[[outcome.variable]] <- .mint.excluded(out0, excl)
    }
  }

  weights <- switch(weight.type,
                    `RDS-I` = .rds.I.weights.core(rds.data,outcome.variable=outcome.variable,N=N,...),
                    `RDS-I (DS)` = .rds.I.weights.core(rds.data,outcome.variable=outcome.variable,N=N,smoothed=TRUE,...),
                    `RDS-II` = vh.weights(degs = deg, N=N),
                    `Arithmetic Mean` = rep(ifelse(is.null(N),1,N)/n, n),
                    `HCG` = .hcg.weights.core(rds.data, outcome.variable=outcome.variable, N = N,
                      reltol=control$hcg.reltol, max.optim=control$hcg.max.optim, ...),
                    `Gile's SS` = gile.ss.weights(degs = deg, N = N, SS.infinity=control$SS.infinity, ...)
  )
  if(uses.outcome && any(excl))
    weights[excl] <- 0
  if(!is.null(se)){
    a <- weights[subset]
    if(weight.type=="Gile's SS"){
      attr(a,"N") <- attr(weights,"N")
      attr(a,"estimateN") <- attr(weights,"estimateN")
    }
    weights <- a
  }
  return(weights)
}



#' RDS-I weights
#' @param rds.data An rds.data.frame
#' @param outcome.variable The variable used to base the weights on.
#' @param N Population size
#' @param smoothed Should the data smoothed RDS-I weights be computed.
#' @param ... Additional parameters passed to \code{\link{compute.weights}}.
#' @details This function dispatches to \code{\link{compute.weights}}, which
#' handles missing outcome values by retaining the subjects internally (as
#' their own group, preserving recruitment chains) and returning weight 0
#' for them.
#' @export
rds.I.weights<-function(rds.data, outcome.variable, N=NULL,smoothed=FALSE, ...){
  compute.weights(rds.data,
                  weight.type = if(smoothed) "RDS-I (DS)" else "RDS-I",
                  N = N, outcome.variable = outcome.variable, ...)
}

.rds.I.weights.core<-function(rds.data, outcome.variable, N=NULL,smoothed=FALSE, ...){
  if(is.null(rds.data[[outcome.variable]])){
    stop(paste("RDS-I outcome.variable", outcome.variable,"not present in data"))
  }
  .assert.complete.outcome(rds.data[[outcome.variable]])
  tij <- count.transitions(rds.data, outcome.variable)
  network.size <- attr(rds.data, "network.size.variable")
  remvalues <- !is.na(rds.data[[network.size]]) & (rds.data[[network.size]] > 0)
  if(sum(remvalues) < nrow(rds.data)){
    warning(paste(nrow(rds.data)-sum(remvalues),"of",nrow(rds.data),
                  "network size values were missing and were removed."), call. = FALSE)
    rds.data[[network.size]][!remvalues] <- max(rds.data[[network.size]],na.rm=TRUE)+1
  }
  if(ncol(tij)>1){
    smoothed.tij <- matrix(nrow=0,ncol=0)
    markov.mle <- prop.table(tij, margin=1)
    q.hat <- get.stationary.distribution(markov.mle)
    if(smoothed){
      #			Demographic adjustment of raw transition counts
      smoothed.tij <- sum(tij)*q.hat*markov.mle
      smoothed.tij <- 0.5 * (smoothed.tij + t(smoothed.tij))
      markov.mle <- prop.table(smoothed.tij, margin=1)
      q.hat <- get.stationary.distribution(markov.mle)
    }
    h.hat <- get.h.hat(rds.data,outcome.variable,network.size)    
    est <- as.list((q.hat/h.hat)/sum(q.hat/h.hat))
    
    d <- tapply(rds.data[[outcome.variable]],rds.data[[outcome.variable]],length)
    d <- d[match(names(est), names(d))]
    f <- match(rds.data[[outcome.variable]],names(est))
    weights <- rep(0,length(f))
    weights[!is.na(f)] <- as.numeric(unlist(est[f]))
    weights <- as.numeric(weights/d[f])
    weights[is.na(weights)] <- 0
  }else{
    weights <- rep(1/length(rds.data[[outcome.variable]]),length(rds.data[[outcome.variable]]))
  }
  if(!is.null(N)){
    weights <- N * weights / sum(weights)
  }
  weights
}


#' Volz-Heckathorn (RDS-II) weights
#' @param degs The degrees (i.e. network sizes) of the sample units.
#' @param N Population size
#' @export
vh.weights<-function(degs, N=NULL){
  if(is.null(degs)){
    return(NULL)
  }
  degs[degs==0]<-1
  isnadegs <- is.na(degs)
  degs <- sum(!isnadegs)*(degs)/sum(degs,na.rm=TRUE)
  weights <- (1/degs)
  weights[is.na(weights)] <- 0
  if(!is.null(N)){
    weights <- N * weights / sum(weights)
  }
  weights
}

#was spps_weights
#' Weights using Giles SS estimator
#' @param degs subjects' degrees (i.e. network sizes).
#' @param N Population size estimate.
#' @param number.ss.samples.per.iteration The number of samples to use to estimate inclusion probabilities
#' in a probability proportional to size without replacement design.
#' @param number.ss.iterations number of iterations to use in giles SS algorithm.
#' @param hajek Should the hajek estiamtor be used. If false, the HT estimator is used.
#' @param SS.infinity The sample proportion, \code{n/N}, below which the computation of the SS weights should simplify to that of the \code{RDS-II} weights.
#' @param se Should covariances be included.
#' @param ... unused
#' @export
gile.ss.weights<-function(degs,N,number.ss.samples.per.iteration=500,number.ss.iterations=5,
                          hajek=TRUE,SS.infinity=0.04,se=FALSE,...){
  if(is.null(degs)){
    return(NULL)
  }
  degs[degs==0]<-1
  isnadegs <- is.na(degs)
  n <- length(degs[!isnadegs])
  if(n==0){
    return(NULL)
  }
  if(n==1){
    weights <- rep(0,length(degs))
    weights[!isnadegs] <- N
    return(weights)
  }
  if(N < n){
    stop(sprintf("Population size N (%s) cannot be smaller than the sample size (%s).", N, n))
  }
  mapping<-getestCstacked(samp=degs[!isnadegs],n=N,nit=number.ss.iterations,
                          nsampsamp=number.ss.samples.per.iteration,trace=FALSE,hajek=hajek,SS.infinity=SS.infinity)

  weights <- degs
  if(length(mapping$classes) == 1)
    pis <- degs[!isnadegs]
  else
    pis <- stats::approx(x=mapping$classes,y=1/mapping$probs,xout=degs[!isnadegs],rule=2)$y
  weights[!isnadegs] <- pis
  weights[is.na(weights)] <- 0
  weights <- mapping$n*weights/sum(weights,na.rm=TRUE)
  attr(weights,"N") <- N
  attr(weights,"estimateN") <- mapping$n
  if(se){
    # Set up for covariances
    theta <- mapping$props
    covtheta <- -outer(theta, theta) / mapping$n
    diag(covtheta) <- (1-theta)*theta / mapping$n
    cov <- covtheta
    cov[is.na(cov)] <- 0
    cov <- cov / outer(mapping$n*theta*mapping$probs, mapping$n*theta*mapping$probs)
    #   covariances of the (normalized) weights
    md <- match(degs,mapping$classes)
    wcov <- outer(md,md,function(i,j){cov[cbind(i,j)]})
    attr(weights,"cov") <- wcov
    attr(weights,"classes") <- degs
  }
  weights
}




#' homophily configuration graph weights
#' @param rds.data An rds.data.frame
#' @param outcome.variable The variable used to base the weights on.
#' @param N Population size
#' @param small.fraction should a small sample fraction be assumed
#' @param reltol Relative convergence tolerance for the HCG estimator.  The algorithm stops if
#' it is unable to reduce the log-likelihood by a factor of \code{reltol * (abs(log-likelihood) + reltol)}
#' at a step. Defaults to \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.
#' @param max.optim The number of iterations on the likelihood optimization for the HCG estimator.
#' @param theta.start The initial value of theta used in the likelihood optimization for the HCG estimator. If NULL, the default, it is the margin of the table of counts for the transitions.
#' @param weights.include.seeds logical Should the weights be computed including the influence of the seeds?
#' @param hcg.wave.fallback Logical. Internal use: allow recruitment-wave ordering
#' when recruitment times are absent (used by the bootstrap machinery). Direct
#' weight computations require recruitment time unless \code{small.fraction=TRUE}.
#' @param ... Additional parameters passed to \code{\link{compute.weights}}.
#' @details This function dispatches to \code{\link{compute.weights}}, which
#' handles missing outcome values by retaining the subjects internally (as
#' their own group, preserving recruitment chains) and returning weight 0
#' for them.
#' @examples 
#' data(fauxtime)
#' hcg.weights(fauxtime,"var1",N=3000)
#' fauxtime$NETWORK[c(1,100,40,82,77)] <- NA
#' @export
hcg.weights<-function(rds.data, outcome.variable, N=NULL,small.fraction=FALSE,
                      reltol=sqrt(.Machine$double.eps), max.optim=500, theta.start=NULL,
                      weights.include.seeds = TRUE, hcg.wave.fallback=FALSE, ...){
  control <- control.rds.estimates(hcg.reltol=reltol, hcg.max.optim=max.optim)
  compute.weights(rds.data, weight.type="HCG", N=N, control=control,
                  outcome.variable=outcome.variable,
                  hcg.wave.fallback=hcg.wave.fallback,
                  small.fraction=small.fraction, theta.start=theta.start,
                  weights.include.seeds=weights.include.seeds, ...)
}

.hcg.weights.core<-function(rds.data, outcome.variable, N=NULL,small.fraction=FALSE,
                      reltol=sqrt(.Machine$double.eps), max.optim=500, theta.start=NULL,
                      weights.include.seeds = TRUE, ...){
  if(is.null(rds.data[[outcome.variable]])){
    stop(paste("HCG outcome.variable", outcome.variable,"not present in data"))
  }
  out <- rds.data[[outcome.variable]]
  .assert.complete.outcome(out)
  time <- get.recruitment.time(rds.data, to.numeric=FALSE, wave.fallback = TRUE)

  degree <- get.net.size(rds.data)
  obs.levels <- setdiff(unique(as.character(out)), .rds.excluded)
  if(length(obs.levels) <= 1){
    if(is.null(N)){
      weights <- vh.weights(degree, N=N)
    }else{
      weights <- gile.ss.weights(degs = degree, N = N)
    }
    return(weights)
  }
  degree <- .impute.degree.group.median(degree, out)

  hcg.result <- hcg.estimate(get.id(rds.data), get.rid(rds.data), time,
                      degree, out, N, small.fraction, reltol=reltol,
                      max.optim=max.optim, theta.start=theta.start, weights.include.seeds=weights.include.seeds)

  weights <- hcg.result$weights
  if(!is.null(N)){
    weights <- N * weights / sum(weights)
  }
  weights
}
