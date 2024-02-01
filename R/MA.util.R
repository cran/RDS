getnetest <- function(dissample,
                      degsample,
                      N,
                      todisall,
                      tonodisall,
                      nit = 5,
                      M1 = 50,
                      M2 = 10,
                      seed = 1,
                      short = FALSE,
                      theta0 = -0.79,
                      parallel = 0,
                      parallel.type = "PSOCK",
                      Qstart = NULL,
                      MPLEsamplesize = 50000,
                      SAN.nsteps = 2^19,
                      SAN.maxit = 3,
                      fixinitial = 1,
                      nsamp0 = 10,
                      coupons = 2,
                      interval = 10000,
                      smooth = FALSE,
                      maxdeg = max(degsample),
                      crossStat = NULL,
                      net = NULL,
                      use.net = FALSE,
                      rds.data = NULL,
                      trait.variable = "disease",
                      verbose = TRUE) {
  set.seed(seed)
  n <- length(dissample)
  
  if (!is.null(net)) {
    if (is.null(net %v% "disease")) {
      net %v% "disease" <- as.numeric(net %v% "group")
    }
    crossStat <-
      ergm::summary_formula(net ~ nodemix("disease", levels2 = c(1, 3)))
    degpop <- sapply(net$iel, length) + sapply(net$oel, length)
    maxdeg <- max(degpop)
    dispop <-
      as.numeric(network::get.vertex.attribute(net, "disease"))
    dummy <-
      sort(10 * rep((0:maxdeg), 2) + rep(c(0, 1) , each = (maxdeg + 1)))
    rawCounts <-
      tototab.RDS(degpop, dispop, dummy) #actual network counts
    cat(sprintf("rawCounts:\n"))
    print(rawCounts)
    cat(sprintf("crossStat: %f\n", crossStat))
    cat(sprintf(
      "Population diff. activity= %f\n",
      network.differential.activity(net, 'disease')
    ))
    cat(sprintf(
      "Population homophily= %f\n",
      network.homophily(net, 'disease')
    ))
    a = ergm::summary_formula(net ~ nodefactor('disease', base = 0))
    cat(
      sprintf(
        "ergm::summary_formula(net~nodefactor('disease',base=0)): 0=%f 1=%f\n",
        a[1],
        a[2]
      )
    )
    # Next to stop sampling from passed network
    if (!use.net) {
      net <- NULL
      rawCounts <- NULL
    }
    #     if(!use.net){ net <- NULL}
  } else{
    rawCounts <- NULL
  }
  
  dummy <-
    sort(10 * rep((0:maxdeg), 2) + rep(c(0, 1) , each = (maxdeg + 1)))
  # listing of all possible index values
  
  dummydis <- round(10 * (dummy / 10 - round(dummy / 10)))
  dummydegs <- (dummy - dummydis) / 10
  
  #this is the index of each sampled node within the ordered list of equivalence classes in dummy
  indices <- toindexref.RDS(degsample, dissample, dummy)
  
  #this is the tabulated sample value - counts of nodes of each type
  sampleTab <-
    tototab.RDS(degsample, dissample, dummy) #observed sample counts
  
  seed.indices <-
    toindexref.RDS(degsample[1:nsamp0], dissample[1:nsamp0], dummy)
  
  #Initialize sampling probabilities
  Qstart0 <- round(dummy / 10) * sum(1 / degsample) / N #initial
  
  last <- round(length(degsample) / 3):length(degsample)
  s01s <-
    (dissample == 1) * tonodisall + (dissample == 0) * todisall
  
  if (!is.null(rds.data) & length(degsample) < 0.5 * N) {
    network.size <- attr(rds.data, "network.size.variable")
    if (is.null(network.size)) {
      network.size <- "degree"
    }
    tij <- count.transitions(rds.data, trait.variable)
    markov.mle <- prop.table(tij, margin = 1)
    q.hat <- get.stationary.distribution(markov.mle)
    h.hat <- get.h.hat(rds.data, trait.variable, network.size)
    est <- as.list((q.hat / h.hat) / sum(q.hat / h.hat))
    
    d = tapply(rds.data[[trait.variable]], rds.data[[trait.variable]], length)
    f = match(rds.data[[trait.variable]], names(est))
    weights <- rep(0, length(f))
    weights[!is.na(f)] = as.numeric(unlist(est[f]))
    weights <- as.numeric(weights / d[f])
    weights[is.na(weights)] <- 0
    SHw <- N * weights / sum(weights)
    if(verbose)
      cat(sprintf(
        "SW prev: %f\n",
        HT.estimate(weights = SHw, outcome = dissample)
      ))
    
    mean.s01 <-
      HT.estimate(weights = SHw[last], outcome = s01s[last])
    if (is.null(crossStat))
      crossStat <- 0.5 * N * mean.s01
    if(verbose){
      print(HT.estimate(weights = SHw[last], outcome = dissample[last]))
      cat(sprintf("SW crossStat: %f\n", 0.5 * N * mean.s01))
    }
    #
    a = tapply(SHw, indices, mean)
    Qstart.all <- sampleTab
    Qstart.all[as.numeric(names(a))] <- 1 / a
  } else{
    #Initialize sampling probabilities as proportional to degree
    #note:  do I want to restrict to be less than 1?
    SSw <-
      gile.ss.weights(
        degsample,
        N = N,
        number.ss.iterations = nit,
        se = FALSE
      )
    a = tapply(SSw, indices, mean)
    Qstart.all <- sampleTab
    Qstart.all[as.numeric(names(a))] <- 1 / a
  }
  
  SSw <-
    gile.ss.weights(
      degsample,
      N = N,
      number.ss.iterations = nit,
      se = FALSE
    )
  a = tapply(SSw, indices, mean)
  
  Qstart.all <- sampleTab
  Qstart.all[as.numeric(names(a))] <- 1 / a
  
  #also, initialize transformation of todis and tonodis according to dummy categories
  todistab <- rep(0, length(dummy))
  tonodistab <- rep(0, length(dummy))
  nonnullindices <- tapply(indices, indices, stats::median)
  todistab[nonnullindices] <- tapply(todisall, indices, sum)
  tonodistab[nonnullindices] <- tapply(tonodisall, indices, sum)
  
  if (is.null(Qstart)) {
    Qstart <- Qstart.all
  }
  
  
  
  istuff <-
    list() # to keep track of what happens in each iteration
  previt <- c() # to keep track of what happens in each iteration
  
  M1small <- ceiling(M1 / 5)
  s01 <-
    (dissample == 1) * tonodisall + (dissample == 0) * todisall
  
  for (itnum in 1:nit) {
    Qstart.all <- Qstart[indices]
    prop.outcome <-
      HT.estimate(weights = Qstart.all, outcome = dissample)
    mean.degree.outcome0 <-
      HT.estimate(weights = Qstart.all, outcome = degsample * (dissample == 0)) / (1 -
                                                                                     prop.outcome)
    mean.degree.outcome1 <-
      HT.estimate(weights = Qstart.all, outcome = degsample * (dissample == 1)) / prop.outcome
    mean.degree <-
      HT.estimate(weights = Qstart.all, outcome = degsample)
    mean.s01 <- HT.estimate(weights = Qstart.all, outcome = s01)
    homophily <-
      2 * prop.outcome * mean.degree.outcome0 * (1 - prop.outcome) * mean.degree.outcome1 / (mean.s01 *
                                                                                               mean.degree)
    if(verbose){
      cat(sprintf("Est. homophily = %f\n", homophily))
      cat(
        sprintf(
          "Est. diff. activity= %f\n",
          mean.degree.outcome1 / mean.degree.outcome0
        )
      )
    }
    
    if (verbose) {
      cat(paste("!!!!Starting on iteration number ", itnum, ".\n", sep = ""))
    }
    istuff[[itnum]] <-
      oneiteration(
        sampleTab,
        indices,
        dummy,
        N,
        Qstart,
        todistab,
        tonodistab,
        dummydis,
        dummydegs,
        M1,
        M2,
        short = short,
        theta0 = theta0,
        parallel = parallel,
        parallel.type = parallel.type,
        MPLEsamplesize = MPLEsamplesize,
        SAN.maxit = SAN.maxit,
        SAN.nsteps = SAN.nsteps,
        interval = interval,
        crossStat = crossStat,
        rawCounts = rawCounts,
        net = net,
        homophily = homophily,
        fixinitial = fixinitial,
        nsamp0 = nsamp0,
        coupons = coupons,
        seed = seed + itnum,
        seed.indices = seed.indices,
        verbose = verbose
      )
    
    Qstart <- istuff[[itnum]]$Q
    theta0 <- (istuff[[itnum]])$coef
    prev <- sum(istuff[[itnum]]$popests * dummydis) / N
    istuff[[itnum]]$prev <- prev
    istuff[[itnum]]$homophily <- homophily
    istuff[[itnum]]$differential.activity <-
      mean.degree.outcome1 / mean.degree.outcome0
    previt <- c(previt, prev)
    
    crossStat <- NULL
    
    if (verbose) {
      cat(sprintf("Prevalence estimate %f\n", prev))
    }
  }
  finalstuff <- istuff[[nit]]
  if (smooth & nit > 1) {
    ns <- min(3, nit)
    ns <- min(2, nit)
    finalstuff$popests <- istuff[[nit]]$popests / ns
    finalstuff$Q       <- istuff[[nit]]$Q / ns
    finalstuff$QQ      <- istuff[[nit]]$QQ / ns
    finalstuff$bsdegsamples <- istuff[[nit]]$bsdegsamples
    finalstuff$bsdissamples <- istuff[[nit]]$bsdissamples
    for (i in ((nit - ns + 1):(nit - 1))) {
      finalstuff$popests <- finalstuff$popests + istuff[[i]]$popests / ns
      finalstuff$Q <-  finalstuff$Q  + istuff[[i]]$Q / ns
      finalstuff$QQ <- finalstuff$QQ + istuff[[i]]$QQ / ns
      finalstuff$bsdegsamples <- rbind(finalstuff$bsdegsamples,
                                       istuff[[i]]$bsdegsamples)
      finalstuff$bsdissamples <- rbind(finalstuff$bsdissamples,
                                       istuff[[i]]$bsdissamples)
      
    }
  }
  
  prev <- sum(finalstuff$popests * dummydis) / N
  
  varest <-
    varEL2(finalstuff$Q, finalstuff$QQ, indices, N, dissample)
  
  varestBS <-
    list(
      BSest = finalstuff$bsprev,
      var = stats::var(finalstuff$bsprev),
      BSmean = finalstuff$bsmean
    )
  
  admindata <-
    list(dummydis = dummydis,
         dummydegs = dummydegs,
         indices = indices)
  return(
      list(
        prev = prev,
        varest = varest$var,
        varest1 = varest$var1,
        coef = (istuff[[nit]])$coef,
        M1 = M1,
        M2 = M2,
        nit = nit,
        theta = theta0,
        nsamp0 = nsamp0,
        coupons = coupons,
        fixinitial = fixinitial,
        MPLEsamplesize = MPLEsamplesize,
        prev.by.it = varest$est,
        weights = (istuff[[nit]])$Q[indices],
        varestBS = varestBS
        
      )
    )
  
}

varEL2 <- function(Q, QQ, indices, N, dissample) {
  # This is the EL estimator based o nthe pairwize inclusion probabilities
  # Q: the inclusion probabilities of the 40 classes
  # QQ: the joint inclusion probabilities of the 40 classes
  # indices: the classes of the 500 samples
  # N: the population size
  # dissample: the disease status of the sample
  # mu: GHT estimate of the prevalence
  Qi <- Q[indices]
  # Qi: vector of inclusion probabilities
  QQi <- QQ[indices,]
  QQi <- QQi[, indices]
  diag(QQi) <- Qi
  # QQi: matrix of joint inclusion probabilities
  #   Compute the EL standard error (if available)
  outcome <- dissample
  wi <- (1 / Qi)
  wi[is.na(wi) | is.infinite(wi)] <- 0
  wij <- (1 / QQi)
  wij[is.na(wij) | is.infinite(wij)] <- 0
  
  iLj <- row(wij) > col(wij)
  oij <- 0.5 * outer(outcome, outcome, "+")
  
  estij <- sum(oij * wij * iLj) / sum(wij * iLj)
  esti  <- sum(outcome * wi) / sum(wi)
  varesti <-
    sum(((outcome - esti) * wi * wi) ^ 2) / (sum(wi * wi) ^ 2)
  varesti <- (N - length(outcome)) * varesti / (N - 1)
  #
  Gi <- apply((wij * wij * (oij - estij)) * iLj, 1, sum)
  wvi <- apply((wij * wij) * iLj, 1, sum)
  varest <- 4 * sum(Gi ^ 2) / (sum(wvi) ^ 2)
  varest <- (N - length(outcome)) * varest / (N - 1)
  list(est = estij, var = varest, var1 = varesti)
}

varGHT <- function(Q, QQ, indices, N, dissample) {
  # This is the unbiased HT estimator
  # Q: the inclusion probabilities of the 40 classes
  # QQ: the joint inclusion probabilities of the 40 classes
  # indices: the classes of the 500 samples
  # N: the population size
  # dissample: the disease status of the sample
  # mu: GHT estimate of the prevalence
  Qi <- Q[indices]
  # Qi: vector of inclusion probabilities
  QQi <- QQ[indices,]
  QQi <- QQi[, indices]
  diag(QQi) <- Qi
  # QQi: matrix of joint inclusion probabilities
  QtQi <- Qi %*% t(Qi)
  Ne = sum(1 / Qi)
  w = 1 / (Qi * Ne)
  mu = sum(dissample * w)
  dev <- dissample - mu
  tmp <- (QQi - QtQi) / QtQi * (dev %*% t(dev)) / QQi
  #   Note that this appears to double count the diagonal!!!
  #   And is not the formula in Thompson or the paper!!!
  varest <- sum((1 - Qi) / (Qi ^ 2) * dev ^ 2) + sum(tmp)
  varest <- varest / (Ne * Ne)
  varest
}
varJGHT <- function(Q, QQ, indices, dissample) {
  # This is the jacknife estimator
  # Q: the inclusion probabilities of the 40 classes
  # QQ: the joint inclusion probabilities of the 40 classes
  # indices: the classes of the 500 samples
  # dissample: the disease status of the sample
  # prev: GHT estimate of the prevalence
  Qi <- Q[indices]
  # Qi: vector of inclusion probabilities
  QQi <- QQ[indices,]
  QQi <- QQi[, indices]
  diag(QQi) <- Qi
  # QQi: matrix of joint inclusion probabilities
  #
  e <- rep(0, length(Qi))
  Ne = sum(1 / Qi)
  w = 1 / (Qi * Ne)
  mu = sum(dissample * w)
  for (j in seq(along = Qi)) {
    Nj = sum(1 / (Qi[-j]))
    wj = 1 / ((Qi[-j]) * Nj)
    muj = sum(dissample[-j] * wj)
    e[j] = (1 - w[j]) * (mu - muj)
  }
  QtQi <- Qi %*% t(Qi)
  tmp <- (QQi - QtQi) * (e %*% t(e)) / QQi
  # Double count the diagonal!!! This is wrong
  # But it seems to work!!! Why!!!
  # And it is *identical* to varGHT!!!
  # when that is double counted too
  varest <- sum(tmp) + sum(diag(tmp))
  varest
}


varBSGHT <- function(Q, dummy, bsdegsamples, bsdissamples) {
  # This is the MS BS estimator
  # Q: the inclusion probabilities of the 40 classes
  # indices: the classes of the 500 samples
  # BSdissamples: the disease status of the BS samples
  # BSdegsamples: the degrees of the BS samples
  e <- rep(0, nrow(bsdegsamples))
  for (j in 1:nrow(bsdegsamples)) {
    indices <-
      toindexref.RDS(bsdegsamples[j,], bsdissamples[j,], dummy)
    #   this is the index of each sampled node within the ordered list of equivalence classes in dummy
    Qi <- Q[indices]
    #   Qi: vector of inclusion probabilities
    Ne = sum(1 / Qi)
    w = 1 / (Qi * Ne)
    e[j] = sum(bsdissamples[j,] * w)
  }
  list(var = stats::var(e),
       BSest = e,
       median = stats::median(e))
}

varNBSGHT <-
  function(degsample,
           dissample,
           bsdegsamples,
           bsdissamples,
           N,
           M1,
           M2,
           popCounts,
           dummy,
           dummydis,
           prev) {
    # This is the BS estimator that scales up the within-network variance
    # Q: the inclusion probabilities of the 40 classes
    # indices: the classes of the 500 samples
    # BSdissamples: the disease status of the BS samples
    # BSdegsamples: the degrees of the BS samples
    if (M1 * M2 != nrow(bsdegsamples)) {
      stop("M1*M2 must equal nrow(bsdegsamples)")
    }
    sampleTab <- tototab.RDS(degsample, dissample, dummy)
    e <- rep(0, M1)
    index <- 0
    for (j in 1:M1) {
      Msamples <- rep(0, length(dummy))
      for (i in 1:M2) {
        index <- index + 1
        #counts of nodes of each type
        tabsamp <-
          tototab.RDS(bsdegsamples[index,], bsdissamples[index,], dummy)
        Msamples <- Msamples + tabsamp
      }
      Q <- (Msamples + 1) / (M2 * popCounts + 1)
      popests <- sampleTab / Q
      popests[sampleTab > 0 & Q == 0] <- 1
      popests[sampleTab == 0 & Q == 0] <- 0
      popests <- N * popests / sum(popests)
      e[j] <- sum(popests * dummydis) / N
    }
    
    list(
      var = stats::var(e),
      mse = mean((e - prev) ^ 2),
      bsprev = mean(e),
      prev = prev
    )
  }

oneiteration <- function(sampleTab,
                         indices,
                         dummy,
                         N,
                         Qstart,
                         todistab,
                         tonodistab,
                         dummydis,
                         dummydegs,
                         M1,
                         M2,
                         short = FALSE,
                         theta0 = -0.79,
                         SAN.maxit = 3,
                         SAN.nsteps = 2^19,
                         interval = 10000,
                         parallel = 0,
                         parallel.type = "PSOCK",
                         crossStat = NULL,
                         rawCounts = NULL,
                         net = NULL,
                         homophily = NULL,
                         MPLEsamplesize = 50000,
                         fixinitial = 1,
                         nsamp0 = 10,
                         coupons = 2,
                         seed = 1,
                         seed.indices = NULL,
                         verbose = TRUE) {
  maxdeg <- length(dummy) / 2 - 1
  prev <-
    sum((sampleTab / Qstart)[dummydis == 1 &
                               Qstart > 0], na.rm = T) / sum((sampleTab / Qstart)[Qstart > 0], na.rm =
                                                               T)
  #
  crossStat <- NULL
  #
  if (is.null(crossStat)) {
    crossStat <-
      getCrossStat(
        sampleTab,
        todistab,
        tonodistab,
        dummydis,
        dummydegs,
        Qstart,
        N = N,
        homophily = homophily
      )
  }
  newPop <-
    getPopulation(
      sampleTab = sampleTab,
      Qstart = Qstart,
      N = N,
      dummydegs = dummydegs,
      dummydis = dummydis,
      homophily = homophily,
      todistab = todistab,
      tonodistab = tonodistab,
      rawCounts = rawCounts,
      verbose=verbose
    )
  if (verbose) {
    cat(sprintf(
      "newPop$crossStat: %f, crossStat: %f\n",
      newPop$crossStat,
      crossStat
    ))
  }
  if (verbose) {
    a <-
      ergm::summary_formula(newPop$initialNet ~ nodefactor('disease', base =
                                                             0))
    cat(
      sprintf(
        "ergm::summary_formula(newPop$initialNet~nodefactor('disease',base=0)): 0=%f 1=%f\n",
        a[1],
        a[2]
      )
    )
  }
  n <- length(indices)
  if (is.null(net)) {
    fit <- fitnet.RDS(
      target.stats = crossStat,
      theta0 = theta0,
      initialNet = newPop$initialNet,
      short = short,
      parallel = parallel,
      parallel.type = parallel.type,
      MPLEsamplesize = MPLEsamplesize,
      SAN.maxit = SAN.maxit,
      SAN.nsteps = SAN.nsteps,
      maxdeg = maxdeg,
      verbose = verbose
    )
    initialNet <- fit$newnetwork
    if (verbose) {
      cat(
        sprintf(
          "fitted differential.activity= %f\n",
          network.differential.activity(initialNet, 'disease')
        )
      )
      cat(sprintf(
        "fitted homophily= %f\n",
        network.homophily(initialNet, 'disease')
      ))
      cat(
        sprintf(
          "ergm::summary_formula(fit~nodemix('disease'))= %f\n",
          ergm::summary_formula(initialNet ~ nodemix("disease", levels2 =
                                                       c(1, 3)))
        )
      )
      a = ergm::summary_formula(initialNet ~ nodefactor('disease', base =
                                                          0))
      cat(
        sprintf(
          "ergm::summary_formula(initialNet~nodefactor('disease',base=0)): 0=%f 1=%f\n",
          a[1],
          a[2]
        )
      )
    }
  }
  if (!short) {
    if (parallel > 1) {
      sequential <- FALSE
    } else{
      sequential <- TRUE
    }
    if (is.null(net)) {
      sfn <- function(i, fit, interval) {
        as_edgelist_compressed_rds(fit$newnetwork)
      }
    } else{
      fit <- net
      sfn <- function(i, fit, interval) {
        as_edgelist_compressed_rds(fit)
      }
    }
    if (F) {
      # Standard version (based on simulations)
      if (parallel > 1) {
        #      Run the jobs with rpvm or Rmpi
        cl <-
          beginparallel(parallel = parallel,
                        type = parallel.type,
                        seed = seed)
        if (verbose) {
          cat(paste("Starting parallel network simulations\n"))
        }
        sim <-
          parallel::clusterApplyLB(cl, as.list(1:M1), sfn, fit, interval)
        if (verbose) {
          cat(paste("Finished parallel network simulations\n"))
        }
        endparallel(cl, type = parallel.type)
      } else{
        sim <- vector(mode = "list", length = M1)
        if (verbose) {
          cat(paste("Starting non-parallel network simulations\n"))
        }
        for (j in 1:M1) {
          sim[[j]] <- sfn(j, fit, interval)
          initialNet = as_network_uncompressed(sim[[j]])
          if (verbose) {
            cat(
              sprintf(
                "ergm::summary_formula(sim~nodemix('disease'))= %f\n",
                ergm::summary_formula(initialNet ~ nodemix("disease", levels2 =
                                                             c(1, 3)))
              )
            )
            cat(
              sprintf(
                "simulated differential.activity = %f\n",
                network.differential.activity(initialNet, "disease")
              )
            )
            cat(sprintf(
              "simulated homophily = %f\n",
              network.homophily(initialNet, "disease")
            ))
          }
        }
        if (verbose) {
          cat(paste("Finished non-parallel network simulations\n"))
        }
      }
    } else{
      # SAN version
      if (verbose) {
        cat(paste("Using SAN network (rather than simulations)\n"))
      }
      sim <- vector(1, mode = "list")
      sim[[1]] <- sfn(1, fit, interval)
    }
  } else{
    #
    sim <- stats::simulate(
      fit,
      nsim = M1,
      control = ergm::control.simulate.ergm(
        seed = seed,
        MCMC.burnin = 10000,
        MCMC.interval = 10000
      ),
      verbose = verbose
    )
  }
  if (verbose) {
    cat("Starting RDS using C code...\n")
  }
  samps <-
    getsamples.RDS.C(
      sim,
      dummy,
      M1,
      M2,
      N,
      nsamp0,
      fixinitial,
      n,
      coupons,
      dummydis = dummydis,
      parallel = parallel,
      parallel.type = parallel.type,
      seed = seed,
      seed.indices = seed.indices,
      verbose=verbose
    )
  
  net <- as_network_uncompressed_rds(sim[[1]])
  deg <- sapply(net$iel, length) + sapply(net$oel, length)
  deg[deg > maxdeg] <- maxdeg
  disease <- net %v% "disease"
  disease <- 1 * (disease == max(disease))
  popCounts <- tototab.RDS(deg, disease, dummy)
  if (verbose) {
    cat("Ended RDS...\n")
  }
  M12 <- M1 * M2
  Q <- (samps$Msamples + 1) / (M12 * popCounts + 1)
  QQ <-
    (samps$bisamples + (sampleTab %*% t(sampleTab)) / N) / (M12 * (popCounts %*% t(popCounts)) + (popCounts %*% t(popCounts)) /
                                                              N)
  Q[is.na(Q)] <- 0
  QQ[is.na(QQ)] <- 0
  diag(QQ) <- Q
  
  popests <- sampleTab / Q
  popests[sampleTab > 0 & Q == 0] <- 1
  popests[sampleTab == 0 & Q == 0] <- 0
  popests <- N * popests / sum(popests)
  list(
    Q = Q,
    coef = fit$coef,
    crossStat = crossStat,
    popests = popests,
    QQ = QQ,
    bisamples = samps$bisamples,
    popCounts = newPop$rawpopCounts,
    sampbynet = samps$sampbynet,
    bsprev = samps$bsprev,
    bsmean = mean(disease),
    bsQ = samps$bsQ,
    bsdegsamples = samps$bsdegsamples,
    bsdissamples = samps$bsdissamples
  )
}


fitnet.RDS <- function(target.stats,
                       theta0,
                       initialNet,
                       short = FALSE,
                       parallel = 0,
                       parallel.type = "PSOCK",
                       MPLEsamplesize = 50000,
                       SAN.maxit = 3,
                       SAN.nsteps = 2^19,
                       deltadis = 0.05,
                       maxdeg = 30,
                       verbose = TRUE) {
  if (!short) {
    # Generate a sample network using MCMC
    srun <- 0
    obs <- target.stats - target.stats
    #
    # I wish I knew why this was necessary, but it is and I can't seem to fix it.
    #
    formula <- initialNet ~ nodemix("disease", levels2 = c(1, 3))
    #
    formula <- statnet.common::nonsimp_update.formula(formula, initialNet ~ .)
    dform <- statnet.common::nonsimp_update.formula(formula,
                                       stats::as.formula(
                                         paste('. ~ . + degree(0:', maxdeg + 2, ',by="disease")', sep = "")
                                       ))
    dis <- network::get.vertex.attribute(initialNet, "disease")
    if (mean(dis == 0) > 0.99 |
        mean(dis == 1) > 0.99 | target.stats < 1) {
      cat(
        sprintf(
          "Network too extreme... mean(dis==0)= %f; target.stats= %f\n",
          mean(dis == 0),
          target.stats
        )
      )
      fit <-
        structure(
          list(
            coef = -5,
            newnetwork = initialNet,
            formula = formula,
            constraints =  ~ degrees,
            control = NULL
          ),
          class = "ergm"
        )
      return(fit)
    }
    obs <- ergm::summary_formula(formula)
    obs.prev <- obs + 1
    obs.prev2 <- obs.prev + 1
    if (verbose) {
      cat(sprintf(
        "number of diseased= %d; prevalence= %f\n",
        sum(network::get.vertex.attribute(initialNet, "disease")),
        mean(network::get.vertex.attribute(initialNet, "disease"))
      ))
      cat(
        sprintf(
          "ergm::summary_formula(initialNet~nodemix('disease'))= %f\n",
          ergm::summary_formula(initialNet ~ nodemix("disease", levels2 =
                                                       c(1, 3)))
        )
      )
    }
    while (isTRUE(all.equal(obs, obs.prev2)) &
           (srun < SAN.maxit) &
           (sum((obs - target.stats) ^ 2, na.rm = TRUE) > 5)) {
      dd <-
        ergm::summary_formula(stats::as.formula(
          paste(
            'initialNet ~ nodemix("disease", levels2 = c(1, 3)) + degree(0:',
            maxdeg + 2,
            ',by="disease")',
            sep = ""
          )
        ))
      dd[1] <- target.stats
      initialNet <- ergm::san(
        object=dform,
        target.stats = dd,
        control = ergm::control.san(
          SAN.init.maxedges = ergm::summary_formula(initialNet ~ edges) + 100,
          SAN.nsteps = SAN.nsteps,
          parallel = 0
        ),
        verbose = verbose,
        nsim = 1,
        constraints =  ~ edges
      )
      if (verbose) {
        cat(
          sprintf(
            "ergm::summary_formula(initialNet~nodemix('disease'))= %f\n",
            ergm::summary_formula(initialNet ~ nodemix("disease", levels2 =
                                                         c(1, 3)))
          )
        )
      }
      obs.prev2 <- obs.prev
      obs.prev <- obs
      obs <- ergm::summary_formula(formula)
      obs[is.na(obs)] <- 0
      srun <- srun + 1
      if (verbose) {
        cat(paste("Finished SAN run", srun, "\n"))
      }
      if (verbose) {
        cat(sprintf("SAN summary statistics: %f\n", obs))
        cat(sprintf("Meanstats Goal: %f\n", target.stats))
        cat(sprintf(
          "Difference: SAN - target.stats = %f\n",
          round(obs - target.stats, 0)
        ))
      }
      dis <- network::get.vertex.attribute(initialNet, "disease")
      if (all(obs / target.stats < 0.75) & mean(dis == 1) < 0.99) {
        if (verbose) {
          cat("Increasing the number of diseased (and hence the number of cross-ties)\n")
        }
        mindelta <- min(floor(sum(dis == 0) * deltadis),
                        floor(sum(dis == 1) * deltadis))
        a <- sample(seq(along = dis)[dis == 0], size = mindelta)
        cat(sprintf(
          "Number of diseased: %d; Proportion %f\n",
          sum(dis == 1),
          mean(dis)
        ))
        dis[a] <- 1
        if (mean(dis == 1) < 0.99) {
          cat(sprintf(
            "New number of diseased: %d; Proportion %f\n",
            sum(dis == 1),
            mean(dis)
          ))
          initialNet <-
            network::set.vertex.attribute(initialNet, "disease", dis)
        }
      }
      if (all(obs / target.stats > 1 / 0.75) &
          mean(dis == 1) > 0.01) {
        if (verbose) {
          cat("Decreasing the number of diseased (and hence the number of cross-ties)\n")
        }
        mindelta <- min(floor(sum(dis == 0) * deltadis),
                        floor(sum(dis == 1) * deltadis))
        a <- sample(seq(along = dis)[dis == 1], size = mindelta)
        cat(sprintf(
          "Number of diseased: %d; Proportion %f\n",
          sum(dis == 1),
          mean(dis)
        ))
        dis[a] <- 0
        if (mean(dis == 1) > 0.01) {
          cat(sprintf(
            "New number of diseased: %d; Proportion %f\n",
            sum(dis == 1),
            mean(dis)
          ))
          initialNet <-
            network::set.vertex.attribute(initialNet, "disease", dis)
        }
      }
    }
    if (verbose) {
      cat("Starting MPLE...\n")
    }
    fit <-
      structure(
        list(
          coef = -5,
          newnetwork = initialNet,
          formula = formula,
          constraints =  ~ degrees,
          control = NULL
        ),
        class = "ergm"
      )
    fit
  } else{
    fit <- ergm::ergm(
      formula,
      target.stats = target.stats,
      constraints =  ~ degrees,
      verbose = verbose,
      eval.loglik = FALSE,
      estimate = "MPLE",
      control = ergm::control.ergm(
        MCMC.samplesize = 10000,
        loglik.control = ergm::control.logLik.ergm(warn.dyads =
                                                     FALSE)
      )
    )
    fit$coef[is.na(fit$coef)] <- 0
    cat(paste("Final theta: ", round(fit$coef, 4), "\n", sep = ""))
    fit
  }
  
}



getPopulation <- function(sampleTab,
                          Qstart,
                          N,
                          dummydegs,
                          dummydis,
                          homophily = NULL,
                          todistab = NULL,
                          tonodistab = NULL,
                          rawCounts = NULL,
                          verbose = TRUE) {
  if (is.null(rawCounts)) {
    aaa <-
      sampleTab / Qstart # this is the estimate of the population counts
    aaa[sampleTab == 0 & Qstart == 0] <- 0
    aaa[sampleTab > 0 &
          Qstart == 0] <-
      0  # so far, not integers.  need to make it so.
    #   VIP uncomment this to have a normalized networked population size at N
    #   rather than an estimated one
    rawCounts <- N * aaa / sum(aaa)
    rawCounts[rawCounts < sampleTab] <-
      sampleTab[rawCounts < sampleTab]
  }
  if (is.null(homophily)) {
    homophily = 1
  }
  if (verbose) {
    cat(sprintf("Homophily estimated to be %f\n", homophily))
  }
  if (homophily < 5) {
    done <- FALSE
    while (!done) {
      #want to use reedmolloy to find a suitable network.  iterate finding of counts at the same time.
      Nhat <- try(roundstoc.RDS(rawCounts))
      newdegs <- rep(dummydegs, times = Nhat)
      newdis <- rep(dummydis, times = Nhat)
      if (sum(newdegs) == 2 * round(sum(newdegs) / 2)) {
        simr <- try(reedmolloy(newdegs[newdegs > 0]), silent = TRUE)
        if (!inherits(simr, "try-error")) {
          done <- TRUE
        }
      }
    }
    if (any(newdegs == 0)) {
      add.vertices(simr, sum(newdegs == 0))
      newdis <- c(newdis[newdegs > 0], newdis[newdegs == 0])
      simr %v% "vertex.names" <- paste(1:network.size(simr))
    }
    # Assign the disease states. Note that reedmolloy reverses the degree
    # sequence, so we need to reverse the disease states too
    # simr<- network::set.vertex.attribute(simr,"disease",rev(newdis))
    simr <- network::set.vertex.attribute(simr, "disease", newdis)
  } else{
    # So try the high homophily case
    done <- FALSE
    while (!done) {
      #want to use reedmolloy to find a suitable network.  iterate finding of counts at the same time.
      Nhat <- roundstoc.RDS(rawCounts)
      newdegs <- rep(dummydegs, times = Nhat)
      newdis <- rep(dummydis, times = Nhat)
      #  sel is the non diseased
      sel <- (newdis == 0)
      if (verbose) {
        cat(
          sprintf(
            "Working network number of non-diseased %d, and diseased %d\n",
            sum(newdis == 0),
            sum(newdis == 1)
          )
        )
      }
      if (sum((newdis == 1) &
              (newdegs == 1)) > sum((newdis == 1) &
                                    (newdegs > 1))) {
        #
        sele <- seq(along = newdis)[(newdis == 1) & (newdegs == 1)]
        # sele are those with
        while (inherits(try(igraph::degree.sequence.game(newdegs[!sel], method =
                                                         "simple"))
                        , "try-error")) {
          sel[sample(sele, size = ceiling((sum((newdis == 1)
          ) - sum((newdis == 1) & (newdegs > 1)
          )) * 0.1), replace = FALSE)] <- TRUE
        }
      }
      if (sum(newdegs[sel]) == 2 * round(sum(newdegs[sel]) / 2)
          &
          sum(newdegs[!sel]) == 2 * round(sum(newdegs[!sel]) / 2)) {
        done <- TRUE
        simr <-
          network.initialize(length(newdegs), directed = FALSE)
        if (any(sel)) {
          # Generate ties between the non-diseased
          sm0 <-
            try(igraph::get.edgelist(igraph::degree.sequence.game(newdegs[sel], method =
                                                                    "vl")))
          if (inherits(sm0, "try-error")) {
            done <- FALSE
          } else{
            simr <- network::add.edges(x = simr,
                                       tail = as.list(sm0[, 1]),
                                       head = as.list(sm0[, 2]))
            if (any(!sel)) {
              sm1 <-
                try(igraph::get.edgelist(igraph::simplify(
                  igraph::degree.sequence.game(newdegs[!sel], method = "simple")
                )))
              if (inherits(sm1, "try-error")) {
                done <- FALSE
              } else {
                simr <- network::add.edges(
                  x = simr,
                  tail = as.list(length(newdegs[sel]) +
                                   sm1[, 1]),
                  head = as.list(length(newdegs[sel]) +
                                   sm1[, 2])
                )
              }
            }
          }
        }
        simr <-
          network::set.vertex.attribute(simr, "disease", sort(newdis))
      }
    }
  }
  # This is the new way. MSH is unsure which is correct
  sampleTab[sampleTab == 0] <- NA
  crossStatAll = 0.5 * sum(todistab[dummydis == 0] * Nhat[dummydis == 0] /
                             sampleTab[dummydis == 0], na.rm = TRUE) + 0.5 * sum(tonodistab[dummydis ==
                                                                                              1] * Nhat[dummydis == 1] / sampleTab[dummydis == 1], na.rm = TRUE)
  #
  # This is the new new way. MSH is unsure which is better
  #
  prop.outcome <- sum(Nhat[dummydis == 1]) / sum(Nhat)
  mean.degree.outcome0 <-
    sum((Nhat * dummydegs)[dummydis == 0]) / sum(Nhat[dummydis == 0])
  mean.degree.outcome1 <-
    sum((Nhat * dummydegs)[dummydis == 1]) / sum(Nhat[dummydis == 1])
  mean.degree <- sum((Nhat * dummydegs)) / sum(Nhat)
  mean.s01 <-
    2 * prop.outcome * mean.degree.outcome0 * (1 - prop.outcome) * mean.degree.outcome1 / (homophily *
                                                                                             mean.degree)
  crossStat <- 0.5 * sum(Nhat) * mean.s01
  list(
    popCounts = Nhat,
    rawpopCounts = rawCounts,
    initialNet = simr,
    newDis = newdis,
    newDegs = newdegs,
    crossStat = crossStatAll
  )
}



getCrossStat <-
  function(sampleTab,
           todistab,
           tonodistab,
           dummydis,
           dummydegs,
           Qstart,
           N = sum(1 / Qstart),
           homophily) {
    # mtodis is vector of number of ties to infected for each sampled node
    # mtonodis is vector of number of ties to uninfected for each sampled node
    # dissample is vector of observed diseases
    # Q is computed sampling probability.
    aaa <-
      sampleTab / Qstart # this is the estimate of the population counts
    aaa[sampleTab == 0 & Qstart == 0] <- 0
    aaa[sampleTab > 0 &
          Qstart == 0] <-
      0  # so far, not integers.  need to make it so.
    Nhat <- N * aaa / sum(aaa)
    Nhat[Nhat < sampleTab] <- sampleTab[Nhat < sampleTab]
    #
    Qstart[Qstart == 0] <- NA
    
    cs <-
      0.5 * sum(todistab[dummydis == 0] / Qstart[dummydis == 0], na.rm = TRUE) + 0.5 *
      sum(tonodistab[dummydis == 1] / Qstart[dummydis == 1], na.rm = TRUE)
    # This is the new new way. MSH is unsure which is better
    #
    prop.outcome <- sum(Nhat[dummydis == 1]) / sum(Nhat)
    mean.degree.outcome0 <-
      sum((Nhat * dummydegs)[dummydis == 0]) / sum(Nhat[dummydis == 0])
    mean.degree.outcome1 <-
      sum((Nhat * dummydegs)[dummydis == 1]) / sum(Nhat[dummydis == 1])
    mean.degree <- sum((Nhat * dummydegs)) / sum(Nhat)
    mean.s01 <-
      2 * prop.outcome * mean.degree.outcome0 * (1 - prop.outcome) * mean.degree.outcome1 / (homophily *
                                                                                               mean.degree)
    cs
  }


tototab.RDS <- function(deggs, diss, dummy) {
  index <- deggs * 10 + diss
  totabulate.RDS(index, dummy)
}

totabulate.RDS <- function(index, dummy) {
  aaa <-
    tabulate(c(index, dummy) + 1) # intermediate variable to get to counts
  tab <- aaa[aaa > 0] - 1 #
  tab
}

roundstoc.RDS <- function(vec) {
  # takes a vector and makes it integers, keeping the total the same
  target <- sum(vec)
  temp <- floor(vec)
  needed <- target - sum(temp)
  away <- vec - temp
  while (needed > .5) {
    toget <- sample(c(1:length(vec)), size = 1, prob = away)
    temp[toget] <- temp[toget] + 1
    away <- vec - temp
    away[away < 0] <- 0
    needed <- needed - 1
  }
  temp
}

##  Need a function to translate the sample to indices of dummy
toindexref.RDS <- function(degsample, dissample, dummy) {
  aaa <- rep(0, length(degsample))
  temp <- degsample * 10 + dissample
  for (i in 1:length(degsample)) {
    aaa[i] <- which(dummy == temp[i])
  }
  aaa
}

# krista.to.rds.data.frame.RDS <- function(x, population.size = NULL) {
#   source.data.frame <- data.frame(
#     wave = x$wsample,
#     degree = x$degsample,
#     disease = x$dissample,
#     todiseased = x$todis,
#     tonondiseased = x$tonodis,
#     id = x$nsample,
#     recruiter.id = x$nominators
#   )
#   as.rds.data.frame(source.data.frame,
#     population.size = population.size,
#     network.size = "degree"
#   )
# }


#
# Find a graph from a given degree sequence
#
reedmolloy <- function(deg,
                       maxit = 10,
                       verbose = TRUE) {
  sm <-
    .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg, method =
                                                                     "vl")))
  iter <- 0
  if (!is.null(sm$error)) {
    sm <-
      .catchToList(igraph::get.edgelist(
        igraph::degree.sequence.game(deg, method = "simple.no.multiple")
      ))
    while (!is.null(sm$error) & iter < maxit) {
      jitterdeg <- sample(seq_along(deg), size = 2, prob = deg)
      deg[jitterdeg] <- deg[jitterdeg] + 2 * (runif(2) > 0.5) - 1
      deg[deg == 0] <- 2
      sm <-
        .catchToList(igraph::get.edgelist(
          igraph::degree.sequence.game(deg, method = "simple.no.multiple")
        ))
      iter <- iter + 1
    }
  }
  if (iter >= maxit) {
    stop(
      'The reedmolloy function failed to form a valid network from the passed degree sequence.'
    )
  }
  smn <- network::network.initialize(length(deg), directed = FALSE)
  smn <- network::add.edges.network(x = smn,
                                    tail = as.list(sm$value[, 1]),
                                    head = as.list(sm$value[, 2]))
  return(smn)
}


network.differential.activity <-
  function(y, group.variable = "disease") {
    y <- as_network_uncompressed_rds(y)
    gv <- network::get.vertex.attribute(y, group.variable)
    u <- sort(unique(gv))
    gv <- match(gv, u)
    y <- network::set.vertex.attribute(y, group.variable, paste(gv))
    mm <- sprintf("y~nodefactor('%s',base=0)", group.variable)
    mm <- ergm::summary_formula(stats::as.formula(mm)) / table(paste(gv))
    da <- mm[2] / mm[1]
    da[is.infinite(da) | table(paste(gv))[1] == 0] <- NA
    names(da) <- "differential activity"
    da
  }

#Force the input into edgelist form.  Network size, directedness, and vertex
#names are stored as attributes, since they cannot otherwise be included
#A copy of as.edgelist.sna
as_edgelist_compressed_rds <-
  function(x,
           attrname = NULL,
           force.bipartite = FALSE) {
    #In case of lists, process independently
    if (is.list(x) && (!(is(x,"network"))))
      return(
        lapply(
          x,
          as_edgelist_compressed_rds,
          attrname = attrname,
          force.bipartite = force.bipartite
        )
      )
    #Begin with network objects
    if (is(x,"network")) {
      out <- network::as.matrix.network.edgelist(x, attrname = attrname)
      if (NCOL(out) == 2)
        #If needed, add edge values
        out <- cbind(out, rep(1, NROW(out)))
      attr(out, "n") <- network::network.size(x)
      attr(out, "directed") <- network::is.directed(x)
      attr(out, "vnames") <- network::network.vertex.names(x)
      van <- network::list.vertex.attributes(x)
      if (length(van) > 0) {
        va <- vector(mode = "list", length(van))
        for (i in (1:length(van))) {
          va[[i]] <- network::get.vertex.attribute(x, van[i], unlist = TRUE)
        }
        names(va) <- van
        attr(out, "vertex.attributes") <- va
      }
      if (is.bipartite(x))
        attr(out, "bipartite") <-
        network::get.network.attribute(x, "bipartite")
      else if (force.bipartite)
        out <-
        as_edgelist_compressed_rds(out, attrname = attrname, force.bipartite = force.bipartite)
    } else{
      warning(
        "as_edgelist_compressed_rds input must be network, or list thereof.\n Returning the original object.\n"
      )
      return(x)
    }
    #Return the result
    out
  }

###############################################################################
# The <as_network_uncompressed> function is basically the inverse of the above
# <as_edgelist_compressed> function
#
# --PARAMETERS--
#   x         : a compressed network or a network
#   edge.check: whether computationally expensive checks of the legality
#               of submitted edges should be performed (T or F); default=FALSE
#
# --IGNORED PARAMTERS--
#   na.rm:  whether NA valuse should be removed for ??; default=FALSE
#   ...  :  additional parameters for flexibility
#
# --RETURNED--
#   x: the original network if it is already uncompressed or if 'x' is neither
#      a compressed or uncompressed network
#   g: the uncompressed version of x
#
###############################################################################
as_network_uncompressed <- function(x,
                                    na.rm = FALSE,
                                    edge.check = FALSE,
                                    ...) {
  #Initialize the network object
  if (inherits(x, "network")) {
    return(x)
  }
  if (is.null(attr(x, "vnames"))) {
    warning(
      "as_network_uncompressed input must be a compressed network, or a network.\n Returning the original object.\n"
    )
    return(x)
  }
  n <- attr(x, "n")
  directed <- attr(x, "directed")
  g <- network::network.initialize(n, directed = directed)
  #Call the specific coercion routine, depending on matrix type
  g <-
    network::add.edges(g, as.list(x[, 1]), as.list(x[, 2]), edge.check = edge.check)
  va <- attr(x, "vertex.attributes")
  if (length(va) > 0) {
    for (i in (1:length(va))) {
      g <- network::set.vertex.attribute(g, names(va)[i], va[[i]])
    }
  }
  #Return the result
  g
}

as_network_uncompressed_rds <- function(x,
                                        na.rm = FALSE,
                                        edge.check = FALSE) {
  #Initialize the network object
  if (is(x,"network")) {
    return(x)
  }
  if (is.null(attr(x, "vnames"))) {
    warning(
      "as_network_uncompressed_rds input must be a compressed network, or a network.\n Returning the original object.\n"
    )
    return(x)
  }
  n <- attr(x, "n")
  directed <- attr(x, "directed")
  g <- network::network.initialize(n, directed = directed)
  #Call the specific coercion routine, depending on matrix type
  # g<-network.edgelist(x,g,na.rm=na.rm,edge.check=edge.check)
  g <-
    network::add.edges(g, as.list(x[, 1]), as.list(x[, 2]), edge.check = edge.check)
  va <- attr(x, "vertex.attributes")
  if (length(va) > 0) {
    for (i in (1:length(va))) {
      g <- network::set.vertex.attribute(g, names(va)[i], va[[i]])
    }
  }
  #Return the result
  g
}


network.homophily <- function(net,
                              group.variable = "disease",
                              return.finite = TRUE,
                              alpha = FALSE) {
  net <- as_network_uncompressed_rds(net)
  gv <- network::get.vertex.attribute(net, group.variable)
  net <-
    network::set.vertex.attribute(net, group.variable, paste(gv))
  mf <- ergm::summary_formula(stats::as.formula(sprintf(
    "net~nodefactor('%s',base=0)", group.variable
  )))
  mm = sprintf("net~nodemix('%s')", group.variable)
  mm <- diag(mf)
  mm[col(mm) >= row(mm)] <-
    ergm::summary_formula(stats::as.formula(sprintf("net~nodemix('%s')", group.variable)))
  mt = mm + t(mm)
  mm = mm + t(mm) - diag(diag(mm))
  ff <- apply(mt, 1, sum)
  alpha = sum(ff) * mt / outer(ff, ff, "*")
  a = outer(mf, mf, "*") / outer(mf, mf, "+")
  h = sum(a[col(a) < row(a)]) / sum(mm[col(mm) < row(mm)])
  if (return.finite &
      (is.infinite(h) | is.nan(h) | is.na(h))) {
    h = sum(a[col(a) < row(a)])
  }
  names(h) = "homophily"
  dh <-
    mf * (sum(mf) - mf) / (sum(mf) * (apply(mm, 2, sum) - diag(mm)))
  bad.dh <- is.infinite(dh) | is.nan(dh) | is.na(dh)
  if (return.finite & any(bad.dh)) {
    dh[bad.dh] <- (mf * (sum(mf) - mf) / sum(mf))[bad.dh]
  }
  names(dh) <-
    paste(sort(unique(
      network::get.vertex.attribute(net, attrname = group.variable)
    )))
  attr(h, "differential") <- dh
  attr(h, "alpha") <- alpha
  # alpha
  class(h) = "network.homophily"
  h
}
