#' Create RDS samples with given characteristics
#' 
#' @param net the network object from which to draw a sample
#' @param nnodes the number of nodes in the network [at least as default]
#' @param nsamp0 the number of seeds to be drawn (i.e. the size of the 0th wave of sampling)
#' @param fixinitial a variable that indicates the distribution from which to draw the initial seeds, if the seeds
#' variable is NULL and the seed.distribution variable is NULL
#' @param nsamp number of individuals in each RDS sample
#' @param replace sampling with replacement
#' @param coupons number of coupons
#' @param select not used
#' @param bias not used
#' @param rds.samp not used
#' @param seed.distribution a variable [what kind?] that indicates the distribution from which to draw the initial seeds
#' @param attrall Whether all the information about the sample should be returned [??]
#' @param trait.variable attribute of interest
#' @param nsims number of RDS samples to draw
#' @param seeds an array of seeds. Default is NULL, in which case the function draws the seeds from the nodes of the network.
#' @param prob.network.recall simulates the probability that an individual will remember any particular link
#' @param verbose Print verbose output
#' @return A list with the following elements:
#' nsample: vector of indices of sampled nodes
#' wsample: vector of waves of each sampled node
#' degsample: vector of degrees of sampled nodes
#' attrsample: vector of attrs of sampled nodes
#' toattr: vector of numbers of referrals to attrsd nodes
#' tonoattr: vector of number of referrans to unattrsd
#' nominators: recruiter of each sample
#' @export
rdssampleC <- function(net,
                       nnodes = network.size(net),
                       nsamp0,
                       fixinitial,
                       nsamp,
                       replace,
                       coupons,
                       select = NULL,
                       bias = NULL,
                       rds.samp = NULL,
                       seed.distribution = NULL,
                       attrall = FALSE,
                       trait.variable = "disease",
                       nsims = 1,
                       seeds = NULL,
                       prob.network.recall = 1,
                       verbose = TRUE) {
  #via Coupons to the C function and minor changes there.
  net <- as_network_uncompressed_rds(net)
  attrs <- network::get.vertex.attribute(net, trait.variable)
  if (is.null(attrs) || all(is.na(attrs))) {
    stop(sprintf("No variable called %s appears in the data.", trait.variable))
  }
  degs <- sapply(net$iel, length) + sapply(net$oel, length)
  
  if (is.logical(fixinitial[1]) &
      fixinitial[1] == FALSE) {
    fixinitial <- -1
  }
  
  #next piece of code seems to be modifying seed distribution
  if (is.null(seed.distribution)) {
    seed.distribution <- switch(
      as.character(fixinitial[1]),
      "-3" = {
        rep(1, length = nnodes)
      },
      "-2" = {
        #        So attr(fixinitial,"table.seeds")=table(outcome[wave==0])
        if (any(attrs == 1) &
            any(attrs == 0)) {
          table.seeds = attr(fixinitial, "table.seeds")
          attrs * table.seeds[2] / sum(attrs ==
                                         1) + (1 - attrs) * table.seeds[1] / sum(attrs == 0)
        } else{
          rep(1, length = nnodes)
        }
      },
      "-1" = {
        degs
      },
      "1" = {
        attrs
      },
      "2" = {
        attrs * degs
      },
      "3" = {
        degs
      },
      degs
    )
  }
  
  # Select seeds according to the seed distribution
  if (is.null(seeds)) {
    my.seeds <- c()
    for (i in 1:nsims) {
      temp.seeds <-
        sample(
          x = seq(along = seed.distribution),
          size = nsamp0,
          replace = FALSE,
          prob = seed.distribution
        )
      my.seeds <- c(my.seeds, temp.seeds)
    }
  } else{
    my.seeds <- seeds
  }
  if (sum(seed.distribution > 0) < nsamp0) {
    if (verbose == TRUE) {
      print(paste(
        "No individuals available with positive probability of",
        "being a seed."
      ))
    }
  }
  
  n <- network.size(net)
  
  target_num_recruits <- rep(coupons, n)
  
  #code i wrote integrated into new function minorly modified.
  Clist <- list(
    n = n,
    dir = FALSE,
    bipartite = FALSE,
    ndyads = n * (n - 1) / 2
  )
  
  e <-
    as.edgelist(net, attrname = NULL) # Ensures that for undirected networks, tail<head
  if (length(e) == 0) {
    stop("The network is empty!")
  }
  if (!is.matrix(e)) {
    e <- matrix(e, ncol = 2)
  }
  
  Clist$nedges <- dim(e)[1]
  Clist$tails <- e[, 1]
  Clist$heads <- e[, 2]
  
  my.recruitedSampleAll <- c(rep(0, nsims * nsamp * 10))
  my.recruitersAll <-  c(rep(0, nsims * nsamp * 10))
  my.recruitTimesAll <- c(rep(0, nsims * nsamp * 10))
  
  my.recruitedSample <- c(rep(0, nsamp * 10))
  my.recruiters <- c(rep(0, nsamp * 10))
  my.recruitTimes <- c(rep(0, nsamp * 10))
  max.log <- 15000
  my.log <-
    c(rep("                                                          ", max.log))
  my.seeds0 <-
    c(rep(0, nsamp0)) #just added this (Maggie), also changed call to .C
  log.length <- 0
  #call to the C function
  
   result <- .C(
    "CRDSSample",
    as.integer(Clist$tails),
    as.integer(Clist$heads),
    as.integer(Clist$nedges),
    as.integer(Clist$n),
    as.integer(Clist$dir),
    as.integer(Clist$bipartite),
    as.double(seed.distribution),
    as.integer(my.recruitedSample),
    as.integer(my.recruitedSampleAll),
    as.integer(my.recruiters),
    as.integer(my.recruitersAll),
    as.double(my.recruitTimes),
    as.double(my.recruitTimesAll),
    as.character(my.log),
    as.integer(target_num_recruits),
    as.integer(my.seeds0),
    as.integer(my.seeds),
    as.integer(nsamp0),
    as.integer(nsamp),
    as.integer(nsims),
    as.integer(verbose),
    as.integer(max.log),
    as.integer(log.length),
    PACKAGE = "RDS"
   )
  init.deg <- degs[result[[9]]]
  #simulates imperfect recall using the binomial distribution
  if (prob.network.recall == 1) {
    degree <- init.deg
  } else {
    degree <- c()
    for (i in 1:length(init.deg)) {
      degree <- c(degree, rbinom(1, init.deg[i], prob.network.recall))
    }
  }
    allsamples <- list(
      nsample = result[[9]],
      nominators = result[[11]],
      degsample = degree,
      attrsample = attrs[result[[9]]],
      time.recruited = result[[13]]
    )
  return(allsamples)
}

getsamples.RDS.C <- function(sim,
           dummy,
           M1,
           M2,
           N,
           nsamp0,
           fixinitial,
           n,
           coupons,
           dummydis,
           parallel = 0,
           parallel.type = "PSOCK",
           seed = 1,
           seed.indices = NULL,
           verbose = TRUE) {
    maxdeg <- length(dummy) / 2 - 1
    Msamples <- rep(0, length(dummy))
    bisamples <- matrix(0, length(dummy), length(dummy))
    sampbynet <- vector(M1, mode = "list")
    bsdegsamples <- matrix(0, ncol = n, nrow = M1 * M2)
    bsdissamples <- matrix(0, ncol = n, nrow = M1 * M2)
    bsprev <- rep(0, M1 * M2)
    parallel <- min(parallel, M1)
    #       parallel=1
    if (!is.null(fixinitial)) {
      if (is.logical(fixinitial[1]) &
          fixinitial[1] == "FALSE") {
        fixinitial <- -1
      } #line added 0527
    }
    data <- list(
      sim = sim,
      N = N,
      nsamp0 = nsamp0,
      fixinitial = fixinitial,
      M2 = M2,
      dummy = dummy,
      dummydis = dummydis,
      maxdeg = maxdeg,
      n = n,
      coupons = coupons,
      seed = seed,
      seed.indices = seed.indices
    )
    sfn <- function(j, data) {
      # SAN Next for sim version
      #         net<-as_network_uncompressed_rds(data$sim[[j]])
      net <- as_network_uncompressed_rds(data$sim[[1]])
      deg <- sapply(net$iel, length) + sapply(net$oel, length)
      deg[deg > data$maxdeg] <- data$maxdeg
      disease <- network::get.vertex.attribute(net, "disease")
      popCounts <- tototab.RDS(deg, disease, data$dummy)
      biCounts <- popCounts %*% t(popCounts)
      seed.distribution <- switch(
        as.character(data$fixinitial[1]),
        "-3" = {
          rep(1, length = network.size(net))
        },
        "-2" = {
          #                So attr(data$fixinitial,"table.seeds")=table(outcome[wave==0])
          if (any(data$dissample == 1) &
              any(disease == 0)) {
            table.seeds = attr(data$fixinitial, "table.seeds")
            disease * table.seeds[2] / sum(disease ==
                                             1) + (1 - disease) * table.seeds[1] / sum(disease == 0)
          } else{
            rep(1, length = network.size(net))
          }
        },
        "-1" = {
          deg
        },
        "1" = {
          disease
        },
        "2" = {
          disease * deg
        },
        "3" = {
          deg
        },
        deg
      )
      indices <- toindexref.RDS(deg, disease, data$dummy)
      u <- sort(unique(data$seed.indices))
      for (i in u) {
        seed.distribution[indices == i] <- 1000
      }
      si <- rep(0, data$M2 * data$nsamp0)
      for (k in 1:data$M2) {
        si[(k - 1) * data$nsamp0 + (1:data$nsamp0)] <-
          sample(
            seq(along = seed.distribution),
            size = data$nsamp0,
            replace = FALSE,
            prob = seed.distribution
          )
        for (i in u) {
          if (any(indices == i)) {
            if (sum(data$seed.indices == i) <= sum(indices == i)) {
              si[(k - 1) * data$nsamp0 + (1:data$nsamp0)][data$seed.indices == i] <-
                sample(
                  seq(along = indices)[indices == i],
                  size = sum(data$seed.indices == i),
                  replace = FALSE
                )
            } else{
              si[(k - 1) * data$nsamp0 + (1:data$nsamp0)][data$seed.indices == i] <-
                sample(
                  seq(along = indices)[abs(indices - i) < 2],
                  size = sum(data$seed.indices == i),
                  replace = FALSE
                )
            }
          }
        }
      }
      temp <-
        rdssampleC(
          net,
          trait.variable = "disease",
          nsamp0 = data$nsamp0,
          fixinitial = data$fixinitial,
          nsamp = data$n,
          replace = FALSE,
          coupons = data$coupons,
          nsims = data$M2,
          select = NULL,
          bias = NULL,
          seed.distribution = seed.distribution,
          seeds = si
        )
      temp$degsample[temp$degsample > data$maxdeg] <- data$maxdeg
      #counts of nodes of each type
      disCount <- rep(0, data$N + 1)
      Msamples <- rep(0, length(data$dummy))
      bisamples <- matrix(0, length(data$dummy), length(data$dummy))
      for (k in 1:data$M2) {
        ind <- (k - 1) * data$n + (1:data$n)
        tabsamp <-
          tototab.RDS(temp$degsample[ind], temp$attrsample[ind], data$dummy)
        disCount[sum(temp$attrsample[ind]) + 1] <-
          disCount[sum(temp$attrsample[ind]) + 1] + 1
        Msamples <- Msamples + tabsamp
        mbym <- tabsamp %*% t(tabsamp) - diag(tabsamp)
        bisamples <- bisamples + mbym
      }
      #         Q  <-  (Msamples+data$n/data$N)/(data$M2*popCounts+1)
      #         QQ <- (bisamples+(data$n/data$N)^2)/((data$M2)^2 * biCounts+1)
      Q  <-  (Msamples + 1) / (data$M2 * popCounts + 1)
      #         Q  <- estpi(n,Msamples,popCounts,tabsamp,dummydis,disease,M2,Q)
      QQ <- (bisamples + 1) / (data$M2 * biCounts + 1)
      #          Q[popCounts==0] <- 0
      #         QQ[ biCounts==0] <- 0
      bsprev <- rep(0, data$M2)
      temp2 <- matrix(0, length(data$dummy), length(data$dummy))
      bsdegsamples <-
        matrix(
          temp$degsample,
          ncol = data$n,
          nrow = data$M2,
          byrow = TRUE
        )
      bsdissamples <-
        matrix(
          temp$attrsample,
          ncol = data$n,
          nrow = data$M2,
          byrow = TRUE
        )
      for (k in 1:data$M2) {
        tabsamp <- tototab.RDS(bsdegsamples[k, ], bsdissamples[k, ], data$dummy)
        popests <- tabsamp / Q
        popests[Q == 0] <- 0
        popests <- data$N * popests / sum(popests, na.rm = TRUE)
        bsprev[k] <- sum(popests * data$dummydis, na.rm = TRUE) / data$N
      }
      #if(any(is.na(Q)) | any(is.na(bsprev))) {browser()}
      #if(any(is.na(Q)) | any(is.na(bsprev))) {print(Q);print(bsprev)}
      
      list(
        Msamples = Msamples,
        bisamples = bisamples,
        temp2 = temp2,
        Q = Q,
        QQ = QQ,
        bsprev = bsprev,
        disCount = disCount,
        bsdegsamples = bsdegsamples,
        bsdissamples = bsdissamples
      )
    }
    if (parallel > 1) {
      #      Run the jobs with rpvm or Rmpi
      cl <-
        beginparallel(parallel = parallel,
                              type = parallel.type,
                              seed = seed)
      if (verbose) {
        cat(paste("Starting parallel RDS samples in C\n"))
      }
      out <- parallel::clusterApplyLB(cl, as.list(1:M1), sfn, data)
      if (verbose) {
        cat(paste("Finished parallel RDS samples\n"))
      }
      endparallel(cl, type = parallel.type)
    } else{
      # non-parallel
      out <- vector(mode = "list", length = M1)
      if (verbose) {
        cat(paste("Starting non-parallel RDS samples in C\n"))
      }
      for (j in 1:M1) {
        out[[j]] <- sfn(j, data)
      }
      if (verbose) {
        cat(paste("Finished non-parallel RDS samples in C\n"))
      }
    }
    index <- 0
    for (j in 1:M1) {
      temp <- out[[j]]
      bsdegsamples[index + (1:M2), ] <- temp$bsdegsamples
      bsdissamples[index + (1:M2), ] <- temp$bsdissamples
      bsprev[index + (1:M2)] <- temp$bsprev
      if (j == 1) {
        Q <- temp$Q / M1
        QQ <- temp$QQ / M1
        disCount <- temp$disCount
      } else{
        Q <- Q + temp$Q / M1
        QQ <- QQ + temp$QQ / M1
        disCount <- disCount + temp$disCount
      }
      index <- index + M2
      Msamples <- Msamples + temp$Msamples
      bisamples <- bisamples + temp$bisamples
      sampbynet[[j]] <- temp$temp2
    }
    list(
      Msamples = Msamples,
      bisamples = bisamples,
      sampbynet = sampbynet,
      bsdegsamples = bsdegsamples,
      bsdissamples = bsdissamples,
      bsprev = bsprev,
      Q = Q,
      QQ = QQ,
      disCount = disCount
    )
  }
