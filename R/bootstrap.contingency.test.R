#' Performs a bootstrap test of independance between two categorical variables
#' @param rds.data an rds.data.frame
#' @param row.var the name of the first categorical variable
#' @param col.var the name of the second categorical variable
#' @param number.of.bootstrap.samples The number of simulated boootstrap populations
#' @param table.only only returns the weighted table, without bootstrap.
#' @param verbose level of output
#' @param weight.type The type of weighting to use for the contningency table. Only large sample methods are allowed.
#' @param ... Additional parameters for compute_weights
#' @details 
#' This function first estimates a Homophily Configuration Graph model for the underlying
#' network under the assumption that the two variables are independant and that the population size is large. 
#' It then draws bootstrap
#' RDS samples from this population distribution and calculates the chi.squared statistic on
#' the weighted contingency table. Weights are calculated using the HCG estimator assuming a large population size.
#' @examples
#' data(faux)
#' bootstrap.contingency.test(rds.data=faux, row.var="X", col.var="Y",
#'   number.of.bootstrap.samples=50, verbose=FALSE)
#' @export
bootstrap.contingency.test <- function(rds.data,
                                       row.var,
                                       col.var,
                                       number.of.bootstrap.samples = 1000,
                                       weight.type=c("HCG", "RDS-II","Arithmetic Mean"),
                                       table.only=FALSE,
                                       verbose = TRUE, ...) {
  if(missing(weight.type)){
    weight.type <- "HCG"
  }
  weight.type <- match.arg(weight.type, c("RDS-II","Arithmetic Mean","HCG"))
  
  network.size <- attr(rds.data, "network.size.variable")
  remvalues <- rds.data[[network.size]] == 0 | is.na(rds.data[[network.size]])
  if (any(remvalues)) {
    warning(
      paste(
        sum(remvalues),
        "of",
        nrow(rds.data),
        "network sizes were missing or zero. The estimator will presume these are",
        max(rds.data[[network.size]], na.rm = TRUE)
      ),
      call. = FALSE
    )
    rds.data[[network.size]][remvalues] <- max(rds.data[[network.size]], na.rm = TRUE)
  }
  
  #Set up
  N <- nrow(rds.data) * 100000
  wave <- get.wave(rds.data)
  if (!has.recruitment.time(rds.data)) {
    time <- rep(1, nrow(rds.data))
  } else{
    time <- get.recruitment.time(rds.data)
  }
  rds <- rds.data[order(time, wave),]
  if (!has.recruitment.time(rds)) {
    time <- rep(1, nrow(rds.data))
  } else{
    time <- get.recruitment.time(rds)
  }
  
  # Create combination variable of both row and column
  row.fact <- as.factor(rds[[row.var]])
  row.nms <- levels(row.fact)
  row <- as.numeric(row.fact)
  nr <- max(row, na.rm=TRUE)
  col.fact <- as.factor(rds[[col.var]])
  col.nms <- levels(col.fact)
  col <- as.numeric(col.fact)
  nc <- max(col, na.rm=TRUE)
  eg <- expand.grid(1:nc, 1:nr)
  v <-
    factor(paste(row, col, sep = "_"), levels = paste(eg$Var2, eg$Var1, sep =
                                                        "_"))
  out.miss <- is.na(row) | is.na(col)
  varname <- .make.unique(names(rds), "variable")
  rds[[varname]] <- v
  
  if(length(row.nms)<2){
    stop(paste(row.var,"contains less than 2 unique values"))
  }

  if(length(col.nms)<2){
    stop(paste(col.var,"contains less than 2 unique values"))
  }  
  # Construct population null distribution where row and column are independant using HCG
  rid <- get.rid(rds)
  id <- get.id(rds)
  degree <- get.net.size(rds)
  hcg.row <-
    hcg.estimate(id, rid, time, degree, row, N, small.fraction = TRUE)
  hcg.col <-
    hcg.estimate(id, rid, time, degree, col, N, small.fraction = TRUE)
  theta <- matrix(0, nrow = nr * nc, ncol = nr * nc)
  yhat <- rep(0, nr * nc)
  for (i in 1:nr) {
    for (j in 1:nr) {
      for (k in 1:nc) {
        for (l in 1:nc) {
          theta[(i - 1) * nc + k, (j - 1) * nc + l] <-
            hcg.row$theta[i, j] * hcg.col$theta[k, l]
          yhat[(i - 1) * nc + k] <-
            hcg.row$yhat[i] * hcg.col$yhat[k]
        }
      }
    }
  }
  #dbar <- tapply(degree,v,function(x) length(x) / sum(1/x))
  
  # Calculate weights from HCG model
  names(yhat) <-
    colnames(theta) <-
    row.names(theta) <- paste(eg$Var2, eg$Var1, sep = "_")
  weights <- rep(0, nrow(rds))
  for (i in 1:nr) {
    for (k in 1:nc) {
      d <- degree[row == i & col == k & !out.miss]
      weights[row == i & col == k & !out.miss] <- yhat[(i - 1) * nc + k] * (1 / d) / sum(1 / d)
    }
  }
  hcg.est <- list(theta = theta,
                  yhat = yhat,
                  weights = weights)
  
  #Calculates the chi squared statistic on the weighted contingency table for RDS dataset 'x'
  chi.squared.func <- function(x, with.names=FALSE) {
    vspl <- strsplit(as.character(x[[varname]]), "_", fixed = TRUE)
    v1 <- as.numeric(sapply(vspl, function(y)
      y[1]))
    v2 <- as.numeric(sapply(vspl, function(y)
      y[2]))
    x[[varname]] <- factor(x[[varname]])
    wts <- compute.weights(
      x,
      outcome.variable = varname,
      weight.type = weight.type,
      N = N,
      small.fraction = TRUE
    )
    wts <- wts * nrow(x) / sum(wts)
    tab <- tapply(wts, list(v1, v2), sum, simplify = TRUE)
    if(with.names){
      rownames(tab) <- row.nms[as.numeric(rownames(tab))]
      colnames(tab) <- col.nms[as.numeric(colnames(tab))]
    }
    tab[is.na(tab)] <- 0
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop=FALSE]
    if(is.null(dim(tab)) || ncol(tab) < 2 || nrow(tab) < 2)
      result <- structure(NA, .Names = "X-squared")
    else
      result <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
	  attr(result,"table") <- tab
	  result
  }
  
  stat <- chi.squared.func(rds, with.names=TRUE)
  tbl <- attr(stat,"table")
  attr(stat,"table") <- NULL
  if(table.only){
    return(tbl)
  }
  
  # Perform bootstrap
  boot.stats <-
    unlist(
      HCG.bootstrap(
        rds,
        varname,
        number.of.bootstrap.samples,
        N = N,
        fun = chi.squared.func,
        hcg.est = hcg.est,
        small.fraction = TRUE,
        verbose = verbose
      )
    )

  pvalue <- mean(boot.stats > stat, na.rm=TRUE)
  result <- list(
    table = tbl,
    p.value = pvalue,
    statistic = stat,
    method = paste("RDS Bootstrap Test of", row.var, "versus", col.var),
    row.var = row.var,
    col.var = col.var,
    boot.stats = boot.stats
  )
  class(result) <- c("rds.contin.bootstrap", "htest")
  result
}

#' Displays an rds.contin.bootstrap
#' @param x an rds.contin.bootstrap object
#' @param show.table Display weighted contingency table
#' @param ... additional parameters passed to print.matrix.
#' @export
#' @method print rds.contin.bootstrap
print.rds.contin.bootstrap <- function(x, show.table=FALSE, ...){
  if(show.table){
    cat("Weighted table:", x$row.var, " by ", x$col.var,"\n")
    print(x$table, ...)
  }
  class(x) <- "htest"
  print(x)
}

.make.unique <- function(names, name) {
  i <- 0
  if (!(name %in% names))
    return(name)
  while (TRUE) {
    i <- i + 1
    newname <- paste0(name, i)
    if (!(newname %in% names))
      return(newname)
  }
}
