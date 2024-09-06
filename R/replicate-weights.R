#' HCG parametric bootstrap replicate weights
#' @param rds.data An rds.data.frame
#' @param outcome.variable The column name of the variable defining the groups for the homophily configuration graph
#' @param number.of.bootstrap.samples The number of bootstrap replicate weights to be generated
#' @param include.sample.weights If TRUE, the first column of the returned frame are the HCG weights for the sample
#' @param N The population size
#' @param small.fraction If TRUE, the sample size is assumed to be small compared to the population size
#' @details
#' This function generates bootstrap replicate weights which may be used to
#' analyze RDS data in other packages or software systems (e.g. the survey package
#' with svrepdesign).
#' 
#' @returns A data.frame of replicate weights. If include.sample.weights is TRUE, the
#' first column are the HCG weights for the observed sample.
#'
#' @examples
#' data("fauxmadrona")
#' set.seed(1)
#' # Generate replicate weights
#' result <- hcg.replicate.weights(fauxmadrona, "disease", 500, TRUE)
#' # Analyze with survey package and compare to internal function
#' if(require(survey)){
#'   set.seed(1)
#'   design <- svrepdesign(fauxmadrona, type = "bootstrap", 
#'     weights= result[[1]], repweights = result[-1])
#'   svymean(~disease, design) |> print()
#'   RDS.bootstrap.intervals(fauxmadrona, "disease", "HCG", "HCG", 
#'   number.of.bootstrap.samples = 500) |> print()
#' }
#' 
#' @export
hcg.replicate.weights <- function(rds.data, outcome.variable, 
                                  number.of.bootstrap.samples = 500, 
                                  include.sample.weights = FALSE,
                                  N=NULL, 
                                  small.fraction = FALSE){
  if(is.null(N)){
    N <- attr(rds.data, "population.size.mid")
    if(is.null(N)){
      small.fraction <- TRUE
    }else{
      cat("\nNote: Using the data's mid population size estimate: N =",
          N,
          "\n")
    }
  }
  
  # Handle missing degrees
  net.size.variable <- attr(rds.data, "network.size.variable")
  degree <- get.net.size(rds.data)
  mdeg <- is.na(degree)
  out <- rds.data[[outcome.variable]]
  if(any(mdeg)){
    warning(paste(sum(mdeg),"subjects missing degree. Imputing within group median."))
    levs <- na.omit(unique(out))
    med <- median(degree, na.rm = TRUE)
    for(i in 1:length(levs)){
      group.median <- median(degree[out==levs[i]], na.rm=TRUE)
      degree[out==levs[i] & mdeg] <- group.median
    }
    mdeg <- is.na(degree)
    if(any(mdeg)){
      degree[mdeg] <- med
    }
    rds.data[[net.size.variable]] <- degree
  }
  
  # Calculate replicate weights from bootstrap sample
  replicate.weights <- function(x){
    g <- rds.data[[outcome.variable]]
    d <- get.net.size(rds.data)
    gb <- x[[outcome.variable]]
    db <- get.net.size(x)
    wts <- hcg.weights(x, outcome.variable, N=N, small.fraction = small.fraction, weights.include.seeds = TRUE)
    wts <- wts / sum(wts)
    wts2 <- rep(0,length(g))
    for(i in 1:length(gb)){
      if(wts[i]  < .Machine$double.eps)
        next
      ids <- which(g == gb[i] & d == db[i])
      id <- ids[sample.int(length(ids), 1)]
      wts2[id] <- wts2[id] + wts[i]
    }
    wts2
  }
  
  res <- as.data.frame(HCG.bootstrap(rds.data, outcome.variable, 
                                    number.of.bootstrap.samples, 
                                    fun = replicate.weights,
                                    N=N))
  names(res) <- paste0("rep_wts",1:ncol(res))
  if(include.sample.weights){
    wts <- hcg.weights(rds.data, outcome.variable, N=N, small.fraction = small.fraction, 
                       weights.include.seeds = TRUE)
    wts <- wts / sum(wts)
    res <- cbind(wts=wts, res)
  }
  res
}
