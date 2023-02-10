#' @keywords internal
beginparallel <- function(parallel = 1,
           type = "PSOCK",
           seed = NULL,
           packagenames = c("RDS"),
           verbose = TRUE) {
    ### Snow commands to set up cluster
    cl <- parallel::makeCluster(parallel, type = type)
    ### initialize parallel random number streams
    if (is.null(seed)) {
      parallel::clusterSetRNGStream(cl)
    } else{
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    ### start each virtual machine with libraries loaded
    for (pkg in packagenames) {
      attached <-
        parallel::clusterCall(cl,
                              require,
                              package = pkg,
                              character.only = TRUE)
    }
    return(cl)
  }

#' @keywords internal
endparallel <- function(cl,
                        type = "PSOCK",
                        finalize = FALSE,
                        verbose = TRUE) {
  parallel::stopCluster(cl)
  invisible()
}
