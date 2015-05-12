#' Good-Fellows weights
#' @param rds.data An rds.data.frame
#' @param outcome.variable The variable used to base the weights on.
#' @param N Population size
#' @param smoothed Should the data smoothed Good-Fellows weights be computed.
#' @param ... Unused
#' @export
gf.weights<-function(rds.data, outcome.variable, N=NULL,smoothed=FALSE,...){
	if(is.null(rds.data[[outcome.variable]])){
		stop(paste("Good-Fellows outcome.variable", outcome.variable,"not present in data"))
	}
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
		
		d=tapply(rds.data[[outcome.variable]],rds.data[[outcome.variable]],length)
		f=match(rds.data[[outcome.variable]],names(est))
		weights <- rep(0,length(f))
		weights[!is.na(f)]=as.numeric(unlist(est[f]))
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

