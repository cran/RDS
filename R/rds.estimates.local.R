

RDS.estimates.local <- function(rds.data,outcome.variable,subset=NULL,
		empir.lik=TRUE, weight.type, N=NULL, ...){
	
	if(is(rds.data,"rds.data.frame")){
		
		if(!(outcome.variable %in% names(rds.data))){
			stop(sprintf("No variable called %s appears in the data.",
							outcome.variable))
		}
		network.size <- attr(rds.data, "network.size.variable")
		
	}
	else{
		stop("rds.data must be of type rds.data.frame")
	}
	
	if(weight.type %in% c("RDS-I","RDS-I (DS)"))
		rds.data[[outcome.variable]] <- factor(rds.data[[outcome.variable]])
	
	
	if(is.null(N)){
		N <- attr(rds.data, "population.size.mid")
	}
	
	###########################################################################################
	# Check for missing values and warn the user if any are removed.   This should really taken
	# care of elsewhere.  NB: It is also worth considering the semantics of the message 
	# "values were missing and were removed".    
	remvalues <- rds.data[[network.size]]==0 | is.na(rds.data[[network.size]])
	if(any(remvalues)){
		warning(paste(sum(remvalues),"of",nrow(rds.data),
						"network sizes were missing or zero. The estimator will presume these are",max(rds.data[[network.size]],na.rm=TRUE)), call. = FALSE)
		
		rds.data[[network.size]][remvalues] <- max(rds.data[[network.size]],na.rm=TRUE)
	}
	
	
	rds.data.nomiss <- rds.data
	
	#########################################################################################
	# The simple data management tasks have been taken care of, so it's now time to compute 
	# the RDS estimates.   The cases for numeric and categorical outcomes are handled 
	# separately.
	
	if(is.null(subset)){
		subset <- rep(TRUE,length=nrow(rds.data.nomiss))
	}else{
		subset[is.na(subset)] <- FALSE
		if(!is.null(N)){
			#use VH estimator to adjust population size to sub-population
			tmp.wts <- vh.weights(rds.data.nomiss[[network.size]])
			tmp.wts <- tmp.wts/sum(tmp.wts)
			prop <- sum(tmp.wts*subset)
			N <- N * prop
			if(N < sum(subset))
				stop("Estimated sub-population size smaller than subset.")
		}
		
		#This subsets setting orphaned children to be seeds
		#in order to maintain a valid recruitment tree.
		rds.data.nomiss <- rds.data.nomiss[subset,,warn=FALSE]
		
		#drop 0 count levels
		if(is.factor(rds.data[[outcome.variable]]))
			rds.data.nomiss[[outcome.variable]] <- factor(rds.data.nomiss[[outcome.variable]])
	}
	weights.all <- compute.weights(rds.data.nomiss,
			weight.type=weight.type,outcome.variable=outcome.variable, N=N, ...)
	
	
	#########################################################################################
	# The simple data management tasks have been taken care of, so it's now time to compute 
	# the RDS-II estimates.   The cases for numeric and categorical outcomes are handled 
	# separately.
	outcome <- rds.data.nomiss[[outcome.variable]]
	toutcome <- table(outcome)
	if(is.numeric(outcome)){
		estimate <- HT.estimate(weights=weights.all,outcome=outcome)
		if(empir.lik)
			attr(estimate,"se") <- EL.se(weights.all,outcome,N=N)
	}else{
		outcome <- as.factor(outcome)
		
		# Calculate the within group totals of the *inverse* of the weights, again we exclude the seeds. 
		wtd.degree.totals <- xtabs(weights.all~outcome,
				data=rds.data.nomiss)
		
		estimate.xtabs <- wtd.degree.totals/sum(wtd.degree.totals)
		# Now convert to a numeric vector ('estimate.xtabs' is currently of class 'xtabs')
		estimate <- as.numeric(estimate.xtabs)
		names(estimate) <- names(estimate.xtabs)
		names(estimate) <- names(estimate.xtabs)
		if(empir.lik){
			aaa <- levels(outcome)
			bbb = rep(NA,length(aaa))
			for(i in seq_along(aaa)){
				bbb[i] <- EL.se(weights.all,outcome==aaa[i],N=N)
			}
			attr(estimate,"se") <- bbb
		}
	}
	
	
	if(is.null(attr(estimate,"se"))){
		result <- rds.weighted.estimate(estimate,weights.all,outcome.variable,weight.type,subset=subset)
	}else{ 
		
		nsamples <- sum(!is.na(as.vector(outcome)))
		estimate <- cbind(
				estimate,
				estimate-1.96*attr(estimate,"se"),
				estimate+1.96*attr(estimate,"se")
		)
		colnames(estimate) <- c("point", "lower", "upper")		
		if(is.numeric(outcome)){
			#
			# Design effect (for numeric outcomes)
			#
			varoutcome <- var(outcome)
			if(is.null(N)){
				varsrs <- varoutcome/nsamples
			}else{
				varsrs <- (((N - nsamples)/(N - 1)) * varoutcome/nsamples)
			}
			de <- ((estimate[, 3] - estimate[, 2])/(1.96 * 2))^2/varsrs
			estimate <- cbind(estimate, de, (estimate[, 3] - estimate[, 2])/(1.96 * 2),nsamples)
			colnames(estimate)[c(4, 5, 6)] <- c("Design Effect", "s.e.", "n")
			rownames(estimate) <- outcome.variable
			names(estimate) <- outcome.variable
		}else{
			nsamplesbyoutcome <- table(suppressWarnings(as.numeric(outcome)))
			#
			# Design effect (for factor outcomes)
			#
			varoutcome <- estimate[, 1] * (1 - estimate[, 1])
			if(is.null(N)){
				varsrs <- varoutcome/nsamples
			}else{
				varsrs <- (((N - nsamples)/(N - 1)) * varoutcome/nsamples)
			}
			de <- ((estimate[, 3] - estimate[, 2])/(1.96 * 2))^2/varsrs
			estimate <- cbind(estimate, de, (estimate[, 3] - estimate[, 2])/(1.96 * 2), nsamplesbyoutcome)
			colnames(estimate)[c(4, 5, 6)] <- c("Design Effect", "s.e.", "n")
			names(estimate) <- levels(outcome)
		}
		result <- rds.interval.estimate(estimate, outcome.variable, 
				weight.type=weight.type, uncertainty="EL", weights=weights.all,N=N)		
	}
	return(result)
}
