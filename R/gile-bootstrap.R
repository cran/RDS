
SSBS.estimates <- function(rds.data,trait.variable,
		number.of.bootstrap.samples=NULL,
		number.ss.samples.per.iteration=NULL,
		ties.to.trait=NULL,
		N=NULL, 
		confidence.level=NULL,
		fast=TRUE,
		weight.type="Gile's SS", verbose=TRUE)
{
	
	if(!(trait.variable %in% names(rds.data))){
		stop(sprintf("No variable called %s appears in the data.",
						trait.variable))
	}
	network.size <- attr(rds.data, "network.size.variable")
	
	if(is.null(N)){
		N <- attr(rds.data, "population.size.mid")
		if(is.null(N)){
			N <- ceiling(nrow(rds.data)/0.04)
			warning(paste("Parameter N missing, with no default for this data set. Using N =",N,"\n"))
		}
	}
	
	if(is.null(confidence.level)){
		confidence.level <- 0.95
	}
	if(is.null(number.of.bootstrap.samples)){ 
		number.of.bootstrap.samples <- 500
	}
	if(is.null(number.ss.samples.per.iteration)){
		number.ss.samples.per.iteration <- 500
	}
	
	trait.min.value <- NULL
	tv <- rds.data[[trait.variable]]
	#rds.data[[trait.variable]] <- tv     ###??? why? Are elements of trait.variable all numbers?
	remvalues <- !is.na(rds.data[[trait.variable]])
	if(sum(remvalues) < nrow(rds.data)){
		warning(paste(nrow(rds.data)-sum(remvalues),"of",nrow(rds.data),
						"trait variable values were missing and were removed."), call. = FALSE)
		if(is.numeric(rds.data[[trait.variable]])){
			trait.min.value <- min(rds.data[[trait.variable]],na.rm=TRUE)-1
			rds.data[[trait.variable]][!remvalues] <- trait.min.value
		}else{
			a=as.vector(rds.data[[trait.variable]])
			a[!remvalues] <- "NA"
			rds.data[[trait.variable]] <- factor(a,exclude=NULL)
			trait.min.value <- "NA"
		}
	}
	

	recruiter.id <- get.rid(rds.data)#as.character(rds.data[[attr(rds.data,"recruiter.id")]])
	seed.rid <- get.seed.rid(rds.data)
	outcome <- factor(rds.data[[trait.variable]],exclude=NULL)
	outclasses <- levels(outcome)
	deg <- rds.data[[network.size]]
	g <- length(outclasses)
	
	if(is.null(ties.to.trait)){
		# matrix of numbers of referrals to diseased or non-diseased nodes
		ties.to.trait<-matrix(0,nrow=nrow(rds.data),ncol=g)
		for(j in 1:g){
			# b is the number of referrals for each recruiter to those of group j
			b <- tapply(outcome==outclasses[j],recruiter.id,sum,na.rm=TRUE)[-1]
			# match the names of the recruiters to the id
			#d=match(names(b),as.character(rds.data[[attr(rds.data,"id")]]))
			d <- match(names(b),get.id(rds.data))
			# match the names of the recruiters to the id
			ties.to.trait[d[!is.na(d)],j] <- b[!is.na(d)]
		}
	}else{
		if(is.character(ties.to.trait)){
			TT<-matrix(0,nrow=nrow(rds.data),ncol=g)
			for(j in 1:g){
				TT[,j] = as.vector(rds.data[[ties.to.trait[j]]])
			}
			ties.to.trait<-TT
			rm(TT)
		}
	}
	
	n0 <- sum(recruiter.id==seed.rid, na.rm=TRUE)
	

	
	###########################################################################################
# Check for missing values and warn the user if any are removed.   This should really taken
# care of elsewhere.  NB: It is also worth considering the semantics of the message 
# "values were missing and were removed".    
	
	remvalues <- !is.na(as.vector(outcome))
	remvalues[as.vector(outcome)=="NA"] <- FALSE
	if(sum(remvalues) < nrow(rds.data)){
		warning(paste(nrow(rds.data)-sum(remvalues),"of",nrow(rds.data),
						"values were missing and were removed."))
		deg <- deg[remvalues]
		outcome <- factor(as.vector(outcome[remvalues]),exclude=NULL)
		outclasses <- levels(outcome)
		g <- length(outclasses)
		ties.to.trait <- ties.to.trait[remvalues,]
		recruiter.id <- recruiter.id[remvalues]
		n0 <- sum(recruiter.id==seed.rid, na.rm=TRUE)
	}
	
	if(length(unique(outcome))>1){
		result<-sppsboot4ppsreal(degs=deg,dis=outcome,n=N, n0=n0,
				ties.to.trait=ties.to.trait,
				nit=number.of.bootstrap.samples,
				number.ss.samples.per.iteration=number.ss.samples.per.iteration,
				fast=fast,weight.type=weight.type,verbose=verbose) 
		names(result$point_estimate) <- outclasses
		names(result$se_estimate) <- outclasses
		colnames(result$bsests) <- outclasses
	}else{
		result<-list(point_estimate=c(0,1),
				bsests=cbind(rep(1,number.of.bootstrap.samples),
						rep(1,number.of.bootstrap.samples)),
				se_estimate=c(0,0))
		names(result$point_estimate) <- c("0",outclasses)
		names(result$se_estimate) <- c("0",outclasses)
		colnames(result$bsests) <- c("0",outclasses)
	}
	
	observed.estimate <- result$point_estimate
	crit <- qnorm((1-confidence.level)/2,lower.tail=FALSE)
	a <- cbind(observed.estimate,
			observed.estimate-crit*apply(result$bsests, 2, sd,na.rm=TRUE),
			observed.estimate+crit*apply(result$bsests, 2, sd,na.rm=TRUE))
	colnames(a)[1] <- trait.variable
	
	attr(a,"bsresult") <- result
	return(a)
}



#was vh_est
# computes horvitz-thompson estimator
vh.est<-function(degs,dis,wstart=0,wsample=NULL){
	degs[degs==0]<-1 #added 0707
	if(wstart>0){
		degsnew<-degs[wsample>=wstart]
		disnew<-dis[wsample>=wstart]		
	}else{
		degsnew<-degs
		disnew<-dis
	}
	if(is.factor(dis)){
		num <- rep(0,length=length(levels(disnew)))
		a <- tapply(1/degsnew,as.numeric(disnew),sum)
		num[as.numeric(names(a))] <- a
	}else{
		num<-sum(disnew/degsnew)
	}
	den<-sum(1/degsnew)
	num/den
}
#was spps_est_keep
spps.est.keep<-function(degs,dis,nguess,wstart=0,wsample=NULL,number.ss.iterations=5,number.ss.samples.per.iteration=2000){
	degs[degs==0]<-1 
	mapping<-getestCstacked(degs,n=nguess,nit=number.ss.iterations,nsampsamp=number.ss.samples.per.iteration,trace=FALSE)
	weights=approx(x=mapping$classes,y=1/mapping$probs,xout=degs,rule=2)$y
	pis <- 1/weights
	if(is.factor(dis)){
		num <- rep(0,length=length(levels(dis)))
		a <- tapply(1/degs,as.numeric(dis),sum)
		num[as.numeric(names(a))] <- a
	}else{
		num<-sum(dis/degs)
	}
	est<-num/sum(1/degs)
	
	list(samplewts=pis,
			est=est, 
			classes=mapping$classes, pvec=mapping$probs,
			mapping=mapping) 
}

#was ssps_est
spps.est<-function(degs,dis,nguess,wstart=0,wsample=NULL,nsampsamp=500,hajek=TRUE, mapping=NULL){
	degs[degs==0]<-1 #added 070708
	if(is.null(mapping)){
		mapping<-getestCstacked(degs,n=nguess,nit=3,nsampsamp=nsampsamp,trace=FALSE,hajek=hajek)
	}
	pis=1/approx(x=mapping$classes,y=1/mapping$probs,xout=degs,rule=2)$y
	vh.est(pis,dis,wstart=wstart,wsample=wsample)
}

######################################
#####  Main Function #################
######################################

# compute the bootstrap.  Returns the point estimate and bootstrap samples. We then take the sd of the bootstrap estimates as the estimate of the se (in simulations, this performs better than the quantiles of the bootstrap estimates)

#############  Arguments:
# degs = degrees of sample units
# n = population size,                         
# n0 = number of seeds, 
# dis is a nsamp vector with elements 1,...,g 
# ties.to.trait is a nsamp*g matrix, where ties.to.trait[i,j]=number of referrels from sample i to recruit with z=j.
sppsboot4ppsreal<-function(degs,dis,n,n0,ties.to.trait,nit,                                   # ignore fullcup
		fullcup=rep(TRUE,length(degs)),
		number.ss.samples.per.iteration=500,
		fast=TRUE,
		weight.type="Gile's SS",verbose=TRUE){ 
	fdis <- dis
	dis <- as.numeric(dis)
	outclasses <- levels(fdis)
	g=length(outclasses)
	nsamp<-length(degs)
	nrefs<-apply(ties.to.trait,1,sum)[which(fullcup)]
	if(n*0.04 > nsamp){n <- round(nsamp/0.04)}  
	aaa<-spps.est.keep(degs,fdis,n,number.ss.samples.per.iteration=number.ss.samples.per.iteration) 
	wts<-aaa$samplewts
	est<-aaa$est         
	pvec<-aaa$pvec
	if(fast){
		data.mapping<-aaa$mapping
	}else{
		data.mapping<-NULL
	}
	classes<-sort(unique(degs))
	
	K=max(degs)
	classesbotha<-classes # Create unique labels for each (deg, dis) class
	for(i in 2:g)
		classesbotha<-c(classesbotha,classes+(i-1)*K)
	#g times the length of the vector of degree classes, to allow for degree-(g)disease categories. (We can assume the disease has g status ^_^)
	
	samplecountsa<-c(table(degs,dis))
	wtsbotha<-rep(pvec,g)
	
	CC=matrix(0,g,g)
	for(i in 1:nsamp){
		for(j in 1:g){
			CC[dis[i],j]=CC[dis[i],j]+ties.to.trait[i,j]
		}
	}
	# CC is a g*g matrix, where CC[i,j]=number of referrals from recruiter with z=i to recruit with z=j.
	
	deg=numeric()
	for(i in 1:g){
		deg[i]<-vh.est(wts[dis==i],degs[dis==i])
	}
	# deg[i] is the estimated mean degrees of nodes with z=i
	
	R = sweep(CC,1,apply(CC,1,sum),"/")
	R[is.na(R)] <- 0
	# R is a g*g matrix, where R[i,j]=observed rates of referral of nodes with z=i for z=j.
	
	tiemtx = sweep(R,1,deg*est*n,"*")
	# tiemtx is a matrix with the estimated number of edges between each pair of classes
	
	classesboth<-classesbotha[samplecountsa>0]
	samplecounts<-samplecountsa[samplecountsa>0]   
	wtsboth<-wtsbotha[samplecountsa>0]
	
	manynewests<-matrix(0,nrow=nit,ncol=g)
	manysamp<-vector("list",length=nit)
	
	if(weight.type=="Gile's SS"){
		effn <- n
	}else{
		effn <- 10^8
	}
	
	data<-list(n=n, effn=effn,
			classesboth=classesboth,samplecounts=samplecounts, wtsboth=wtsboth,
			n0=n0, nsamp=nsamp, tiemtx=tiemtx, nrefs=nrefs,
			mapping=data.mapping)  
	
	bsfn <- function(i, data){
		newprops<-probtodist(data$classesboth,data$samplecounts, data$wtsboth, data$n)$props
		props2<-newprops/sum(newprops)
		nbyclass<-round(data$n*props2)  
		nbyclass[nbyclass==0]<-1
		offby<-data$n-sum(nbyclass)
		if(is.na(offby)){print("Error in getincl: offby"); } 
		
		if(offby>0){
			for(ii in 1:offby){
				dif<-props2-nbyclass/data$n
				dif[dif<0]<-0
				tochange<-sample(1:length(dif),1,prob=dif)
				nbyclass[tochange]<-nbyclass[tochange]+1
			}
		}else{if(offby<0){
				for(ii in 1:abs(offby)){                
					dif<-nbyclass/data$n - props2
					dif[dif<0]<-0
					dif[nbyclass==1]<-0 
					if(sum(dif==0)){
						dif<-1-(props2-nbyclass/data$n)
						dif[nbyclass==1]<-0
					}       
					tochange<-sample(1:length(dif),1,prob=dif)
					nbyclass[tochange]<-nbyclass[tochange]-1
				}
			}}   
		
		popclass<-rep(data$classesboth,times=nbyclass)    # Create a population of labels (for each (deg, dis) class)
		newdis<-((popclass-1)%/%K)+1  
		newdeg<-round(popclass-K*(popclass%/%K))
		newdeg[newdeg==0]<-K
		
		nsamplea<-sample(1:length(newdeg),data$n0,prob=abs(newdeg),replace=FALSE)
		todisnew=matrix(0,nrow=data$nsamp,ncol=g)
		
		tempties <- (data$tiemtx+t(data$tiemtx))/2
		
		nprop <- tabulate(newdis,nbins=g)
		nprop <- nprop / sum(nprop)
		for(k in 1:g){
			if(all(tempties[k,]==0)){
				tempties[k,] <- nprop
			}
			if(all(tempties[,k]==0)){
				tempties[,k] <- nprop
			}
		}
		
		nrefs<-data$nrefs 
		i<-sample(1:length(nrefs),1)                  
		numberfrom<-nrefs[i] 
		while(numberfrom==0){
			i<-sample(1:length(nrefs),size=1) 
			numberfrom<-nrefs[i]
		} 
		
		nsample<-c(nsamplea,rep(0,(data$nsamp-data$n0)))    ### Vector W !!!
		activnode<-1
		countrefs<-0
		temptiess <- apply(tempties,1,sum)
		for(mm in (data$n0+1):data$nsamp){                 ### loop begins!!!  mm is not m in the paper!
			k = newdis[nsample[activnode]]
			p <- tempties[k,] / temptiess[k]
#        ## p is the vector of probability

			nextdis<-sample(x=1:g,size=1,replace=FALSE,prob=tempties[k,])
			
#       numsfrom is the index of the "nextdis" disease status left .
			numsfrom<-intersect(
					which(newdis==nextdis),
					c(1:length(newdeg))[-nsample[1:(mm-1)]]
			)
			
			if(length(numsfrom)==0){
				numsfrom<-c(1:length(newdeg))[-nsample[1:(mm-1)]] 
				warning(paste("Ran out of",nextdis))
			}
			if(length(numsfrom)==1){
				nsample[mm]<-numsfrom
			}else{     
				# Draw a node those with newdis left proportional to degree
				nsample[mm]<-sample(numsfrom,size=1,prob=newdeg[numsfrom])
			}
			
			totaltmp<-sum(newdeg[numsfrom])  # the sum of degrees left with newdis
			
			for(k in 1:g){
				if(newdis[nsample[mm]]==k){ 
					todisnew[activnode,k]=todisnew[activnode,k]+1
					tdeg<-newdeg[nsample[mm]] # the degree of the sample
					tempties[,k]<-tempties[,k]*(totaltmp-tdeg)/totaltmp       ###??? Shouldn't tempties be symmetric?
					
				}
			}
			countrefs<-countrefs+1
			if((mm<nsamp)&(countrefs==numberfrom)){
				activnode<-activnode+1                     ### move to the next seed (or node)!!! i=i+1
				countrefs<-0 
				i<-sample(1:length(nrefs),size=1) 
				numberfrom<-nrefs[i] 
				while(numberfrom==0){
					i<-sample(1:length(nrefs),size=1) 
					numberfrom<-nrefs[i] 
				}
			}
		}                                            ### loop ends!!!
		
		degsample<-newdeg[nsample]
		dissample<-factor(newdis[nsample],
				levels=1:length(outclasses),labels=outclasses,exclude=NULL)
		manysamp<-list(deg=degsample,dis=dissample,ties.to.trait=todisnew)
#THE ORIGINAL ONE:
		manynewests<-spps.est(degsample,
				dissample,
				data$effn,nsampsamp=ceiling(number.ss.samples.per.iteration/4), 
				mapping=data$mapping)
  
		list(samp=manysamp, newests=manynewests)
		
	}
	if(verbose)  cat(paste('Computation 0% completed ...\n',sep=""))
	for(i in 1:nit){
		out <- bsfn(i,data)
		manysamp[[i]] <- out[["samp"]]
		manynewests[i,] <- out[["newests"]]
		
		if(verbose) {
			if(i == trunc(i/(nit/10))*(nit/10)){
				cat(paste((100/10)*trunc(i/(nit/10)),'% completed ...\n',sep=""))
			}
		}
	}

	list(point_estimate=aaa$est,bsests=manynewests,
			se_estimate=apply(manynewests,2,sd))
}





