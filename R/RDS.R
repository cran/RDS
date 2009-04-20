# Skeleton Package Generator.
        
   
# This function assumes that recruitment information is
# given by recruiter id. 
        
getWaveByID <- Vectorize(function(rds.data,recruitment.id)
   {
   row <- which(rds.data$recruitment.id == recruitment.id)
   recruiter.id <- rds.data[row,'recruiter.id']  
   if(recruiter.id == 0)
       {
       return(0)	
       }
    else
       {
       return(1 + getWaveByID(rds.data[,c('recruitment.id','recruiter.id')],recruiter.id))
       }
    		
   },vectorize.args = c("recruitment.id"))


# This function makes the standard assumptions.
        
count.transitions <- function(rds.data, group.variable)
   {
   	
   # Stop right now and cause an error if the group 
   # variable does not appear in the data frame.
   	
   stopifnot(group.variable %in% names(rds.data))
   	
   # The levels in the group variable are used to identify rows and
   # columns of the matrix of transition counts.  Therefore we
   # need to turn them into character strings.
   	
   group.names <- unique(as.character(rds.data[,group.variable]))
   
   # The transition matrix is the object that we will return.
   transition.matrix <- matrix(0,nrow=length(group.names),ncol=length(group.names))
   	
   colnames(transition.matrix) <- group.names
   rownames(transition.matrix) <- group.names	
   	
   #  The rest is easy.  For every row in the 
   # data frame we identify the recruiter id
   # if this id is not zero, then we increment
   # the count.  
   	
   for(idx in 1:nrow(rds.data))
     {	
     if(rds.data[idx,'recruiter.id'] != 0)
        {
        recruiter.row <- which(rds.data[idx,'recruiter.id'] == rds.data$id)[1]
        
        recruiter.group <- as.character(rds.data[recruiter.row,group.variable])
        recruit.group <- as.character(rds.data[idx,group.variable])
        transition.matrix[recruiter.group,recruit.group] <- 1 + transition.matrix[recruiter.group,recruit.group] 
        }	
     	
     }	
   	
   return(transition.matrix)	
   }
        


# This function assumes that transition counts can be turned
# into MLEs.
transition.counts.to.Markov.mle <- function(transition.counts)
    {
    tij <- transition.counts
    ti <- apply(tij,1,sum)
    stopifnot(ti>0)
    
    mle <- diag(1/ti) %*% tij    
    rownames(mle) <- rownames(tij)
    return(mle)	
    }


get.stationary.distribution <- function(mle)
    {

    	
    stat.dist <- NULL	
    	
    eigen.analysis <- eigen(t(mle))	

   # If -1 is an eigenvalue then we have a periodic chain ... which is bad.
   values.close.to.minus.one <- abs(eigen(t(mle))$values + 1.0) < sqrt(.Machine$double.eps)
   if(sum(values.close.to.minus.one) == 1)
          {
          cat(" Markov chain apears to be periodic.\n")
          return(stat.dist)
          }
   
       


    # The stationary distribution will be given by
    # the eigenvector with eigenvalue equal to one.    
    values.close.to.one <- abs(eigen(t(mle))$values - 1.0) < sqrt(.Machine$double.eps)
    
    # If there is exactly one of these, then we can take it as eigenvalue corresponding to the
    # stationary distribution.
    
    if(sum(values.close.to.one) == 1)
       {
       stationary.dist.eigen.vector <- eigen.analysis$vectors[,which(values.close.to.one)]
       stat.dist <- stationary.dist.eigen.vector/sum(stationary.dist.eigen.vector)
       names(stat.dist) <- colnames(mle)	
       }
   
    	
    return(stat.dist)	
    }




get.h.hat <- function(rds.data,group.variable,network.var='network.size')
   {
   	
   harmonic.mean <- function(x)
    {
    length(x)/sum(1/x) 	
    }
   	
   degrees <- rds.data[,network.var]
   temp <- function(x)
       {  
       deg <- degrees[rds.data[,group.variable] == x]
       return(harmonic.mean(deg))
       }	
       
   group.names <- unique(as.character(rds.data[,group.variable]))
   h.hat <- sapply(group.names,temp)
   names(h.hat) <- group.names
   return(h.hat)    	
   }

RDS.I.estimates <- function(rds.data,group.variable,network.variable)
    {
    tij <- count.transitions(rds.data,group.variable)
    markov.mle <- transition.counts.to.Markov.mle(tij)
    q.hat <- get.stationary.distribution(markov.mle)
    h.hat <- get.h.hat(rds.data,group.variable,network.variable)    
    
    return((q.hat/h.hat)/sum(q.hat/h.hat))	
    }


RDS.I.DS.estimates <- function(rds.data,group.variable,network.variable)
    {
    tij <- count.transitions(rds.data,group.variable)
    
    tij.ds <- 0.5*(tij + t(tij))
    
    markov.mle <- transition.counts.to.Markov.mle(tij.ds)
    q.hat <- get.stationary.distribution(markov.mle)
    h.hat <- get.h.hat(rds.data,group.variable,network.variable)    
    
    return((q.hat/h.hat)/sum(q.hat/h.hat))	
    }



RDS.II.estimates <- function(rds.data,trait.variable,network.variable)
   {
   outcome <- rds.data[,trait.variable]
   deg <- rds.data[,network.variable]	
   if(is.numeric(outcome))
      {
      result <- sum(outcome/deg)/sum(1/deg)	
      }
   else
      {
      group.names <- unique(as.character(outcome))
      numerator <- sapply(group.names,function(g){sum(1/deg[outcome == g])})
      names(numerator) <- group.names
      denominator <- sum(1/deg)
      result <- numerator/denominator
      }
   	
   return(result) 	
   }




