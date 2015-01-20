memle <-function(network.size,num.recruits,recruit.time,max.coupons=3,unit.scale=NULL,

          unit.model="nbinom",

          cutoff=0,cutabove=1000,

          guess=c(3,0,0,0.6,1.5,0),

          method="BFGS", hessian=TRUE, K=max(100,network.size), verbose=FALSE){

 logit <- function(p){log(p/(1-p))}

 if(is.null(unit.scale)){

   unit.scale <- -1

   vtrans <- function(l){c((exp(l[1])+1)/exp(l[1]),l[2],l[3],exp(l[4]),exp(l[5])+1,exp(l[6]))}

   ltrans <- function(v){c(-log(v[1]-1),v[2],v[3],log(v[4]),log(v[5]-1),log(v[6]))}

   gtrans <- function(l){c(-exp(-l[1]),1,1,exp(l[4]),exp(l[5]),exp(l[6]))}

 }else{

   unit.scale <- 1

   vtrans <- function(l){c((exp(l[1])+1)/exp(l[1]),l[2],l[3],exp(l[4]),exp(l[5]))}

   ltrans <- function(v){c(-log(v[1]-1),v[2],v[3],log(v[4]),log(v[5]))}

   gtrans <- function(l){c(-exp(-l[1]),1,1,exp(l[4]),exp(l[5]))}

 }

 n <- length(network.size)

#

 llme <- switch(unit.model,

  "cmp"=llcmpme, "nbinom"=llnbme)

 if(sum(network.size>=cutoff & network.size <= cutabove) > 0){

  aaa <- optim(par=ltrans(guess),fn=llme,

   method=method,

   hessian=hessian,control=list(fnscale=-10),

   n=n,

   recruit.time=recruit.time,

   num.recruits=num.recruits,network.size=network.size,

   max.coupons=max.coupons,K=K,

   unit.scale=unit.scale, verbose=verbose)

   if(unit.scale<0){

    if(unit.model=="cmp"){

     names(aaa$par) <- c("CMP mean","Recruitment Odds","Recruitment Odds Time","Error log-s.d.", "CMP skew","Optimism")

    }else{

     names(aaa$par) <- c("Neg.Bin. mean","Recruitment Odds","Recruitment Odds Time","Error log-s.d.", "Neg.Bin scale","Optimism")

    }

   }else{

    if(unit.scale<0){

     names(aaa$par) <- c("CMP mean","Recruitment Odds","Recruitment Odds Time","Error log-s.d.","Optimism")

    }else{

     names(aaa$par) <- c("Neg.Bin. mean","Recruitment Odds","Recruitment Odds Time","Error log-s.d.","Optimism")

    }

   }

   v=vtrans(aaa$par)

# if(is.psd(-aaa$hessian)){

  if(TRUE){

   g=gtrans(aaa$par)

   asycov = diag(g) %*% robust.inverse(-aaa$hessian) %*% diag(g)

   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))

   asyse <- sqrt(diag(asycov))

   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)

   dimnames(asycor) <- dimnames(asycov)

   names(asyse) <- names(aaa$par)

   ccc <- list(coef=v,iterations=as.numeric(aaa$counts[1]),

               covar=asycov,se=asyse,asycor=asycor,

               df=length(network.size),loglik=aaa$value,

               loglik.null=-length(network.size)*(log(max(network.size))+log(max.coupons+1)))

  }else{

   ccc <- list(theta=vtrans(aaa$par))

  }

 }else{

  ccc <- list(theta=rep(NA,length=6))

 }

 class(ccc) <- "me"

 ccc

}

memle.nb <-function(network.size,num.recruits,recruit.time,max.coupons=3,nb.scale=NULL,cutoff=0,cutabove=1000,

          guess=c(3,0,0,0.6,1,0),

          method="BFGS", hessian=TRUE, K=max(100,network.size), verbose=FALSE){

 logit <- function(p){log(p/(1-p))}

 if(is.null(nb.scale)){

   nb.scale <- -1

   vtrans <- function(l){c((exp(l[1])+1)/exp(l[1]),l[2],l[3],exp(l[4]),exp(l[5])+1,exp(l[6]))}

   ltrans <- function(v){c(-log(v[1]-1),v[2],v[3],log(v[4]),log(v[5]-1),log(v[6]))}

   gtrans <- function(l){c(-exp(-l[1]),1,1,exp(l[4]),exp(l[5]),exp(l[6]))}

 }else{

   nb.scale <- 1

   vtrans <- function(l){c((exp(l[1])+1)/exp(l[1]),l[2],l[3],exp(l[4]),exp(l[5]))}

   ltrans <- function(v){c(-log(v[1]-1),v[2],v[3],log(v[4]),log(v[5]))}

   gtrans <- function(l){c(-exp(-l[1]),1,1,exp(l[4]),exp(l[5]))}

 }

 n <- length(network.size)

#

 if(sum(network.size>=cutoff & network.size <= cutabove) > 0){

  aaa <- optim(par=ltrans(guess),fn=llnbme,

   method=method,

   hessian=hessian,control=list(fnscale=-10),

   n=n,

   recruit.time=recruit.time,

   num.recruits=num.recruits,network.size=network.size,

   max.coupons=max.coupons,K=K,

   nb.scale=nb.scale, verbose=verbose)

   if(nb.scale<0){

     names(aaa$par) <- c("Neg.Bin. mean","Recruitment Odds","Recruitment Odds Time","Error log-s.d.", "Neg.Bin scale","Optimism")

   }else{

     names(aaa$par) <- c("Neg.Bin. mean","Recruitment Odds","Recruitment Odds Time","Error log-s.d.","Optimism")

   }

   v=vtrans(aaa$par)

# if(is.psd(-aaa$hessian)){

  if(TRUE){

   g=gtrans(aaa$par)

   asycov = diag(g) %*% robust.inverse(-aaa$hessian) %*% diag(g)

   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))

   asyse <- sqrt(diag(asycov))

   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)

   dimnames(asycor) <- dimnames(asycov)

   names(asyse) <- names(aaa$par)

   ccc <- list(coef=v,iterations=as.numeric(aaa$counts[1]),

               covar=asycov,se=asyse,asycor=asycor,

               df=length(network.size),loglik=aaa$value,

               loglik.null=-length(network.size)*(log(max(network.size))+log(max.coupons+1)))

  }else{

   ccc <- list(theta=vtrans(aaa$par))

  }

 }else{

  ccc <- list(theta=rep(NA,length=6))

 }

 class(ccc) <- "me"

 ccc

}

llnbme <- function(v,n,network.size,num.recruits,recruit.time,

                  max.coupons,K=100,unit.scale=-1,verbose=FALSE){

 Cret <- .C("gllnbmeC",

           v=as.double(v),

           n=as.integer(n),

           srd=as.double(network.size),

           numrec=as.double(num.recruits),

           rectime=as.double(recruit.time),

           maxcoupons=as.integer(max.coupons),

           K=as.integer(K),

           nb.scale=as.double(unit.scale),

           llik=as.double(0),

           verbose=as.integer(verbose), PACKAGE="RDS")

 out=Cret$llik

 if(verbose){cat(sprintf("llik=%f\n", out))}

 if(is.infinite(out)){out <- -10^10}

 if(is.na(out)){out <- -10^10}

 out

}



dnbmepdf <- function(v,network.size,num.recruits,recruit.time,max.coupons=3,K=100,nb.scale=NULL,verbose=FALSE){

 if(is.null(nb.scale)){

   nb.scale <- -1

 }else{

   nb.scale <- 1

 }

 n <- length(network.size)

 network.size[is.na(network.size)] <- -1

 Cret <- .C("gnbmepdfC",

           v=as.double(v),

           n=as.integer(n),

           srd=as.double(network.size),

           numrec=as.double(num.recruits),

           rectime=as.double(recruit.time),

           maxcoupons=as.integer(max.coupons),

           K=as.integer(K),

           nb.scale=as.double(nb.scale),

           pdf=double(n*K),

           verbose=as.integer(verbose), PACKAGE="RDS")



  return(matrix(Cret$pdf,ncol=n,nrow=K))

}

llcmpme <- function(v,n,network.size,num.recruits,recruit.time,

                  max.coupons,K=100,unit.scale=-1,verbose=FALSE){

 Cret <- .C("gllcmpmeC",

           v=as.double(v),

           n=as.integer(n),

           srd=as.double(network.size),

           numrec=as.double(num.recruits),

           rectime=as.double(recruit.time),

           maxcoupons=as.integer(max.coupons),

           K=as.integer(K),

           cmp.scale=as.double(unit.scale),

           llik=as.double(0),

           verbose=as.integer(verbose), PACKAGE="RDS")

 out=Cret$llik

 if(verbose){cat(sprintf("llik=%f\n", out))}

 if(is.infinite(out)){out <- -10^10}

 if(is.na(out)){out <- -10^10}

 out

}



dmepdf <- function(v,network.size,num.recruits,recruit.time,max.coupons=3,K=100,unit.scale=NULL,unit.model="nbinom",verbose=FALSE){

 if(is.null(unit.scale)){

   unit.scale <- -1

 }else{

   unit.scale <- 1

 }

 n <- length(network.size)

 network.size[is.na(network.size)] <- -1

 if(unit.model=="cmp"){

   Cret <- .C("gcmpmepdfC",

             v=as.double(v),

             n=as.integer(n),

             srd=as.double(network.size),

             numrec=as.double(num.recruits),

             rectime=as.double(recruit.time),

             maxcoupons=as.integer(max.coupons),

             K=as.integer(K),

             cmp.scale=as.double(unit.scale),

             pdf=double(n*K),

             verbose=as.integer(verbose), PACKAGE="RDS")

 }else{

   Cret <- .C("gnbmepdfC",

             v=as.double(v),

             n=as.integer(n),

             srd=as.double(network.size),

             numrec=as.double(num.recruits),

             rectime=as.double(recruit.time),

             maxcoupons=as.integer(max.coupons),

             K=as.integer(K),

             nb.scale=as.double(unit.scale),

             pdf=double(n*K),

             verbose=as.integer(verbose), PACKAGE="RDS")

 }



  return(matrix(Cret$pdf,ncol=n,nrow=K))

}

