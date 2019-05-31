/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "llmeC.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void gllcmpmeC (double *v,
                int *n, 
                int *srd, 
                double *numrec, 
                double *rectime,
                int *maxcoupons,
                int *K,
                double *cmpscale, 
                double *llik, 
                int *verbose
) {
  int i, iu, ni, Ki, maxc;
  double dmean,beta0,beta1,lnsd,scale,opt;
  int give_log0=0;
  double temp, loglik, sint, u, maxpdf;
  double eps=0.00001;
  
  Ki=(*K);
  ni=(*n);
  maxc=(*maxcoupons);
  
  double *pdf = (double *) malloc(sizeof(double) * (Ki+1));
  double *rtprob = (double *) malloc(sizeof(double) * ni);
  
  dmean=(exp(v[0])+1.0)/exp(v[0]);
  beta0=v[1];
  beta1=v[2];
  lnsd=exp(v[3]);
  if((*cmpscale)<0.0){
    scale=exp(v[4]);
    opt=exp(v[5]);
  }else{
    scale=(*cmpscale);
    opt=exp(v[4]);
  }
  for (i=0; i<ni; i++){
    temp = beta0 + beta1*rectime[i];
    rtprob[i] = exp(temp)/(1.0+exp(temp));
//    Rprintf("%f %f %f %f\n",beta0,beta1,rtprob[i],temp);
  }
  //  for (iu=0; iu<Ki; iu++){
  //    pdf[iu]=exp((iu+1.0)*log(dmean)-scale*lgamma(iu+2.0));
  //    temp+=pdf[iu];
  ////  pdf[iu] = dnbinom_mu(u,scale,dmean,give_log0);
  //  }
  ////dcmp (intu, dmean, scale, Ki, 0.000000001, give_log0, pdf)
  //  for (iu=0; iu<Ki; iu++){
  //    pdf[iu]/=temp;
  //  }
  pdf[0]=0.0;
  pdf[1]=1.0;
  temp=1.0;
  maxpdf=1.0;
  loglik=1.0;
  iu=2;
  while ((loglik > eps*maxpdf) && (iu <= Ki) ){
    u=iu-1.0;
    loglik=exp(u*log(dmean)-scale*lgamma(u+1.0));
    pdf[iu]=loglik;
    if(loglik>maxpdf) maxpdf=loglik;
    temp+=loglik;
    iu++;
//  Rprintf("u %d pdf %f\n",iu,pdf[iu]);
  }
  Ki=iu-1;
  //Rprintf("%f %f %f %f %f\n",dmean,scale,beta0,beta1,temp);
  for (iu=1; iu<=Ki; iu++){
    pdf[iu]/=temp;
    //Rprintf("u %d pdf %f\n",iu,pdf[iu]);
  }
  
  loglik=0.0;
  for (i=0; i<ni; i++){
    sint = 0.0; 
    for (iu=1; iu<=Ki; iu++){
      temp = pdf[iu];
      u = (double)iu;
      if((iu <= maxc)|(numrec[i]<maxc)){
        temp *= dbinom(numrec[i],u,rtprob[i],give_log0);
//        Rprintf("< %f %f %f %f\n",numrec[i],rtprob[i],u,dbinom(numrec[i],u,rtprob[i],give_log0));
      }else{
        temp *= 1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0);
//        Rprintf("< %f\n",1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0));
      }
      //  temp *= dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0);
      if(srd[i]>=0){
        temp *= poilog(srd[i],log(u)-log(opt),lnsd);
      }
      //  Rprintf("dn %f\n", poilog(srd[i],log(u)-log(opt),lnsd));
//      Rprintf("dn %f\n", log(temp));
      if(!isnan(temp)) sint += temp; 
    }
    // if(sint <= 0.0) Rprintf("sint %f\n", sint);
//  Rprintf("log(sint) %f\n", log(sint));
    loglik+=log(sint);
  }
  (*llik)=loglik;
  free(rtprob);
  free(pdf);
}

void gcmpmepdfC (double *v,
                 int *n, 
                 int *srd, 
                 double *numrec, 
                 double *rectime,
                 int *maxcoupons,
                 int *K,
                 double *cmpscale, 
                 double *dpdf, 
                 int *verbose
) {
  int i, iu, ni, Ki, maxc;
  double dmean,beta0,beta1,lnsd,scale,opt;
  int give_log0=0;
  double temp, sint, u, loglik, maxpdf;
  double eps=0.00001;
  
  Ki=(*K);
  ni=(*n);
  maxc=(*maxcoupons);
  
  double *pdf = (double *) malloc(sizeof(double) * (Ki+1));
  //int *intu = (double *) malloc(sizeof(int) * Ki);
  double *rtprob = (double *) malloc(sizeof(double) * ni);
  
  dmean=v[0];
  beta0=v[1];
  beta1=v[2];
  lnsd=v[3];
  if((*cmpscale)<0.0){
    scale=v[4];
    opt=v[5];
  }else{
    scale=(*cmpscale);
    opt=v[4];
  }
  for (i=0; i<ni; i++){
    temp = beta0 + beta1*rectime[i];
    rtprob[i] = exp(temp)/(1.0+exp(temp));
  }
  
  temp=1.0;
  pdf[0]=0.0;
  pdf[1]=1.0;
  maxpdf=1.0;
  loglik=1.0;
  iu=2;
  while ((loglik > eps*maxpdf) && (iu <= Ki) ){
    u=iu-1.0;
    loglik=exp(u*log(dmean)-scale*lgamma(u+1.0));
    pdf[iu]=loglik;
    if(loglik>maxpdf) maxpdf=loglik;
    temp+=loglik;
    iu++;
    //Rprintf("u %d pdf %f\n",iu,pdf[iu]);
  }
  Ki=iu-1;
  for (iu=1; iu<=Ki; iu++){
    pdf[iu]/=temp;
  }
  
  for (i=0; i<ni; i++){
    sint = 0.0; 
    for (iu=1; iu<=Ki; iu++){
      temp = pdf[iu];
      u = (double)iu;
      if((iu <= maxc)|(numrec[i]<maxc)){
        temp *= dbinom(numrec[i],u,rtprob[i],give_log0);
        //  Rprintf("< %f %f %f\n",rtprob[i],u,dbinom(numrec[i],u,rtprob[i],give_log0));
      }else{
        temp *= 1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0);
        //  Rprintf("< %f\n",1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0));
      }
      if(srd[i]>=0){
        //   temp *= dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0);
        temp *= poilog(srd[i],log(u)-log(opt),lnsd);
      }
      //  Rprintf("dn %f\n", dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0));
      sint += temp; 
      dpdf[i*Ki+iu-1] = temp; 
    }
    for (iu=1; iu<=Ki; iu++){
      dpdf[i*Ki+iu-1] /= sint; 
    }
  }
  *K=Ki;
  free(rtprob);
  free(pdf);
}

void gllnbmeC (double *v,
               int *n, 
               int *srd, 
               double *numrec, 
               double *rectime,
               int *maxcoupons,
               int *K,
               double *nbscale, 
               double *llik, 
               int *verbose
) {
  int i, iu, ni, Ki, maxc;
  double dmean,beta0,beta1,lnsd,shape,opt;
  int give_log0=0;
  double temp, loglik, sint, u;
  double eps=0.00001;
  
  Ki=(*K);
  ni=(*n);
  maxc=(*maxcoupons);
  
  double *pdf = (double *) malloc(sizeof(double) * (Ki+1));
  double *rtprob = (double *) malloc(sizeof(double) * ni);
  
  dmean=(exp(v[0])+1.0)/exp(v[0]);
  beta0=v[1];
  beta1=v[2];
  lnsd=exp(v[3]);
  if((*nbscale)<0.0){
    shape=dmean / exp(v[4]);
    opt=exp(v[5]);
  }else{
    shape=(*nbscale);
    opt=exp(v[4]);
  }
  for (i=0; i<ni; i++){
    temp = beta0 + beta1*rectime[i];
    rtprob[i] = exp(temp)/(1.0+exp(temp));
    //  Rprintf("%f %f %f %f\n",beta0,beta1,rtprob[i],temp);
  }
  pdf[0]=0.0;
  temp=0.0;
  iu=1;
  while ((temp < 1.0 - eps) && (iu <= Ki) ){
    u=((double)iu)-1.0;
    loglik = dnbinom_mu(u,shape,dmean,give_log0);
    temp+=loglik;
    pdf[iu] = loglik;
    iu++;
  }
  Ki=iu-1;
  *K=iu-1;
  pdf[Ki]+=1.0-temp;
  //pdf[Ki]+=1.0-temp-dnbinom_mu(0.0,shape,dmean,give_log0);
  //temp=1.0-dnbinom_mu(0.0,shape,dmean,give_log0);
  //for (iu=1; iu<=Ki; iu++){
  //  pdf[iu]/=temp;
  //}
  
  loglik=0.0;
  for (i=0; i<ni; i++){
    sint = 0.0; 
    for (iu=1; iu<=Ki; iu++){
      temp = pdf[iu];
      u = (double)iu;
      if((iu <= maxc)|(numrec[i]<maxc)){
        temp *= dbinom(numrec[i],u,rtprob[i],give_log0);
        //  Rprintf("< %f %f %f\n",rtprob[i],u,dbinom(numrec[i],u,rtprob[i],give_log0));
      }else{
        temp *= 1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0);
        //  Rprintf("< %f\n",1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0));
      }
      //  temp *= dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0);
      if(srd[i]>=0){
        temp *= poilog(srd[i],log(u)-log(opt),lnsd);
      }
      //  Rprintf("pln srd %d %f\n", srd[i], poilog(srd[i],log(u)-log(opt),lnsd));
      //  Rprintf("dn %f\n", dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0));
      sint += temp; 
    }
    loglik+=log(sint);
  }
  (*llik)=loglik;
  free(rtprob);
  free(pdf);
}

void gnbmepdfC (double *v,
                int *n, 
                int *srd, 
                double *numrec, 
                double *rectime,
                int *maxcoupons,
                int *K,
                double *nbscale, 
                double *dpdf, 
                int *verbose
) {
  int i, iu, ni, Ki, maxc;
  double dmean,beta0,beta1,lnsd,shape,opt;
  int give_log0=0;
  double temp, sint, u;
  double zero=0.0, eps=0.00001;
  
  Ki=(*K);
  ni=(*n);
  maxc=(*maxcoupons);
  
  double *pdf = (double *) malloc(sizeof(double) * (Ki+1));
  double *rtprob = (double *) malloc(sizeof(double) * ni);
  
  dmean=v[0];
  beta0=v[1];
  beta1=v[2];
  lnsd=v[3];
  if((*nbscale)<zero){
    shape=dmean / (v[4]-1.0);
    opt=v[5];
  }else{
    shape=(*nbscale);
    opt=v[4];
  }
  for (i=0; i<ni; i++){
    temp = beta0 + beta1*rectime[i];
    rtprob[i] = exp(temp)/(1.0+exp(temp));
  }
  pdf[0]=zero;
  temp=zero;
  iu=1;
  while ((temp < 1.0 - eps) && (iu <= Ki) ){
    u=((double)iu)-1.0;
    sint = dnbinom_mu(u,shape,dmean,give_log0);
    temp+=sint;
    pdf[iu] = sint;
    iu++;
  }
  Ki=iu-1;
  *K=iu-1;
  pdf[Ki]+=1.0-temp;
  //pdf[Ki]+=1.0-temp-dnbinom_mu(0.0,shape,dmean,give_log0);
  //temp=1.0-dnbinom_mu(0.0,shape,dmean,give_log0);
  //for (iu=1; iu<=Ki; iu++){
  //  pdf[iu]/=temp;
  //}
  
  for (i=0; i<ni; i++){
    sint = zero; 
    for (iu=1; iu<=Ki; iu++){
      temp = pdf[iu];
      u = (double)iu;
      if((iu <= maxc)|(numrec[i]<maxc)){
        temp *= dbinom(numrec[i],u,rtprob[i],give_log0);
        //  Rprintf("< %f %f %f\n",rtprob[i],u,dbinom(numrec[i],u,rtprob[i],give_log0));
      }else{
        temp *= 1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0);
        //  Rprintf("< %f\n",1.0-pbinom(maxc-1.0,u,rtprob[i],give_log0,give_log0));
      }
      if(srd[i]>=0){
        //   temp *= dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0);
        temp *= poilog(srd[i],log(u)-log(opt),lnsd);
      }
      //  Rprintf("dn %f\n", dnorm(log(srd[i])-log(u)+log(opt),zero,lnsd,give_log0));
      sint += temp; 
      dpdf[i*Ki+iu-1] = temp; 
    }
    for (iu=1; iu<=Ki; iu++){
      dpdf[i*Ki+iu-1] /= sint; 
    }
  }
  free(rtprob);
  free(pdf);
}
