/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "getinclCstacked.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void getinclCstacked (
    int *nbyclass,
    double *size, 
    double *props2,
    int *offby,
    int *N, 
    int *K, 
    int *n, 
    int *samplesize,
    int *Nk
) {
  int i, ni, Ki, isamp, isamplesize;
  int j, offbi, Ni;
  double rU, temp;
  // nbyclass = the number of members of class k=1,...,K
  // size = the size (i.e., degree) of the kth class k=1,...,K
  // K = number of classes (of degrees)
  // n = RDS sample size
  // samplesize = number of w.o.replacement samples to take
  // Nk = output: the total number of times a member of class k=1:K is sampled.
  
  GetRNGstate();  /* R function enabling uniform RNG */

ni=(*n);
Ki=(*K);
isamplesize=(*samplesize);
// ni = RDS sample size
// Ki = number of classes (of degrees)
// isamplesize = number of w.o.replacement samples to take
offbi=(*offby);
Ni=(*N);

int *perm = (int *) malloc(sizeof(int) * Ki);
int *tperm = (int *) malloc(sizeof(int) * Ki);
double *tsize = (double *) malloc(sizeof(double) * Ki);
int *tnbyclass = (int *) malloc(sizeof(int) * Ki);
int *snbyclass = (int *) malloc(sizeof(int) * Ki);
int *samp = (int *) malloc(sizeof(int) * ni);
double *dif = (double *) malloc(sizeof(double) * Ki);
double *cdif = (double *) malloc(sizeof(double) * Ki);

if(offbi>0){
  for (i=0; i<offbi; i++){
    dif[0]=props2[0]-nbyclass[0]/(1.0*Ni);
    cdif[0] = dif[0];
    if(dif[0]<0.0){dif[0] = 0.0;}
    for (j=1; j<Ki; j++){
      dif[j]=props2[j]-nbyclass[j]/(1.0*Ni);
      if(dif[j]<0.0){dif[j] = 0.0;}
      /* compute cumulative probabilities */
      cdif[j] = cdif[j-1] + dif[j];
    }
    rU = unif_rand()*cdif[Ki-1];
    for (j = 0; j < Ki; j++) {if (rU <= cdif[j]) break;}
    nbyclass[j]++;
  }
}

if(offbi<0){
  for (i=0; i<(-offbi); i++){
    dif[0]=nbyclass[0]/(1.0*Ni) - props2[0];
    if(dif[0]<0.0){dif[0] = 0.0;}
    if(nbyclass[0]==1){dif[0] = 0.0;}
    temp=dif[0];
    cdif[0] = dif[0];
    for (j=1; j<Ki; j++){
      dif[j]=nbyclass[j]/(1.0*Ni) - props2[j];
      if(dif[j]<0.0){dif[j] = 0.0;} if(nbyclass[j]==1){dif[j] = 0.0;}
      temp+=dif[j];
      /* compute cumulative probabilities */
      cdif[j] = cdif[j-1] + dif[j];
    }
    if(temp==0.0){
      dif[0]=1.0-(props2[0]-nbyclass[0]/(1.0*Ni));
      if(nbyclass[0]==1){dif[0] = 0.0;}
      cdif[0] = dif[0];
      for (j=1; j<Ki; j++){
        dif[j]=1.0-(props2[j]-nbyclass[j]/(1.0*Ni));
        if(nbyclass[j]==1){dif[j] = 0.0;}
        /* compute cumulative probabilities */
        cdif[j] = cdif[j-1] + dif[j];
      }
    }
    rU = unif_rand()*cdif[Ki-1];
    for (j = 0; j < Ki; j++) {if (rU <= cdif[j]) break;}
    nbyclass[j]--;
  }
}

temp=0.0;
for (i=0; i<Ki; i++){
  size[i]*=nbyclass[i];
  temp+=size[i];
}
for (i=0; i<Ki; i++){
  size[i]/=temp;
  snbyclass[i]=nbyclass[i];
}

for (i=0; i<Ki; i++){
  Nk[i]=0;
}
/* Record element identities */
for (i = 0; i < Ki; i++)
  perm[i] = i + 1;

/* Sort probabilities into descending order */
/* Order element identities in parallel */
/* perm is the permutation order of the ith element of the pop */
revsort(size, perm, Ki);
/* Order element nbyclass also */
for(i = 0 ; i < Ki ; i++){
  tnbyclass[i]=nbyclass[i];
}
for(i = 0 ; i < Ki ; i++){
  nbyclass[i]=tnbyclass[perm[i]-1];
}
for(isamp = 0 ; isamp < isamplesize ; isamp++){
  /* Draw new sample */
  for(i = 0 ; i < Ki ; i++){
    tnbyclass[i]=nbyclass[i];
    tsize[i]=size[i];
    tperm[i]=perm[i];
  }
  /* Sample ni from population with Ni=sum(nbyclass) elements with the prob */
  /* of the ith pop in tsize[i] (ordered in descending order given */
  /* by the permutation in perm */
  ProbSampleNoReplaceStacked(Ki, tnbyclass, tsize, tperm, ni, samp);
  
  /* Tabulate */
  for(i = 0 ; i < ni ; i++){
    Nk[samp[i]-1]++;
  }
}
PutRNGstate();  /* Disable RNG before returning */
  for (i=0; i<Ki; i++){
    nbyclass[i]=snbyclass[i];
  }
  free(samp);
  free(tsize);
  free(tnbyclass);
  free(tperm);
  free(perm);
}


static void ProbSampleNoReplaceStacked(int n, int *nbyclass, double *p, int *perm, int nans, int *ans)
{
  // n = number of classes (of degrees)
  // nbyclass = the number of members of class k=1,...,K
  // p = class of ith member of the pop i=1,...,N i.e. degree
  // perm = permutation of class in descending order k=1:n
  // nans = RDS sample size
  // ans = sample values drawn in sequential order i=1:nans
  double rT, mass, totalmass;
  int i, j, nby;
  
  
  /* Compute the sample */
  totalmass = 1.0;
  for (i = 0; i < nans; i++) {
    rT = totalmass * unif_rand();
    mass = 0.0;
    for (j = 0; j < n; j++) {
      mass += p[j];
      if (rT <= mass)
        break;
    }
    ans[i] = perm[j];
    /* update the reduced probabilities */
    totalmass -= (p[j] / nbyclass[j]);
    p[j] *= (1.0-1.0/nbyclass[j]);
    nbyclass[j]--;
    if(j < n - 1 && p[j] < p[j+1]){
      perm[j] = perm[j+1];
      perm[j+1] = ans[i];
      mass = p[j];
      p[j] = p[j+1];
      p[j+1] = mass;
      nby = nbyclass[j];
      nbyclass[j] = nbyclass[j+1];
      nbyclass[j+1] = nby;
      
    }
  }
}
