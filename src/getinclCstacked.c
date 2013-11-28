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
            int *K, 
            int *n, 
            int *samplesize,
            int *Nk
		 ) {
	int i, ni, Ki, isamp, isamplesize;
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

	int *perm = (int *) malloc(sizeof(int) * Ki);
	int *tperm = (int *) malloc(sizeof(int) * Ki);
	double *tsize = (double *) malloc(sizeof(double) * Ki);
	int *tnbyclass = (int *) malloc(sizeof(int) * Ki);
	int *samp = (int *) malloc(sizeof(int) * ni);

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
