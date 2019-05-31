/*******************************************************************/
/* Computation of the log-likelihood and marginal posterior of size*/
/*******************************************************************/

#include "bsC.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>

void bsC (
            int *nbyclass,
            int *classesboth,
            int *nrefs,
            double *props2,
            double *tempties,
            double *pis, 
            double *est,
            double *nm,
            int *numsamp,
            int *offby,
            int *K, 
            int *nc, 
            int *g, 
            int *N, 
            int *n, 
            int *n0
		 ) {
	int i, j;
	int isamp, nsamples;
	int numberfrom, nextdis, nextresp, tdeg, thisdis;
	int activnode, countrefs, mm;
	int offbi, Ni, ni, n0i, Ki, gi, ci;
	double rU, temp, totaltmp, den, den2;
	// nbyclass = the number of members of class k=1,...,K
	// size = the size (i.e., degree) of the kth class k=1,...,K
	// K = number of classes (of degrees)
	// n = RDS sample size
	// samplesize = number of w.o.replacement samples to take
	// Nk = output: the total number of times a member of class k=1:K is sampled.

	GetRNGstate();  /* R function enabling uniform RNG */

	nsamples=(*numsamp);
	offbi=(*offby);
	Ki=(*K);
	ci=(*nc);
	gi=(*g);
	Ni=(*N);
	ni=(*n);
	n0i=(*n0);
	// ni = RDS sample size
	// Ki = number of classes (of degrees)
	// isamplesize = number of w.o.replacement samples to take

	int *newnbyclass = (int *) malloc(sizeof(int) * ci);
	double *ntt = (double *) malloc(sizeof(double) * gi*gi);
	double *nprop = (double *) malloc(sizeof(double) * gi);
	double *nprob = (double *) malloc(sizeof(double) * ci);
	double *thistempties = (double *) malloc(sizeof(double) * gi);
	int *ideg = (int *) malloc(sizeof(int) * ci);
	int *idis = (int *) malloc(sizeof(int) * ci);
	int *pclass = (int *) malloc(sizeof(int) * ci);
	double *dif = (double *) malloc(sizeof(double) * ci);
	double *cdif = (double *) malloc(sizeof(double) * ci);
	int *crefs = (int *) malloc(sizeof(int) * ni);
	int *csample = (int *) malloc(sizeof(int) * ni);

//	for (i=0; i<ci; i++){
//	  Rprintf("i %d nbyclass[i] %d\n",i,nbyclass[i]);
//	}
//
	for (j=0; j<(gi*nsamples); j++){
	 est[j]=0.0;
	}
	for (j=0; j<nsamples; j++){
	 nm[j]=0.0;
	}
	for (i=0; i<ci; i++){
	  idis[i]=((classesboth[i]-1)/Ki);
	  ideg[i]=classesboth[i]-Ki*(classesboth[i]/Ki);
	  if(ideg[i]==0){ideg[i]=Ki;}
	}
	for (i=0; i<ni; i++){
	  crefs[i]=0;
	  csample[i]=0;
	}
	for (i=0; i<ni; i++){
	    if(nrefs[i] > 0){crefs[nrefs[i]-1]++;}
	}
	for (i=1; i<ni; i++){
	    crefs[i]+=crefs[i-1];
	}
//	  Rprintf("crefs[0] %d crefs[1] %d crefs[2] %d crefs[3] %d\n",crefs[0], crefs[1], crefs[2], crefs[3]);

	for (isamp=0; isamp<nsamples; isamp++){

	for (j=0; j<ci; j++){
	  newnbyclass[j] = nbyclass[j];
	}
	// ntt is a g x g matrix of ties between response classes
	for (j=0; j<(gi*gi); j++){
	 ntt[j]=tempties[j];
	}

	if(offbi>0){
	  for (i=0; i<offbi; i++){
	    dif[0]=props2[0]-newnbyclass[0]/(1.0*Ni);
	    cdif[0] = dif[0];
	    if(dif[0]<0.0){dif[0] = 0.0;}
	    for (j=1; j<ci; j++){
	      dif[j]=props2[j]-newnbyclass[j]/(1.0*Ni);
	      if(dif[j]<0.0){dif[j] = 0.0;}
	      /* compute cumulative probabilities */
	      cdif[j] = cdif[j-1] + dif[j];
	    }
	    rU = unif_rand()*cdif[ci-1];
	    for (j = 0; j < ci; j++) {if (rU <= cdif[j]) break;}
	    newnbyclass[j]++;
	  }
	}

	if(offbi<0){
	  for (i=0; i<(-offbi); i++){
	    dif[0]=newnbyclass[0]/(1.0*Ni) - props2[0];
	    if(dif[0]<0.0){dif[0] = 0.0;}
	    if(newnbyclass[0]==1){dif[0] = 0.0;}
	    totaltmp=dif[0];
	    cdif[0] = dif[0];
	    for (j=1; j<ci; j++){
	      dif[j]=newnbyclass[j]/(1.0*Ni) - props2[j];
	      if(dif[j]<0.0){dif[j] = 0.0;} if(newnbyclass[j]==1){dif[j] = 0.0;}
	      totaltmp+=dif[j];
	      /* compute cumulative probabilities */
	      cdif[j] = cdif[j-1] + dif[j];
	    }
	    if(totaltmp==0.0){
	     dif[0]=1.0-(props2[0]-newnbyclass[0]/(1.0*Ni));
	     if(newnbyclass[0]==1){dif[0] = 0.0;}
	     cdif[0] = dif[0];
	     for (j=1; j<ci; j++){
	      dif[j]=1.0-(props2[j]-newnbyclass[j]/(1.0*Ni));
	      if(newnbyclass[j]==1){dif[j] = 0.0;}
	      /* compute cumulative probabilities */
	      cdif[j] = cdif[j-1] + dif[j];
	     }
	    }
	    rU = unif_rand()*cdif[ci-1];
	    for (j = 0; j < ci; j++) {if (rU <= cdif[j]) break;}
	    newnbyclass[j]--;
	  }
	}

//	Rprintf("make dif\n");
//	for (i=0; i<ci; i++){
//	  Rprintf("i %d nbyclass[i] %d\n",i,nbyclass[i]);
//	}

	// pclass[i] is the sum of the degrees of those in class i
	// So the class is chosen with probability pclass[i]
	for (i=0; i<ci; i++){
	  pclass[i]=ideg[i]*newnbyclass[i];
	}
//	for (i=0; i<ci; i++){
//	  if(idis[i]<0 | idis[i]>1){Rprintf("Error: i %d idis[i] %d ideg %d\n",i,idis[i],ideg[i]);}
//	}

//	for (i=0; i<ci; i++){
//	  Rprintf("i %d idis[i] %d ideg %d\n",i,idis[i],ideg[i]);
//	}

	// Sample seeds
	for (i=0; i<n0i; i++){
	  cdif[0] = pclass[0];
	  for (j=1; j<ci; j++){
	    /* compute cumulative probabilities */
	    cdif[j] = cdif[j-1] + pclass[j];
	  }
	  rU = unif_rand()*cdif[ci-1];
	  for (j = 0; j < ci; j++) {if (rU <= cdif[j]) break;}
	  csample[i] = j;
	  pclass[j]-=ideg[j];
	}

	for (i=0; i<gi; i++){
	  nprop[i]=0.0;
	}
	for (i=0; i<ci; i++){
//	  Rprintf("i %d ci %d\n",i,ci);
	    nprop[idis[i]]=nprop[idis[i]]+newnbyclass[i];
//	  Rprintf("idis[i] %d nprop[idis[i]] %f nbyclass[i] %d i %d\n",idis[i],nprop[idis[i]],nbyclass[i],i);
//	  if(idis[i]<0 | idis[i]>1){Rprintf("Error: i %d idis[i] %d ideg %d\n",i,idis[i],ideg[i]);}
	  }
	  temp=0.0;
	  for(i = 0 ; i < gi ; i++){
	    temp+=nprop[i];
	  }
	  for(i = 0 ; i < gi ; i++){
	    nprop[i]/=temp;
	  }
	  for(i = 0 ; i < gi ; i++){
	  temp=0.0;
	  for(j = 0 ; j < gi ; j++){
	    temp+=ntt[i+gi*j];
	  }
	  if(temp<=0.0){
	  for(j = 0 ; j < gi ; j++){
	    ntt[i+gi*j] = nprop[j];
	   }
	  }
	  temp=0.0;
	  for(j = 0 ; j < gi ; j++){
	    temp+=ntt[j+gi*i];
	  }
	  if(temp<=0.0){
	   for(j = 0 ; j < gi ; j++){
	    ntt[j+gi*i] = nprop[j];
	   }
	  }
	}

//	
//	Rprintf("make tempties\n");

	rU = unif_rand()*crefs[ni-1];
	for (j = 0; j < ni; j++) {if (rU <= crefs[j]) break;}
	numberfrom = j+1;

//	Rprintf("numberfrom %d\n",numberfrom);

	activnode = 0;
	countrefs = 0;

        for(mm = n0i ; mm < ni ; mm++){ // loop begins!!!  mm is not m in the paper!

//	if(activnode >= mm){Rprintf("Error: activnode >= mm %d %d\n",activnode,mm);}
//	 thisdis is the disease status of the active node
	 thisdis = idis[csample[activnode]];
//	 nextdis<-.Internal(sample(g, 1, FALSE, tempties[thisdis,]))

	 for(i = 0 ; i < gi ; i++){
	  thistempties[i]=ntt[thisdis+gi*i];
	 }
	 for (j=1; j<gi; j++){
	  /* compute cumulative probabilities */
	  thistempties[j]+=thistempties[j-1];
	 }
	 rU = unif_rand()*thistempties[gi-1];
	 for (j = 0; j < gi; j++) {if (rU <= thistempties[j]) break;}
	 // nextdis is the disease status of a random referral for a node of
	 // the same type as the active node
	 nextdis = j;

//	Rprintf("nextdis %d\n",nextdis);

	 temp=0.0;
	 for(i = 0 ; i < ci ; i++){
	  if(idis[i]==nextdis){
	   nprob[i]=pclass[i];
	   temp+=pclass[i];
	  }else{
	   nprob[i]=0.0;
	  }
	 }
	 if(temp==0.0){
	  for(i = 0 ; i < ci ; i++){
	    nprob[i]=pclass[i];
	  }
//	 Rprintf("Ran out of %d\n",nextdis);
	 }

	 cdif[0] = nprob[0];
	 for (j=1; j<ci; j++){
	   /* compute cumulative probabilities */
	   cdif[j] = cdif[j-1] + nprob[j];
	 }
	 rU = unif_rand()*cdif[ci-1];
	 for (j = 0; j < ci; j++) {if (rU <= cdif[j]) break;}
	 nextresp = j;
	 totaltmp=0.0;
	 for (j=0; j<ci; j++){
	   totaltmp+=nprob[j];
	 }
	 pclass[nextresp]-=ideg[nextresp];
	 nextdis = idis[nextresp];
	 tdeg = ideg[nextresp];
	 for(i = 0 ; i < gi ; i++){
	   ntt[i+gi*nextdis]*=(totaltmp-tdeg)/totaltmp;
	 }
//	Rprintf("ntt %f %f %f %f\n",ntt[0],ntt[1],ntt[2],ntt[3]);
	 csample[mm] = nextresp;
	 countrefs++;
//      numberfrom is the number of recruits to get for the current recruiter
	 if((mm<ni)&&(countrefs==numberfrom)){
	  activnode++;                     // move to the next seed (or node)!!! i=i+1
	  countrefs=0;
	  rU = unif_rand()*crefs[ni-1];
	  for (j = 0; j < ni; j++) {if (rU <= crefs[j]) break;}
	  numberfrom = j+1;
	 }
//	Rprintf("make sample %d\n",mm);

//      end mm loop
	}
//	Rprintf("tempties %f %f %f %f\n",ntt[0],ntt[1],ntt[2],ntt[3]);
//	if(est[0]==1017){Rprintf("est %f\n",est[0]);}

	den=0.0;
	den2=0.0;
	for (i=0; i<ni; i++){
	 temp=1.0/pis[csample[i]];
	 den+=temp;
	 den2+=(temp*temp);
	 est[idis[csample[i]]+gi*isamp]+=temp;
	 //	Rprintf("csample %d pix %d pis %f\n",csample[i], csample[i]-UKi*(csample[i]/UKi),pis[csample[i]-UKi*(csample[i]/UKi)]);
//	  if(csample[i]<0 | csample[i]>(ci-1)){Rprintf("Error: i %d csample[i] %d\n",i, csample[i]);}
	}
//	Rprintf("est[0] %f est[1] %f\n",est[0],est[1]);
	for (j=0; j<gi; j++){
	 est[j+gi*isamp]/=den;
//	Rprintf("est %d %f %f\n",j,den,est[j+gi*isamp]);
	}
	nm[isamp]=den*den/den2;

	}

	PutRNGstate();  /* Disable RNG before returning */

	free(newnbyclass);
	free(ntt);
	free(nprop);
	free(nprob);
	free(thistempties);
	free(ideg);
	free(idis);
	free(pclass);
	free(dif);
	free(cdif);
	free(crefs);
	free(csample);
}
