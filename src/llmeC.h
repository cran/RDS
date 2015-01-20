#ifndef LLMEC_H
#define LLMEC_H

void gllnbmeC (double *v,
            int *n, 
            double *srd, 
            double *numrec, 
            double *rectime,
	    int *maxcoupons,
            int *K,
            double *nbscale, 
            double *llik, 
	    int *verbose
			 );

void gnbmepdfC (double *v,
            int *n, 
            double *srd, 
            double *numrec, 
            double *rectime,
	    int *maxcoupons,
            int *K,
            double *nbscale, 
            double *dpdf, 
	    int *verbose
			 );

void gllcmpmeC (double *v,
            int *n, 
            double *srd, 
            double *numrec, 
            double *rectime,
	    int *maxcoupons,
            int *K,
            double *nbscale, 
            double *llik, 
	    int *verbose
			 );

void gcmpmepdfC (double *v,
            int *n, 
            double *srd, 
            double *numrec, 
            double *rectime,
	    int *maxcoupons,
            int *K,
            double *nbscale, 
            double *dpdf, 
	    int *verbose
			 );

#endif /* LLMEC_H */
