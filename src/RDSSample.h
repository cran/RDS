/*
 *File: RDSSample.h
 *-----------------
 *This file implements the RDSSample.h interface.
 *
 */

#ifndef _RDSSample_h
#define _RDSSample_h



#include <math.h>
#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include "edgetree.h"

#include "changestat.h"
#include <errno.h>


typedef struct {
  Vertex Recruited;
  Vertex Recruiter;
  double timeIndex;    
  int Survey;
} RDSEdge;


/*
 * Function: PrintProblems
 * -----------------------
 *
 */

void PrintProblems(char **logfile, int *maxLog);

/*
 * new function: FindEarliest
 */
int FindEarliest(RDSEdge *candidateEdges, int *candidateIndex, int nnodes);



/*
 * Function: CRDSSample
 * ------------------------
 * This function draws a sample of nodes from a social network.  To do so, Respondent driven sampling
 * is simulated.  The first six arguments of the function are used to construct the social network
 * object.  The next three arguments are arrays used to hold the results of the sampling process.
 * The array recruitedSample contains the actual sample.  The corresponding elements of the array
 * recruiters stores the recruiter, if any, of each member of the sample.  (The seeds are given a
 * recruiter equal to 0.)  Finally, the array recruitTimes stores the time each member of the sample
 * was recruited during the simulation.  The seeds are again assigned a default value of 0.
 * The argument seed stores the seeds of the Respondent driven sampling process.  The argument 
 * seedSize stores the number of seeds.  The argument sampleSize stores the size of the sample to be drawn.
 */

void CRDSSample(int *tails, int *heads, int *dnedges, int *dn, int *dflag, int *bipartite, double *seeddist, int *recruitedSample, int *recruitedSampleAll,int *recruiters,  int *recruitersAll, double *recruitTimes, double *recruitTimesAll, char **logfile, int *numCoupons, int *seedssim, int *seeds, int *seedSize, int *sampleSize, int *nsims, int *logResults, int *maxLog, int *events);

/*
 * Function: Reseed
 * -----------------
 *
 */



void Reseed(RDSEdge *candidateEdges, double *seeddist, int *recruitedSample, int *recruiters, double *recruitTimes, int *candidateIndex, int *chosenIndex, char **logfile, int *Used, int* Coupons, int num_nodes, double time, int* fileIndex, int *logResults, int *maxLog);

/*
 * Function:initializeSample
 * -------------------------
 * The function initializes several arrays employed to store information during the sampling process.
 * The arguments candidateEdges, chosenEdges, and used are the arrays that are to be initialized.  
 * The arguments seeds, num_nodes, sampleSize, and seedSize contain data used for 
 * initialization.  The argument log is a file to which debugging information is directed. The final
 * version of this code should remove the log argument.
 */


void initializeSample(RDSEdge * candidateEdges, int *recruitedSample, int *recruiters, double *recruitTimes, int * seeds, int * used, int * coupons, int num_nodes, int num_edges, int *sampleSize, int *seedSize, char **logfile, int* fileIndex, int *logResults, int *maxLog);

/*
 * Function: GetTime
 * ------------------------
 * This function simulates a time by drawing from the exponential distribution.
 */

double GetTime(void);


/*
 * Function: TicketEvent
 * ------------------------
 * [rewrite all of this:]
 * The key element of the sampling process can be thought of as a recruiting event.  At a particular time,
 * a member of the social network recruits a new member by passing a ticket to a neighbor.  The function
 * RecruitOne simulates this event.  The nwp argument indicates the network.  The event is simulated by
 * changing the elements of the arrays candidateEdges and chosenEdges and the values candidateIndex
 * and chosenIndex.  The arguments sampleSize and Used hold information needed for the recruiting event
 * and the argument logfile is included for debugging purposes.
 */

//just changed (10/26/11).  Returns 1 if re-seeding needs to happen and 0 if no re-seeding is necessary

void TicketEvent(Network *nwp, RDSEdge *candidateEdges, double *seeddist, int *recruitedSample, int *recruiters, double *recruitTimes, int *sampleSize, int *candidateIndex, int *chosenIndex, char **logfile, int *nCoupons, int *Used, int *Coupons, int num_nodes, int num_edges, int* fileIndex, int *logResults, int *maxLog, int my_n_nodes);

void TicketEventSimple(Network *nwp, RDSEdge *candidateEdges, double *seeddist, int *recruitedSample, int *recruiters, double *recruitTimes, int *sampleSize, int *candidateIndex, int *chosenIndex, char **logfile, int *nCoupons, int *Used, int *Coupons, int num_nodes, int num_edges, int* fileIndex, int *logResults, int *maxLog, int my_n_nodes);


/*
 * Function: RecruitOne
 * -----------------------
 *
 */



void RecruitOne(RDSEdge *candidateEdges, int *sampleSize, int *candidateIndex, int earliest, char**logfile, int *Used, int *Coupons, int* fileIndex, int *logResults, int *maxLog);

/*
 * Function: CompleteSurvey
 * ---------------------------
 *
 */


void CompleteSurvey(Network *nwp, RDSEdge *candidateEdges, int *recruitedSample, int *recruiters, double *recruitTimes, int earliest, int *candidateIndex,  int *chosenIndex, char**logfile, int *nCoupons, int *Coupons, int* fileIndex, int *logResults, int *maxLog);

/*
 * Function: CompareNodes
 * ------------------------
 * This function is used to sort Respondent objects according to their timeToRecruit field.
 */

int CompareNodes (const void * a, const void * b);




#endif
