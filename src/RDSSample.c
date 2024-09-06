//final coding tasks:

//How to set clock for each alter
//clock for every edge in the network 
//possibly have an array of alters that I sort by time. that would replace the CPickNeighbor funciton.
//one record for each edge.  respondent/alter/clock time.  add multiple entries for each edge [do I mean node?].  decrement coupons.


//ok, maybe this isn't that hard--i.e. I can use the same record even



//once it is done:
//test on fauxmadrona.network: 1000 nodes. RDSSample.  
//fauxmadrona is a rds sample of size 500.

//send it to Mark and see if it works

//latest task: RDSEdge only has "timeIndex" now


/*
 *File: RDSSample.c
 *------------------
 *This file implements the RDSSample.h interface.
 *This interface simulates "respondent driven sampling." It is intended to be called from R and to replace, for the purposes of
 *greater speed, similar (though not identical) code in R.
 * Respondent driven sampling is used to collect samples from hard-to-reach populations, i.e. vulnerable, stigmatized, or
 * otherwise hidden groups.  They may include, for instance, men who have sex with men; sex workers; or right wing extremists.
 * (The first two groups are often labeled "vulnerable" in the bioethics literature.)
 * The code simulates the following sampling procedure: (1) an initial "seed" sample is created by handing out a number of coupons to
 * a few members of the target population.  (2) The sample is created and enlarged by the following two repeated events: (a) unsampled
 * person receives coupons from already sampled person; (b)a person with coupon goes to a survey center, takes a survey and is given
 * coupons of their own to hand out. This two-part procedure continues until the targetted number of surveys have been filled out.  (There may
 * still be people with coupons floating around at this point, but they will not be included in the sample.)
 * Currently the simulation works by assigning a random time to each of the two key events, the person's to the center and passing a ticket
 * to another member of the network.  It does so by using the present time plus an amount of time drawn from an exponential 
 * distribution.
 */


#include "RDSSample.h"
#include <time.h>

/*Global Variables*/
int event;


/*3 column matrix
 recruiter ID
 person recruited
 timestamp*/
/*later I should make it possible to have many tickets; right now there is only one*/

/*this definition is used in the function "select," taken from Numerical Algorithms in C*/



/*
 * Function: CRDSSample
 * ------------------------
 * This is the main function in the RDSSample.h interface.  It is also the function that is called from
 * the R code. The arguments of the function each correspond to R objects.
 * Step 1.   The code first uses the arguments tails, heads, dnedges, dn, dflag, and bipartite to create a network
 * object.
 * Step 2.   Then memory is allocated for a number of pointers that will be passed to the remaining functions in the
 * code.  These include, first, three arrays:  mycandidateEdges, an array for the RDSEdges objects that represent
 * people who have received tickets; myUsed, an array that tracks whether or not a member of the population has been
 * sampled; and myCoupons, an array that tracks the number of coupons that each respondent holds.
 * Then they included ""mycandidateIndex," a variable that keeps track of the number of elements in the array 
 * mycandidateEdges; "mychosenIndex," a variable that keeps track of the number of elements in the array
 *  ;
 */


void CRDSSample(int *tails, int *heads, int *dnedges, int *dn, int *dflag, int *bipartite,
                double *seeddist, 
                int *recruitedSample, 
                int *recruitedSampleAll, 
                int *recruiters, 
                int *recruitersAll,
                double *recruitTimes,
                double *recruitTimesAll,
                char **logfile, int *numCoupons, int *seedssim, int *seeds, int *seedSize, int *sampleSize, int *nsims, int *logResults, int *maxLog, int *events){
  
  
  RDSEdge *mycandidateEdges;
  Network nw;
  Network *my_nwp;
  //  double now;
  int *myUsed;
  int *myCoupons;
  int *mycandidateIndex;
  int *mychosenIndex;
  int i, j;
  //this is 1 if we need to re-seed and 0 otherwise
  
  int *myFileIndex;
  
  GetRNGstate();
  //  clock_t start = clock();
  event = 0;
  myFileIndex = (int *) malloc(sizeof(int));
  *myFileIndex = 0;
  
  /*create the network object*/
  Vertex n_nodes = (Vertex) *dn; 
  Edge n_edges = (Edge) *dnedges;
  int directed_flag = *dflag;
  Vertex bip = (Vertex) *bipartite;
  
  
  nw = NetworkInitialize(tails, heads, n_edges,
                         n_nodes, directed_flag, bip, 0);
  my_nwp = &nw;
  
  mycandidateIndex = (int *) malloc(sizeof(int));
  mychosenIndex = (int *) malloc(sizeof(int));
  mycandidateEdges = (RDSEdge *) malloc((n_nodes + n_edges) * sizeof(RDSEdge)); //this needs to change to edges //huh??
  myUsed = (int *) malloc((int) n_nodes * sizeof(int));
  myCoupons = (int *) malloc((int) n_nodes * sizeof(int));
  
  //the function can now carry out more than one simulation
  for(i=0; i < *nsims; i++){
    for(j = 0; j < *seedSize; j++){
      seedssim[j] = seeds[j+i*(*seedSize)];
      //          Rprintf("%d seed %d\n",j,seedssim[j] );
    }
    //  Rprintf("%d of nsims %d; seed %d\n",i,*nsims,seedssim[0] );
    //event = 0;
    initializeSample(mycandidateEdges, recruitedSample, recruiters, recruitTimes, seedssim, myUsed, myCoupons, n_nodes, n_edges, sampleSize, seedSize, logfile, myFileIndex, logResults, maxLog);
    
    //*mychosenIndex = *seedSize;
    *mychosenIndex = 0;
    *mycandidateIndex = (n_nodes + n_edges) - *seedSize;
    while((*mychosenIndex < *sampleSize)){
      /*TicketEventSimple(my_nwp, mycandidateEdges, seeddist, recruitedSample, recruiters, recruitTimes,sampleSize, mycandidateIndex, mychosenIndex, logfile, numCoupons, myUsed, myCoupons, n_nodes, n_edges, myFileIndex, logResults, maxLog, (int) n_nodes);*/
      
      TicketEvent(my_nwp, mycandidateEdges, seeddist, recruitedSample, recruiters, recruitTimes,sampleSize, mycandidateIndex, mychosenIndex, logfile, numCoupons, myUsed, myCoupons, n_nodes, n_edges, myFileIndex, logResults, maxLog, (int) n_nodes);
      
    } 
    //  Rprintf("%d of nsims %d; needed %d\n",i,*nsims,j);
    //multiple simulations requires a larger array
    for(j=0; j < *sampleSize; j++){
      recruitedSampleAll[j+i*(*sampleSize)] = recruitedSample[j];
    }
    for(j=0; j < *sampleSize; j++){
      recruitersAll[j+i*(*sampleSize)] = recruiters[j];   
    }
    for(j=0; j < *sampleSize; j++){
      recruitTimesAll[j+i*(*sampleSize)] = recruitTimes[j];
    }
    
  }
  *events = event;
  PutRNGstate();
  free(mycandidateEdges);
  free(mycandidateIndex);
  free(mychosenIndex);
  free(myUsed);
  free(myFileIndex);
  free(myCoupons);
}


/*
 * Function: initializeSample
 * --------------------------
 * This function initializes a number of arrays that will be used to hold the results of the simulation. 
 * These are "candidateEdges" which holds representations of people who have been given tickets but have
 * not yet used them; "recruitedSample," which holds representations of people who have returned their tickets
 * to the center, taken the survey, and been given tickets to recruit others; and used, an array which keeps track of
 * what members of the population have already been recruited.  The seeds of the RDS process are passed into the function
 * and first entered into the array of population members recruited.  The recruiters and recruit times of the seeds are
 * set to 0 by convention.  Next the seeds are put in the array CandidateEdges to indicate they need to complete the survey.
 */

void initializeSample(RDSEdge * candidateEdges, int *recruitedSample, int *recruiters, double *recruitTimes, int * seeds, int * used, int * coupons, int num_nodes, int num_edges, int *sampleSize, int *seedSize, char **logfile, int* fileIndex, int *logResults, int *maxLog){
  
  int usedSeed, i;
  //would be nice to check that memory is already allocated
  
  //Initialize arrays candidateEdges for the sampling process. 
  //candidateEdges is bigger so change i
  for(i = 0; i < (num_nodes + num_edges); i++){
    candidateEdges[i].Recruited = 0;
    candidateEdges[i].Recruiter = 0;
    candidateEdges[i].timeIndex = 0;
    candidateEdges[i].Survey = 0;
  }
  
  for(i = 0; i < (int) num_nodes; i++){
    used[i] = 0;
  }
  
  for(i = 0; i < (int) num_nodes; i++){
    coupons[i] = 0;
  }
  
  
  for(i = 0; i < *seedSize; i++){
    //add seeds to the array of chosen edges (removed)
    //add seeds to the array of candidate edges
    candidateEdges[(num_nodes + num_edges) - (*seedSize) + i].Recruited = seeds[i];
    candidateEdges[(num_nodes + num_edges) - (*seedSize) + i].timeIndex = GetTime();
    candidateEdges[(num_nodes + num_edges) - (*seedSize) + i].Survey = 1;
    //log event
    if((*logResults) == 1){
      if(*fileIndex == *maxLog){
        snprintf(logfile[(*maxLog)-1], 58, "STOP");
        (*logResults) = 0;
      } else {
        event++;
        snprintf(logfile[*fileIndex], 58,  "%d 0 %d seed 0 0", event, seeds[i]);
        (*fileIndex)++;
      }
    }
  }
  //indicate seeds have been given tickets
  for(i = 0; i < *seedSize; i++){
    usedSeed = seeds[i];
    used[usedSeed-1] = 1;
  }
  
}

/*
 * Function: GetTime
 * ------------------
 *
 */

double GetTime(void){
  double my_result;
  
  my_result = exp_rand();
  return (my_result);
}




/*test again on flobusiness; make sure candidateIndex is not out of bounds; if it is, end the sampling process*/

/*
 * Function: TicketEvent
 * ---------------------
 * This function represents an "event" in the simulation in which a (potential) member of the sample does something
 * with a ticket. First the data structure that holds candidate members of the sample is re-sorted using the function
 * qsort. Then one of two subsidiary functions is called: either CompleteSurvey, which represents a ticket-holder going
 * to the ticket center and completing a survey for the purposes of data collection;  or Recruit One, which represents a ticket-holder
 * who has taken the survey handing a ticket to a new recruit.
 *
 */

//re-writing this function to include re-seeding

void TicketEvent(Network *nwp, RDSEdge *candidateEdges, double *seeddist, int *recruitedSample, int *recruiters, double *recruitTimes, int *sampleSize, int *candidateIndex, int *chosenIndex, char **logfile, int *nCoupons, int *Used, int *Coupons, int num_nodes, int num_edges, int* fileIndex, int *logResults, int *maxLog, int my_n_nodes){
  
  int edgeSize, earliest;
  double my_time;
  
  edgeSize = num_edges + num_nodes;
  /*Sort the array of recruiters by the time they will recruit*/
  earliest = FindEarliest(candidateEdges, candidateIndex, edgeSize);
  //qsort(candidateEdges, edgeSize, sizeof(RDSEdge), CompareNodes);
  
  if(candidateEdges[earliest].Survey == 1){
    CompleteSurvey(nwp, candidateEdges, recruitedSample, recruiters, recruitTimes, earliest, candidateIndex, chosenIndex, logfile, nCoupons, Coupons, fileIndex,logResults, maxLog);
  } else {
    RecruitOne(candidateEdges, sampleSize, candidateIndex, earliest, logfile, Used, Coupons, fileIndex, logResults, maxLog);
  }
  
  //last part of code indicates whether reseeding needs to happen   
  if((*candidateIndex) > ((num_nodes + num_edges) - 1)){
    my_time = candidateEdges[(*candidateIndex)-1].timeIndex;
    Reseed(candidateEdges, seeddist, recruitedSample, recruiters, recruitTimes, candidateIndex, chosenIndex, logfile, Used, Coupons, my_n_nodes, my_time, fileIndex, logResults, maxLog);
  } 
}

/*
 * Function: Reseed
 * -----------------
 * If there are no more potential recruits from the initial seeds, it may still be possible to generate a new group of seeds and
 * recruit from their links.  This function is used to carry out that task.
 */

void Reseed(RDSEdge *candidateEdges, double *seeddist, int *recruitedSample, int *recruiters, double *recruitTimes, int *candidateIndex, int *chosenIndex, char **logfile, int *Used, int* Coupons, int num_nodes, double time, int* fileIndex, int *logResults, int *maxLog){
  
  //where is the time coming from?
  int n_free_nodes, new_node, i, free_node_i;
  int * new_free_nodes;
  double * cumsum;
  double temp;
  
  n_free_nodes = 0;
  for(i = 0; i < num_nodes ; i++){
    if(Used[i] == 0) n_free_nodes++;
  }
  new_free_nodes = (int *) malloc(n_free_nodes * sizeof(int));
  cumsum = (double *)malloc(n_free_nodes * sizeof(double));
  cumsum[0]=0.0;
  free_node_i = 0;
  for(i = 0; i < num_nodes; i++){
    if(Used[i] == 0){
      new_free_nodes[free_node_i] = i+1;
      if(free_node_i>0){
        cumsum[free_node_i]=cumsum[free_node_i-1]+seeddist[i];
      }else{
        cumsum[0]=seeddist[0];
      }
      free_node_i++;
    }
  }
  //choose new node
  temp = cumsum[n_free_nodes-1] * unif_rand();
  for (i=1; i<=n_free_nodes; i++){
    if(temp <= cumsum[i-1]) break;
  }
  //  Rprintf("i=%d n_free_nodes %d; max %f\n",i,n_free_nodes,cumsum[n_free_nodes-1]);
  //  i = 1 + unif_rand() * (n_free_nodes);  //make sure this is the right code
  //Rprintf("i=%d n_free_nodes %d; max \n",i,n_free_nodes);
  new_node = new_free_nodes[i-1];
  (*candidateIndex)--; 
  //code to put seeds in recruited arrays has been removed from here
  //where are the new seeds marked as used?
  Used[new_node - 1] = 1; 
  //add new node to candidateEdges
  candidateEdges[*candidateIndex].Recruiter = 0;
  candidateEdges[*candidateIndex].Recruited = new_node;
  candidateEdges[*candidateIndex].timeIndex = time + GetTime();
  candidateEdges[*candidateIndex].Survey = 1; 
  //log event
  if((*logResults) == 1){
    if(*fileIndex == *maxLog){
      snprintf(logfile[(*maxLog)-1], 58, "STOP");
      (*logResults) = 0;
    } else {
      event++;
      snprintf(logfile[*fileIndex], 58, "%d %f %d reseed 0 0", event, time, new_node); 
      (*fileIndex)++;
    }
  }
  free(new_free_nodes);
  free(cumsum);
}


/*
 * Function: Complete Survey
 * ---------------------------
 * This function is the first of two functions that can be called from the "TicketEvent" wrapper function.  (The second is 
 * RecruitOne.) This function represents the event of a person who has been given the ticket returning to the center
 * to take the survey.  
 * The event is simulated first by removing the [latest??] member of the array candidateEdges.  (The place in the array 
 * is marked by the variable "candidateIndex.") To simulate giving coupons to the person, the number of coupons is added
 * to the array Coupons.  An instance of type RDSEdge is added to the array CandidateEdges for each member of the network
 * to which the person is linked.  A time is randomly assigned at this point. No attempt is made to check whether these "neighbors" have
 * already been recruited; that
 * task is currently delegated to the function RecruitOne. The survey flag associated with the RDSEdge representing the individual
 * is then set to 0, indicating that the survey has been taken.
 
 * Question: would it be better to put members of the sample into the array at this stage? all don't put seeds in during the InitializeSample
 * portion?
 */


/*
 * new function: FindEarliest
 */

int FindEarliest(RDSEdge *candidateEdges, int *candidateIndex, int nnodes){
  int i, earliest;
  double earliest_time, temp_time;
  
  earliest = *candidateIndex;
  earliest_time = candidateEdges[*candidateIndex].timeIndex;
  for(i = (*candidateIndex); i < nnodes; i++){
    temp_time = candidateEdges[i].timeIndex;
    if(temp_time < earliest_time){
      earliest_time = temp_time;
      earliest = i;
    }    
  }
  return(earliest);
}



void CompleteSurvey(Network *nwp, RDSEdge *candidateEdges, int *recruitedSample, int *recruiters, double *recruitTimes, int earliest, int *candidateIndex,  int *chosenIndex, char**logfile, int *nCoupons, int *Coupons, int* fileIndex, int *logResults, int *maxLog){
  
  //add later, checks for directed network
  double now;
  int n_R_nbrs, i;
  Vertex recruiter, v, parent;
  Edge e;
  Vertex *neighbors;
  
  
  // VIP
  // recruit the earliest person who is given a coupon
  now = candidateEdges[earliest].timeIndex;
  recruiter = candidateEdges[earliest].Recruited;
  parent = candidateEdges[earliest].Recruiter;
  
  
  /*add the edge to the list of edges chosen*/
  
  recruitedSample[*chosenIndex] = recruiter;
  recruiters[*chosenIndex] = parent;
  recruitTimes[*chosenIndex] = now;
  
  
  /*indicate that the sample size has grown by one*/
  (*chosenIndex)++;
  
  
  
  //give coupons to the recruiter
  //do I need some kind of check on nCoupons to make sure it's not 0 or a gazillion
  Coupons[recruiter - 1] = nCoupons[recruiter - 1];
  
  /*add an edge for each of the recruiter's neighbors*/
  //log event
  if((*logResults) == 1){
    if(*fileIndex == *maxLog){
      snprintf(logfile[(*maxLog)-1], 58, "STOP");
      (*logResults) = 0;
    } else {
      event++;
      snprintf(logfile[*fileIndex], 58, "%d %f %d interview 0 0", event, now, recruiter);
      //snprintf(logfile[*fileIndex], 58, "whatever");
      // snprintf(logfile[*fileIndex], 58, "%d", event);
      (*fileIndex)++;
    }
  }
  // n_R_nbrs is the number of neighbors of the recruiter
  n_R_nbrs = nwp->indegree[recruiter]+nwp->outdegree[recruiter];
  if(n_R_nbrs == 0){
    if((*logResults) == 1){
      if(*fileIndex == *maxLog){
        snprintf(logfile[(*maxLog)-1], 58, "STOP");
        (*logResults) = 0;
      } else {
        event++;
        snprintf(logfile[*fileIndex], 58, "%d %f %d no_neighbors 0 0", event, now, recruiter);
        (*fileIndex)++;
      }
    }
    //get rid of the edge i just treated and shrink array by 1 by swapping 
    candidateEdges[earliest].Recruited = candidateEdges[*candidateIndex].Recruited;
    candidateEdges[earliest].Recruiter = candidateEdges[*candidateIndex].Recruiter = recruiter;
    candidateEdges[earliest].timeIndex = candidateEdges[*candidateIndex].timeIndex;
    candidateEdges[earliest].Survey =  candidateEdges[*candidateIndex].Survey ;
    (*candidateIndex)++; 
  } else {
    /*identify the neighbors of the recruiter*/
    neighbors = (Vertex *) malloc(n_R_nbrs*sizeof(Vertex));
    i = 0;
    STEP_THROUGH_OUTEDGES(recruiter,e,v) {
      neighbors[i] = v;
      i++;
    } 
    STEP_THROUGH_INEDGES(recruiter, e, v){
      neighbors[i] = v;
      i++;
    }
    
    //note: looks like the RDSEdge representing visit to the survey center is replaced entirely
    //
    // recruit the first neighbor of the recruiter
    candidateEdges[earliest].Recruited = neighbors[0];
    candidateEdges[earliest].Recruiter = recruiter;
    candidateEdges[earliest].timeIndex = now + GetTime();
    candidateEdges[earliest].Survey = 0;
    for(i = 1; i < n_R_nbrs; i++){
      (*candidateIndex)--;
      candidateEdges[*candidateIndex].Recruited = neighbors[i];
      candidateEdges[*candidateIndex].Recruiter = recruiter;
      candidateEdges[*candidateIndex].timeIndex = now + GetTime();
      candidateEdges[*candidateIndex].Survey = 0;
    }
    free(neighbors);
  }
}

/*
 * Function: RecruitOne
 * -----------------------
 * This function is the second of two functions that can be called from the "TicketEvent" wrapper function.  (The first is 
 * CompleteSurvey.) This function represents the event that follows AFTER a person who has returned his/her ticket to the center, 
 * taken the survey, and picked up the next round of tickets.  The RecruitOne event consists of that person attempting to hand
 * one of the new tickets to someone else in the social network. First, the function checks whether this attempt will be a failure.
 * There are two possible reasons for failure: first, the recruiter is out of coupons; second, the person in the network has already
 * been recruited by someone else.
 * left.  If neither condition holds, the person recruited (with recruiter and time information) are added to the three arrays representing
 * the sample. Then the person becomes a new recruiter in turn.  They is represented by giving the person a time to go to the center, take
 * the survey and get tickets.
 * Note: is this really the most natural time to check for "no more coupons"?
 */

void RecruitOne(RDSEdge *candidateEdges, int *sampleSize, int *candidateIndex, int earliest, char**logfile, int *Used, int *Coupons, int* fileIndex, int *logResults, int *maxLog){
  
  Vertex recruiter, recruited;
  double now;
  
  /*Get the next recruiter*/     
  recruited = candidateEdges[earliest].Recruited;
  recruiter = candidateEdges[earliest].Recruiter;
  now = candidateEdges[earliest].timeIndex;
  
  //need to rewrite this
  //  if(Coupons[recruiter - 1] == 0) {
  if(Used[recruited - 1] == 1) {
    //event++;
    //snprintf(logfile[*fileIndex], 58, "%d %f %d no_coupons %d 0", event, now, recruiter, recruited);
    //(*fileIndex)++;
    recruited = 0;
  }
  if(Coupons[recruiter - 1] == 0){
    //event++;
    //snprintf(logfile[*fileIndex], 58, "%d %f %d node_used %d 0", event, now, recruiter, recruited);
    //(*fileIndex)++;
    recruited = 0;
  }
  if(recruited != 0){
    
    /*now the node that was just recruited gets a time to go take the survey and get its tickets*/
    /*replaces old node*/
    candidateEdges[earliest].Recruiter = recruiter; //changed from chosenIndex
    candidateEdges[earliest].Recruited = recruited; 
    candidateEdges[earliest].timeIndex = now + GetTime();
    candidateEdges[earliest].Survey = 1;
    //log event
    if((*logResults) == 1){
      if(*fileIndex == *maxLog){
        snprintf(logfile[(*maxLog)-1], 58, "STOP");
        (*logResults) = 0;
        //            } else {
        //                event++;
        //                snprintf(logfile[*fileIndex], 58, "%d %f %d exchange %d %d", event, now, recruiter, recruited, Coupons[recruiter - 1]);
        //                (*fileIndex)++;
      }
    }
    /*indicate that the recruiter node has one less coupon*/
    Coupons[recruiter - 1]--;
    /*indicate that the recruited has been given a ticket*/
    Used[recruited - 1] = 1; 
  } else {
    //If the candidate edge is unsuitable, the candidate index increases by 1
    //first swap old and new
    candidateEdges[earliest].Recruited = candidateEdges[*candidateIndex].Recruited;
    candidateEdges[earliest].Recruiter = candidateEdges[*candidateIndex].Recruiter;
    candidateEdges[earliest].timeIndex = candidateEdges[*candidateIndex].timeIndex;
    candidateEdges[earliest].Survey =  candidateEdges[*candidateIndex].Survey ;
    (*candidateIndex)++; 
  }
  /*whereas this needs to be simplified*/
}


