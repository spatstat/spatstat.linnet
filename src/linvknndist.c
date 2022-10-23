#include <R.h>
#include "yesno.h"

/*

  linvknndist.c

  k-th nearest neighbour function at vertices
  (distance from each vertex to the 
  nearest, second nearest, ...  k-th nearest target data point)

  Needs only the sparse representation of the network

  $Revision: 1.8 $  $Date: 2022/10/22 09:53:55 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  ! Data points must be ordered by segment index !

 */

#undef HUH

#define DIST(VERTEX, ORDER) dist[(ORDER) + (VERTEX) * Kmax]
#define WHICH(VERTEX, ORDER) which[(ORDER) + (VERTEX) * Kmax]

#define UPDATE(VERTEX, D, J, EPS)	\
  UpdateKnnList(D, J, \
		dist + (VERTEX) * Kmax, \
		which + (VERTEX) * Kmax, \
		Kmax, \
		EPS)


void linvknndist(
  int *kmax,         /* number of neighbours required */
  /* target data points in local coordinates (ordered by sq) */
  int *nq,
  int *sq,
  double *tq,
  /* network */
  int *nv,           /* number of network vertices */
  int *ns,
  int *from,
  int *to,           /* segments (pairs of vertices) */
  double *seglen,       /* segment lengths */
  double *huge,         /* value taken as infinity */
  double *tol,          /* tolerance for updating distances */
  /* OUTPUT */
  double *dist,         /* distance from each vertex to
			   the nearest, ..., kth nearest data points */
  int *which         /* identifies which data points */
) 
{
  int Nq, Nv, Ns, Kmax, Nout, i, j, k, m;
  int segQj, ivleft, ivright, changed;
  double hugevalue, eps, slen, d, tqj;
  char converged;

  int UpdateKnnList(double d, int j,
		    double *dist, int *which, int Kmax, double eps);

  Kmax = *kmax;
  Nq = *nq;
  Nv = *nv;
  Ns = *ns;
  hugevalue = *huge;
  eps = *tol;

  /* number of values in 'dist' and in 'which' */
  Nout = Nv * Kmax;

#ifdef HUH
  Rprintf("Initialise dist\n");
#endif
  /* initialise to huge value */
  for(i = 0; i < Nout; i++) {
    dist[i] = hugevalue;
    which[i] = -1;
  }


#ifdef HUH
  Rprintf("Run through target points\n");
#endif
  /* assign value to endpoints of segments containing target points */
  for(j = 0; j < Nq; j++) {
    segQj = sq[j];
    tqj = tq[j];
    slen = seglen[segQj];
    ivleft = from[segQj];
    d = slen * tqj;
    UPDATE(ivleft, d, j, (double) 0.0);
    ivright = to[segQj];
    d = slen * (1.0 - tqj);
    UPDATE(ivright, d, j, (double) 0.0);
  }

#ifdef HUH
  Rprintf("Initialised values at vertices:\n");
  Rprintf("\ti\twhich\tdist\n");
  for(i = 0; i < Nv; i++) {
    Rprintf("\t%d", i);
    for(k = 0; k < Kmax; k++) 
      Rprintf(" %d ", WHICH(i, k));
    for(k = 0; k < Kmax; k++) 
      Rprintf(" %lf ", DIST(i, k));
    Rprintf("\n");
  }
#endif

  /* recursively update */
#ifdef HUH
  Rprintf("Recursive update\n");
#endif
  converged = NO;
  while(!converged) {
    converged = YES;
#ifdef HUH
    Rprintf("........... starting new pass ...................... \n");
    Rprintf("Current state:\n");
    Rprintf("\ti\twhich\tdist\n");
    for(i = 0; i < Nv; i++) {
      Rprintf("\t%d", i);
      for(k = 0; k < Kmax; k++) 
	Rprintf(" %d ", WHICH(i, k));
      for(k = 0; k < Kmax; k++) 
	Rprintf(" %lf ", DIST(i, k));
      Rprintf("\n");
    }
#endif
    for(m = 0; m < Ns; m++) {
      ivleft = from[m];
      ivright = to[m];
      slen = seglen[m];

#ifdef HUH
      Rprintf("updating right=%d from left=%d\n", ivright, ivleft);
#endif
      for(k = 0; k < Kmax; k++) {
	changed = UPDATE(ivright, DIST(ivleft, k)+slen, WHICH(ivleft, k), eps);
	converged = converged && !changed;
      }

#ifdef HUH
      Rprintf("updating left=%d from right=%d\n", ivleft, ivright);
#endif
      for(k = 0; k < Kmax; k++) {
	changed = UPDATE(ivleft, DIST(ivright, k)+slen, WHICH(ivright, k), eps);
	converged = converged && !changed;
      }
    }
  }

#ifdef HUH
  Rprintf("Done\nVertex values:\n");
  Rprintf("\ti\twhich\tdist\n");
  for(i = 0; i < Nv; i++) {
    Rprintf("\t%d", i);
    for(k = 0; k < Kmax; k++) 
      Rprintf(" %d ", WHICH(i, k));
    for(k = 0; k < Kmax; k++) 
      Rprintf(" %lf ", DIST(i, k));
    Rprintf("\n");
  }
#endif
}


/* update a list of nearest, second nearest, ..., k-th nearest neighbours */

int UpdateKnnList(
  double d,      /* candidate distance */
  int j,         /* corresponding candidate target point */
  double *dist,  /* pointer to start of vector of length Kmax */
  int *which,    /* pointer to start of vector of length Kmax */
  int Kmax,  
  double eps     /* numerical tolerance, to prevent infinite loops */
) {
  char matched, unsorted, changed;
  int k, Klast, itmp;
  double dtmp, dPlusEps;

  Klast = Kmax - 1;

  dPlusEps = d + eps;
  if(dPlusEps > dist[Klast])
    return(NO);

  changed = NO;

  /* Check whether this data point is already listed as a neighbour */
  matched = NO;
  for(k = 0; k < Kmax; k++) {
    if(which[k] == j) {
      matched = YES;
#ifdef HUH
      Rprintf("\tMatch: which[%d] = %d\n", k, j);
#endif
      if(dPlusEps <= dist[k]) {
	changed = YES;
#ifdef HUH
	Rprintf("\t\tUpdated distance from %lf to %lf\n", dist[k], d);
#endif
	dist[k] = d;
      }
      break;
    }
  }
  if(!matched) {
#ifdef HUH
    Rprintf("\tNo match with current list\n");
    Rprintf("\t\tUpdated distance from %lf to %lf\n", dist[Klast], d);
#endif
    /* replace furthest point */
    changed = YES;
    dist[Klast] = d;
    which[Klast] = j;
  }
  /* Bubble sort entries */
  if(changed) {
#ifdef HUH
    Rprintf("Bubble sort.\nCurrent state:\n\tk\twhich\tdist\n");
    for(k = 0; k <= Klast; k++) 
      Rprintf("\t%d\t%d\t%lf\n", k, which[k], dist[k]);
#endif
    do {
      unsorted = NO;
      for(k = 0; k < Klast; k++) {
	if(dist[k] > dist[k+1]) {
	  unsorted = YES;
	  dtmp = dist[k];   dist[k] = dist[k+1];   dist[k+1] = dtmp;
	  itmp = which[k]; which[k] = which[k+1]; which[k+1] = itmp;
	}
      }
    } while(unsorted);
  }
#ifdef HUH
    Rprintf("Return state:\n\tk\twhich\tdist\n");
    for(k = 0; k <= Klast; k++) 
      Rprintf("\t%d\t%d\t%lf\n", k, which[k], dist[k]);
#endif
  return( (int) changed);
}
