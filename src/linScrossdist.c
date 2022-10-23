#include <R.h>
#include "yesno.h"

/* 
   linScrossdist.c

   Distances between points on a linear network
   One pattern to another pattern

   linScrossdist

   'Sparse version' 

   $Revision: 1.7 $  $Date: 2022/10/22 10:09:51 $

   Works with sparse representation
   Requires point data to be ordered by segment index.

   Macros used:
   VERBOSE    debugging

   ! Data points must be ordered by segment index !

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef VERBOSE

void Clinvdist(int *np, int *sp, double *tp,
	       int *nv, int *ns, int *from, int *to,
	       double *seglen, double *huge, double *tol,
	       double *dist);  /* function from linvdist.c */

void 
linScrossdist(
  /* data points 'from' in local coordinates (ordered by sp) */
  int *np,
  int *sp,
  double *tp,
  /* data points 'to' in local coordinates (ordered by sq) */
  int *nq,
  int *sq,
  double *tq,
  /* network */
  int *nv, /* number of network vertices */
  int *ns,   /* network segments */
  int *from,
  int *to,  
  double *seglen,  /* segment lengths */
  /* parameters */
  double *huge, /* value taken as infinity */
  double *tol, /* tolerance for updating distances */
  /* OUTPUT */
  double *dist  /* matrix of distances from i to j */
)
{
  int Np, Nq, Nv, i, j, ivleft, ivright, spi, sqj;
  double dleft, dright, dij, slen, tpi, tqj;
  double *dminvert;  /* min dist from each vertex */
  int one;

  Np = *np;
  Nq = *nq;
  Nv = *nv;

  one = 1;

  dminvert = (double *) R_alloc(Nv, sizeof(double));

#ifdef VERBOSE
  Rprintf("Start loop through target points j\n");
#endif

  for(j = 0; j < Nq; j++) {
    R_CheckUserInterrupt();
      
    sqj = sq[j];
    tqj = tq[j];
      
#ifdef VERBOSE
    Rprintf("Target point %d\n\t lies on segment %d\n", j, sqj);
    Rprintf("Compute distance to target from each vertex..\n");
#endif
      
    /* First compute min distance to target point j from each vertex */
    Clinvdist(&one, sq+j, tq+j,
	      nv, ns, from, to, seglen, huge, tol,
	      dminvert);

#ifdef VERBOSE
    Rprintf("Run through source points..\n");
#endif
  
    for(i = 0; i < Np; i++) {
      tpi = tp[i];
      spi = sp[i];   /* segment containing this point */
      slen = seglen[spi];
      
      if(spi == sqj) {
	/* target point lies in the same segment */
	dij = slen * fabs(tqj - tpi);
#ifdef VERBOSE
	Rprintf("\tSource %d and target lie on same segment, distance %lf\n",
		i, dij);
#endif
      } else {
	ivleft = from[spi];
	ivright = to[spi];
#ifdef VERBOSE
	Rprintf("\tSource point %d lies on segment %d = [%d,%d]\n", 
		i, spi, ivleft, ivright);
#endif
	dleft  = slen * tpi + dminvert[ivleft];
	dright = slen * (1.0 - tpi) + dminvert[ivright];
	dij = (dleft < dright) ? dleft : dright;
#ifdef VERBOSE
	Rprintf("\tDistance to left endpoint %d is %lf\n", ivleft, dleft);
	Rprintf("\tDistance to right endpoint %d is %lf\n", ivright, dright);
#endif
      }
#ifdef VERBOSE
      Rprintf("\tAssigning distance d[%d, %d] = %lf\n", i, j, dij);
#endif
      dist[i + j * Np] = dij;
    }
  }
}


