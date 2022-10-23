/*

  linvdist.h

  Distance function at vertices
  (shortest distance from each vertex to a data point)

  Function body definitions with macros

  Sparse representation of network

  $Revision: 1.5 $  $Date: 2022/10/21 10:43:01 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  Macros used:
  FNAME   name of function
  WHICH   whether 'nnwhich' is required
  HUH     debugging flag

  ! Data points must be ordered by segment index !

*/

void FNAME(
   /* target data points in local coordinates (ordered by sp) */	   
   int *np,
   int *sp,
   double *tp,
   /* network */
   int *nv,           /* number of network vertices */
   int *ns,           /* segments */
   int *from,
   int *to,
   double *seglen,       /* segment lengths */
   double *huge,         /* value taken as infinity */
   double *tol,          /* tolerance for updating distances */
   /* OUTPUT */
#ifdef WHICH
   double *dist,         /* distance from each vertex to nearest data point */
   int *which         /* identifies nearest data point */
#else 
   double *dist          /* distance from each vertex to nearest data point */
#endif	   
) 
{
  int Np, Nv, Ns, i, j, k, segPj, ivleft, ivright;
  double hugevalue, eps, dleft, dright, slen, d, tpj;
  char converged;

  Np = *np;
  Nv = *nv;
  Ns = *ns;
  hugevalue = *huge;
  eps = *tol;

#ifdef HUH
  Rprintf("Initialise dist\n");
#endif
  /* initialise to huge value */
  for(i = 0; i < Nv; i++) {
    dist[i] = hugevalue;
#ifdef WHICH
    which[i] = -1;
#endif
  }

#ifdef HUH
  Rprintf("Run through target points\n");
#endif
  /* assign correct value to endpoints of segments containing target points */
  for(j = 0; j < Np; j++) {
    segPj = sp[j];
    tpj = tp[j];
    slen = seglen[segPj];
    ivleft = from[segPj];
    d = slen * tpj;
    if(d < dist[ivleft]) {
      dist[ivleft] = d;
#ifdef WHICH
      which[ivleft] = j;
#endif
    }
    ivright = to[segPj];
    d = slen * (1.0 - tpj);
    if(d < dist[ivright]) {
      dist[ivright] = d;
#ifdef WHICH
      which[ivright] = j;
#endif
    }
  }

  /* recursively update */
#ifdef HUH
  Rprintf("Recursive update\n");
#endif
  converged = NO;
  while(!converged) {
    converged = YES;
#ifdef HUH
    Rprintf("........... starting new pass ...................... \n");
#endif
    for(k = 0; k < Ns; k++) {
      ivleft = from[k];
      ivright = to[k];
      slen = seglen[k];
      dleft = (double) dist[ivleft];
      dright = (double) dist[ivright];
      d = (double) (dleft + slen);
      if(d < dright - eps) {
#ifdef HUH
	Rprintf("Updating ivright=%d using ivleft=%d, from %lf to %lf+%lf=%lf\n",
		ivright, ivleft, dright, dleft, slen, d);
#endif
	converged = NO;
	dist[ivright] = d;
#ifdef WHICH
	which[ivright] = which[ivleft];
#endif
      } else {
	d = (double) (dright + slen);
	if(d < dleft - eps) {
#ifdef HUH
	Rprintf("Updating ivleft=%d using ivright=%d, from %lf to %lf+%lf=%lf\n",
		ivleft, ivright, dleft, dright, slen, d);
#endif
	  converged = NO;
	  dist[ivleft] = d;
#ifdef WHICH
	  which[ivleft] = which[ivright];
#endif
	}
      }
    }
  }

#ifdef HUH
  Rprintf("Done\nVertex values:\n");
#ifdef WHICH
  Rprintf("\ti\twhich\tdist\n");
  for(i = 0; i < Nv; i++) 
    Rprintf("\t%d\t%d\t%lf\n", i, which[i], dist[i]);
#else
  Rprintf("\ti\tdist\n");
  for(i = 0; i < Nv; i++) 
    Rprintf("\t%d\t%lf\n", i, dist[i]);
#endif
#endif
}
