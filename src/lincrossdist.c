#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/* 
   lincrossdist.c

   Shortest-path distances between pairs of points in linear network

   $Revision: 1.5 $  $Date: 2022/10/21 10:43:01 $

   lincrossdist

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#define DPATH(I,J) dpath[(I) + Nv * (J)]
#define ANSWER(I,J) answer[(I) + Np * (J)]
#define EUCLID(X,Y,U,V) sqrt(pow((X)-(U),2)+pow((Y)-(V),2))

void 
lincrossdist(
  /* data points from which distances are measured */
  int *np, double *xp, double *yp,
  /* data points to which distances are measured */
  int *nq, double *xq, double *yq,
  /* network vertices */
  int *nv, double *xv, double *yv,
  /* segments of network */
  int *ns, int *from, int *to,
  double *dpath,  /* path distances between vertices */
  int *psegmap, /* map from data points 'p' to segments */
  int *qsegmap, /* map from data points 'q' to segments */
  double *answer /* matrix of path distances from 'p' to 'q' */
) {
  int Np, Nq, Nv, i, j, maxchunk;
  int Psegi, Qsegj, nbi1, nbi2, nbj1, nbj2; 
  double xpi, ypi, xqj, yqj;
  double d, dPiV1, dPiV2, dV1Qj, dV2Qj, d11, d12, d21, d22; 

  Np = *np;
  Nq = *nq;
  Nv = *nv;

  OUTERCHUNKLOOP(i, Np, maxchunk, 1024) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Np, maxchunk, 1024) {
      xpi = xp[i];
      ypi = yp[i];
      Psegi = psegmap[i];
      nbi1 = from[Psegi];
      nbi2 = to[Psegi];
      dPiV1 = EUCLID(xpi, ypi, xv[nbi1], yv[nbi1]);
      dPiV2 = EUCLID(xpi, ypi, xv[nbi2], yv[nbi2]);
      for(j = 0; j < Nq; j++) {
	xqj = xq[j];
	yqj = yq[j];
	Qsegj = qsegmap[j];
	if(Psegi == Qsegj) {
	  /* points i and j lie on the same segment; use Euclidean distance */
	  d = sqrt(pow(xpi - xqj, 2) + pow(ypi - yqj, 2));
	} else {
	  /* Shortest path from i to j passes through ends of segments;
	     Calculate shortest of 4 possible paths from i to j
	  */
	  nbj1 = from[Qsegj];
	  nbj2 = to[Qsegj];
	  dV1Qj = EUCLID(xv[nbj1], yv[nbj1], xqj, yqj);
	  dV2Qj = EUCLID(xv[nbj2], yv[nbj2], xqj, yqj);
	  d11 = dPiV1 + DPATH(nbi1,nbj1) + dV1Qj;
	  d12 = dPiV1 + DPATH(nbi1,nbj2) + dV2Qj;
	  d21 = dPiV2 + DPATH(nbi2,nbj1) + dV1Qj;
	  d22 = dPiV2 + DPATH(nbi2,nbj2) + dV2Qj;
	  d = d11;
	  if(d12 < d) d = d12;
	  if(d21 < d) d = d21;
	  if(d22 < d) d = d22;
	}
	/* write */
	ANSWER(i,j) = d;
      }
    }
  }
}

