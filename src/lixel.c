#include <R.h>
#include <math.h>

/* 
   lixel.c

   divide a linear network into shorter segments

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

void Clixellate(
     /* 
	A linear network with *ns segments and *nv vertices
	is specified by the vectors from, to, xv, yv.

	The i-th segment will be subdivided into nsplit[i] subsegments.

	New data will be added at the end of the vectors 'xv' and 'yv' 
	representing additional vertices in the new network.

	The point pattern data (*np points with local coordinates sp, tp
	in the coarse network) will be mapped to the new 'fine' network. 
	Points are sorted by 'spcoarse' value.
	
	'xv', 'yv', 'svcoarse', 'tvcoarse'
        must each have space for (nv + sum(nsplit-1)) entries.

	'fromfine', 'tofine' must have length = sum(nsplit).

      */
  int *ns,          /* number of segments (input & output) */
  int *fromcoarse,
  int *tocoarse,    /* endpoints of each segment (input) */
  int *fromfine,
  int *tofine,      /* endpoints of each segment (output) */
  int *nv,          /* number of vertices (input & output) */
  double *xv,
  double *yv,       /* cartesian coords of vertices (input & output) */
  int *svcoarse,    /* segment id of new vertex in COARSE network */
  double *tvcoarse, /* location coordinate of new vertex on COARSE network */
  int *nsplit,      /* numbers of pieces into which each segment should split */
  int *np,          /* number of data points */
  int *spcoarse,    /* segment id coordinate in coarse network*/
  double *tpcoarse, /* location coordinate in coarse network */
  int *spfine,      /* segment id coordinate */
  double *tpfine    /* location coordinate */
) {
  int Np, oldNs, oldNv, i, j, k, ll;
  int oldfromi, oldtoi, newlines, newNv, newNs, SegmentForData;
  double xstart, xend, ystart, yend, xincr, yincr, tn, tpk;

  Np = *np;
  newNv = oldNv = *nv;
  oldNs = *ns;
  newNs = 0;

  /* 
     initialise pointer at start of point pattern
     Determine which segment contains first point
  */
  k = 0;
  SegmentForData = (Np > 0) ? spcoarse[0] : -1;
  
  /* loop over line segments in original network */
  for(i = 0; i < oldNs; i++) {

    newlines = nsplit[i]; 
    oldfromi = fromcoarse[i];
    oldtoi   = tocoarse[i];
    
    /* local coordinates of endpoints of segment, in ***coarse*** network */
    svcoarse[oldfromi] = svcoarse[oldtoi] = i;
    tvcoarse[oldfromi] = 0.0;
    tvcoarse[oldtoi] = 1.0;
    
    if(newlines == 1) {
      /* copy existing segment to new segment list */
      fromfine[newNs] = oldfromi;
      tofine[newNs]   = oldtoi;
      /* advance pointer */
      ++newNs;
    } else if(newlines > 1) {
      /* split segment into 'newlines' pieces */
      xstart = xv[oldfromi];
      ystart = yv[oldfromi];

      xend = xv[oldtoi];
      yend = yv[oldtoi];
    
      xincr = (xend-xstart)/newlines;
      yincr = (yend-ystart)/newlines;

      for(j = 1; j < newlines; j++) {
	/* create new vertex, number 'newNv' */
	xv[newNv] = xstart + j * xincr;
	yv[newNv] = ystart + j * yincr;
	/* local coordinates of new vertex relative to ***coarse*** network */
	svcoarse[newNv] = i;
	tvcoarse[newNv] = ((double) j)/((double) newlines);
	/* create new segment, number 'newNs', ending at new vertex */
	fromfine[newNs] = (j == 1) ? oldfromi : (newNv-1);
	tofine[newNs]   = newNv;
	/* advance */
	++newNv;
	++newNs;
      }
      /* create segment from last added vertex to end of old segment */
      fromfine[newNs] = newNv-1;
      tofine[newNs] = oldtoi;
      ++newNs;
    }

    /* handle any data points lying on current segment i */
    while(SegmentForData == i) {
      if(newlines == 1) {
	spfine[k] = spcoarse[k];
	tpfine[k] = tpcoarse[k];
      } else {
	/* map location on coarse segment to fine segment */
	tn = tpcoarse[k] * newlines;
	ll = (int) floor(tn); 
	ll = (ll < 0) ? 0 : (ll >= newlines) ? (newlines - 1): ll;
	tpk = tn - ll;
	if(tpk < 0.0) tpk = 0.0; else if(tpk > 1.0) tpk = 1.0;
	tpfine[k] = tpk;
	spfine[k] = newNs - newlines + ll;
      }
      /* advance to next data point */
      ++k;
      SegmentForData = (k < Np) ? spcoarse[k] : -1;
    }
  }
  *nv = newNv;
  *ns = newNs;
}
