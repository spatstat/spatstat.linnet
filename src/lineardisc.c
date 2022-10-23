#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

/* 
   lineardisc.c

   Disc of radius r in linear network

   Clineardisc    determine the linear disc (NOT USED)
   
   Ccountends     count the number of endpoints

   $Revision: 1.15 $  $Date: 2022/10/21 10:43:01 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#define DPATH(I,J) dpath[(J) + Nv * (I)]

#include "yesno.h"

#undef DEBUG

#ifdef DEBUG
void 
Clineardisc(
  /* local coordinates of centre of disc */
  double *f,  /* fractional distance along segment ('tp') */
  int *seg,
  /* radius of disc */
  double *r,
  /* network vertices */
  int *nv, double *xv, double *yv,
  /* segments of network */
  int *ns, int *from, int *to,       
  /* matrix of shortest path distances between vertices */
  double *dpath,
  /* segment lengths */  
  double *lengths, 
  /* outputs */
  int *allinside,     /* whether each segment is entirely inside disc */
  int *boundary,      /* whether each segment crosses the perimeter */
  double *dxv,        /* distance from centre of disc to each vertex */
  int *nendpoints)
{
  int Nv, Ns;
  double f0, rad;
  int seg0;

  int i, A, B, fromi, toi, allin, bdry, reachable, nends, maxchunk;
  double length0, dxA, dxB, dxAvi, dxBvi, residue;
  double *resid; 
  int *covered;

  Nv = *nv;
  Ns = *ns;

  f0 = *f;
  seg0 = *seg;
  rad = *r;

  /* endpoints of segment containing centre */
  A = from[seg0];
  B = to[seg0];

  /* distances from x to  A and B */
  length0 = lengths[seg0];
  dxA = f0 * length0;
  dxB = (1-f0) * length0;

  /* visit vertices */
  covered = (int *) R_alloc((size_t) Nv, sizeof(int));
  resid = (double *) R_alloc((size_t) Nv, sizeof(double));

  OUTERCHUNKLOOP(i, Nv, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Nv, maxchunk, 16384) {
      /* distance going through A */
      dxAvi = dxA + DPATH(A,i);
      /* distance going through B */
      dxBvi = dxB + DPATH(B,i);
      /* shortest path distance to this vertex */
      dxv[i] = (dxAvi < dxBvi) ? dxAvi : dxBvi;
      /* distance left to 'spend' from this vertex */
      residue = rad - dxv[i];
      resid[i] = (residue > 0)? residue : 0;
      /* determine whether vertex i is inside the disc of radius r */
      covered[i] = (residue >= 0);
    }
  }
  /* 
     Now visit line segments. 
  */
  nends = 0;

  OUTERCHUNKLOOP(i, Ns, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Ns, maxchunk, 16384) {
      /* 
	 Determine which line segments are completely inside the disc,
	 and which cross the boundary.
      */
      if(i == seg0) {
	/* initial segment: disc starts from centre (x, y) */
	allin = covered[A] && covered[B];
	bdry  = !allin;
	if(bdry) {
	  if(!covered[A]) nends++;
	  if(!covered[B]) nends++;
	}
      } else {
	/* another segment: disc extends in from either endpoint */
	fromi = from[i];
	toi   = to[i];
	reachable = (covered[fromi] || covered[toi]);
	if(reachable) {
	  allin = covered[fromi] && covered[toi] && 
                     (resid[fromi] + resid[toi] >= lengths[i]);
	  bdry = !allin;
	} else allin = bdry = NO;
	if(bdry) {
	  if(covered[fromi]) nends++;
	  if(covered[toi]) nends++;
	}
      }
      allinside[i] = allin;
      boundary[i] = bdry;
    }
  }
  *nendpoints = nends;
}
#endif

/* ------------------------------------------------- */
/*   count endpoints of several discs in a network   */
/* ------------------------------------------------- */

void 
Ccountends(
  /* centres of discs (local coords) */
  int *np,
  double *f,
  int *seg,
  /* radii of discs */
  double *r,
  /* network vertices */
  int *nv, double *xv, double *yv,
  /* network segments */
  int *ns, int *from, int *to,
  double *dpath,   /* shortest path distances between vertices */
  double *lengths,   /* segment lengths */
  double *toler, /* tolerance */
  int *nendpoints /* output counts of endpoints */
) {
  int Np, Nv, Ns;
  double f0, rad;
  int seg0;

  int i, m, A, B, fromi, toi, reachable, nends, maxchunk, covfrom, covto, allin;
  double length0, dxA, dxB, dxAvi, dxBvi, dxvi, residue, resfrom, resto, tol;
  double *resid; 
  int *covered, *terminal;

  Np = *np;
  Nv = *nv;
  Ns = *ns;
  tol = *toler;

#ifdef DEBUG
  Rprintf("\nTolerance = %lf\n", tol);
#endif

  covered = (int *) R_alloc((size_t) Nv, sizeof(int));
  terminal = (int *) R_alloc((size_t) Nv, sizeof(int));
  resid = (double *) R_alloc((size_t) Nv, sizeof(double));

  /* loop over centre points */
  OUTERCHUNKLOOP(m, Np, maxchunk, 256) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(m, Np, maxchunk, 256) {

      f0 = f[m];
      seg0 = seg[m];
      rad = r[m];

#ifdef DEBUG
      Rprintf("\nCentre point %d lies in segment %d\n", m, seg0);
#endif

      /* endpoints of segment containing centre */
      A = from[seg0];
      B = to[seg0];

      /* distances from centre to A and B */
      length0 = lengths[seg0];
      dxA = f0 * length0;
      dxB = (1-f0) * length0;

#ifdef DEBUG
      Rprintf("Distances to endpoints: dxA=%lf, dxB=%lf\n", dxA, dxB);
#endif

      nends = 0;

      /* visit vertices */
      for(i = 0; i < Nv; i++) {
#ifdef DEBUG
	Rprintf("\nConsidering vertex %d\n", i);
#endif
	/* distance going through A */
	dxAvi = dxA + DPATH(A,i);
	/* distance going through B */
	dxBvi = dxB + DPATH(B,i);
	/* shortest path distance to this vertex */
	dxvi = (dxAvi < dxBvi) ? dxAvi : dxBvi;
	/* distance left to 'spend' from this vertex */
	residue = rad - dxvi;
#ifdef DEBUG
	Rprintf("dxAvi = %lf; dxBvi = %lf; residue = %lf\n", 
		dxAvi, dxBvi, residue);
#endif
	if(residue > tol) {
	  resid[i] = residue;
	  covered[i] = YES;
	  terminal[i] = NO;
#ifdef DEBUG
	  Rprintf("Vertex is covered\n");
#endif
	} else if(residue < -tol) {
	  resid[i] = 0;
	  covered[i] = terminal[i] = NO;
#ifdef DEBUG
	  Rprintf("Vertex is not covered\n");
#endif
	} else {
	  /* vertex is within 'tol' of an endpoint 
	   - deem it to be one 
	  */
	  resid[i] = 0;
	  covered[i] = terminal[i] = YES;
	  /* vertex is an endpoint of disc */
	  ++nends;  
#ifdef DEBUG
	  Rprintf("Vertex is a terminal endpoint\n");
#endif
	}
      }

#ifdef DEBUG
      Rprintf("%d terminal endpoints\n", nends);
#endif

      /* 
	 Now visit line segments 
	 to count any endpoints that are interior to the segments.
      */

      for(i = 0; i < Ns; i++) {
	/* 
	   Determine which line segments are completely inside the disc,
	   and which cross the boundary.
	*/
	if(i == seg0) {
	  /* initial segment: disc starts from (x0, y0) */
	  if(!covered[A]) nends++;
	  if(!covered[B]) nends++;
#ifdef DEBUG
	  if(!covered[A]) Rprintf("A not covered\n");
	  if(!covered[B]) Rprintf("B not covered\n");
#endif
	} else {
	  /* another segment: disc extends in from either endpoint */
	  fromi = from[i];
	  toi   = to[i];
	  covfrom = covered[fromi];
	  covto   = covered[toi];
	  resfrom = resid[fromi];
	  resto   = resid[toi];
	  reachable = covfrom || covto;
#ifdef DEBUG
	  residue = resfrom + resto - lengths[i];
	  Rprintf("%d: %s %s: %lf + %lf - %lf = %lf sign %s\n", 
		  i,
		  (terminal[fromi]) ? "T" : ((covfrom) ? "Y" : "N"),
		  (terminal[toi]) ? "T" : ((covto) ? "Y" : "N"),
		  resfrom, resto, lengths[i], residue,
		  (residue < 0) ? "-" : ((residue > 0) ? "+" : "0"));
#endif
	  if(reachable) {
	    residue = resfrom + resto - lengths[i];
	    allin = covfrom && covto && (residue >= 0);
#ifdef DEBUG
	    if(allin) {
	      Rprintf("Covered\n"); 
	    } else if((terminal[fromi] || terminal[toi]) &&
		      (residue >= - tol * lengths[i])) {
		Rprintf("Deemed to be covered\n"); 
	    } else Rprintf("Reachable\n");
#endif
	    allin = allin || 
	      ((terminal[fromi] || terminal[toi]) &&
	       (residue >= - tol));
	    if(!allin) {
	      /* segment is not entirely covered by disc
		 - infer endpoint(s) in interior of segment */
	      if(covfrom && !terminal[fromi]) nends++;
	      if(covto && !terminal[toi]) nends++;
#ifdef DEBUG
	      if(covfrom && !terminal[fromi]) Rprintf("fromi => end\n");
	      if(covto && !terminal[toi]) Rprintf("toi => end\n");
#endif
	    }
	  }
	}
      }
      nendpoints[m] = nends;
    }
  }
}
