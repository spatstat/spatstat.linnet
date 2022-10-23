/* 
   linSpairdist.h

   Function body definitions with macros
   (included in linSpairdist.c)

   Pairwise distances

   Sparse representation of network

   $Revision: 1.5 $  $Date: 2022/10/21 10:43:01 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

   Macros used:
   FNAME    name of function
   SYMMETRIC  whether to exploit symmetry of distance matrix
   VERBOSE  debugging

   ! Data points must be ordered by segment index !

   ! Result must be initialised to zero !

*/

void 
FNAME(
  /* data points in local coordinates (ordered by sp) */
  int *np,
  int *sp,
  double *tp,
  /* network */
  int *nv, /* number of network vertices */
  int *ns,  /* segments */
  int *from,
  int *to, 
  double *seglen,  /* segment lengths */
  double *huge, /* value taken as infinity */
  double *tol, /* tolerance for updating distances */
  /* OUTPUT */
  double *dist  /* matrix of distances from i to j , INITIALISED TO ZERO */
)
{
  int Np, Nv, i, j, ivleft, ivright, spi, spj;
  double dleft, dright, dij, slen, tpi, tpj;
  double *dminvert;  /* min dist from each vertex */
  int one;

  Np = *np;
  Nv = *nv;

  if(Np <= 1) return;
  
  one = 1;

  dminvert = (double *) R_alloc(Nv, sizeof(double));

#ifdef VERBOSE
  Rprintf("Start loop through target points j\n");
#endif

#ifdef SYMMETRIC  
  for(j = 1; j < Np; j++)
#else
  for(j = 0; j < Np; j++)
#endif    
  {
    R_CheckUserInterrupt();
      
    spj = sp[j];
    tpj = tp[j];
      
#ifdef VERBOSE
    Rprintf("Target point %d\n\t lies on segment %d\n", j, spj);
    Rprintf("Compute distance to target from each vertex..\n");
#endif
      
    /* First compute min distance to target point j from each vertex */
    Clinvdist(&one, sp+j, tp+j,
	      nv, ns, from, to, seglen, huge, tol,
	      dminvert);

#ifdef VERBOSE
    Rprintf("Run through source points..\n");
#endif

#ifdef SYMMETRIC
    for(i = 0; i < j; i++) 
#else      
    for(i = 0; i < Np; i++)
#endif
    {
      tpi = tp[i];
      spi = sp[i];   /* segment containing this point */
      slen = seglen[spi];
      
      if(spi == spj) {
	/* target point lies in the same segment */
	dij = slen * fabs(tpj - tpi);
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
	Rprintf("\tDistance via left endpoint %d is %lf\n", ivleft, dleft);
	Rprintf("\tDistance via right endpoint %d is %lf\n", ivright, dright);
#endif
      }
#ifdef VERBOSE
      Rprintf("\tAssigning distance d[%d, %d] = %lf\n", i, j, dij);
#endif
      dist[i + j * Np] = dij;
#ifdef SYMMETRIC
      dist[j + i * Np] = dij;
#endif      
    }
  }
}



