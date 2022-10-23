#include <R.h>
#include "yesno.h"

/* 
   linSpairdist.c

   Distances between points on a linear network

   linSpairdist
   linSpairUdist

   'Sparse version' 

   $Revision: 1.3 $  $Date: 2022/10/22 10:09:51 $

   Works with sparse representation
   Requires point data to be ordered by segment index.

   Repeatedly includes 'linSpairdist.h'

   Macros used:
   FNAME      function name
   VERBOSE    debugging
   SYMMETRIC  whether to exploit the symmetry of the distance matrix

   ! Data points must be ordered by segment index !

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef VERBOSE
#undef SYMMETRIC

void Clinvdist(int *np, int *sp, double *tp,
	       int *nv, int *ns, int *from, int *to,
	       double *seglen, double *huge, double *tol,
	       double *dist);  /* function from linvdist.c */

/* This version computes d[i,j] and d[j,i] separately as a check on validity */

#define FNAME linSpairUdist
#include "linSpairdist.h"
#undef FNAME

/* This is the production version which exploits symmetry d[i,j] = d[j,i] */

#define FNAME linSpairdist
#define SYMMETRIC
#include "linSpairdist.h"
#undef FNAME
#undef SYMMETRIC

