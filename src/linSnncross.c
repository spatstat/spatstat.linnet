#include <R.h>
#include "yesno.h"

/* 
   linSnncross.c

   Shortest-path distances between nearest neighbours in linear network
   One pattern to another pattern

   $Revision: 1.5 $  $Date: 2022/10/22 10:09:51 $

   'Sparse version' 

   Works with sparse representation
   Does not allow 'exclusion'
   Requires point data to be ordered by segment index.

   linSnndcross      
   linSnndwhich

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

/* functions from linvdist.c */
void Clinvdist(int *np, int *sp, double *tp,
	       int *nv, int *ns, int *from, int *to,
	       double *seglen, double *huge, double *tol,
	       double *dist); 
void Clinvwhichdist(int *np, int *sp, double *tp,
		    int *nv, int *ns, int *from, int *to,
		    double *seglen, double *huge, double *tol,
		    double *dist, int  *which); 

#undef HUH

/* definition of linSnndcross */
#define FNAME linSnndcross
#undef WHICH
#include "linSnncross.h"

/* definition of linSnndwhich */
#undef  FNAME
#define FNAME linSnndwhich
#define WHICH
#include "linSnncross.h"



