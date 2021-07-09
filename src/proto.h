#include <R.h>
#include <Rinternals.h>

/*
  Prototype declarations for all native routines in spatstat.linnet package

  Automatically generated - do not edit! 

*/

/*
  
                  Functions invoked by .C

*/

void linScrossdist(int *, int *, double *, int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *); 
void lincrossdist(int *, double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, int *, double *);
void heatApprox(int *, double *, double *, double *, double *, int *, int *, int *, double *);
void Ccountends(int *, double *, int *, double *, int *, double *, double *, int *, int *, int *, double *, double *, double *, int *);
void Clinequad(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *, double *, int *, double *, double *, int *, double *, double *, int *); 
void ClineRquad(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *, double *, int *, double *, double *, int *, double *, double *, int *); 
void ClineMquad(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *, double *, int *, int *, double *, double *, int *, double *, double *, int *, int *, double *, double *, int *); 
void ClineRMquad(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *, double *, int *, int *, double *, double *, int *, double *, double *, int *, int *, double *, double *, int *);
void linearradius(int *, int *, int *, double *, int *, double *, double *, double *);
void lintileindex(int *, int *, double *, int *, int *, double *, double *, int *, int *);
void Clixellate(int *, int *, int *, int *, int *, int *, double *, double *, int *, double *, int *, int *, int *, double *, int *, double *);
void linnndist(int *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, double *, double *); 
void linknnd(int *, int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *, int *); 
void linnnwhich(int *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, double *, double *, int *); 
void linknnd(int *, int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *, int *); 
void linknncross(int *, int *, int *, double *, int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *, int *); 
void linSnndwhich(int *, int *, double *, int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *, int *); 
void linndcross(int *, double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, int *, double *, double *, int *); 
void linndxcross(int *, double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, int *, int *, int *, double *, double *, int *);
void Clinvwhichdist(int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *, int *); 
void linvknndist(int *, int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *, int *);
void linSpairUdist(int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *); 
void linSpairdist(int *, int *, double *, int *, int *, int *, int *, double *, double *, double *, double *); 
void linpairdist(int *, double *, double *, int *, double *, double *, double *, int *, int *, double *, int *, double *);
/*

             Functions invoked by .Call

*/
SEXP depthrel(SEXP, SEXP, SEXP, SEXP, SEXP);
