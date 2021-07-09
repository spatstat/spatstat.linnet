#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

/*
  depthrel.c

  Find which segments lie in front of/behind other segments
  
      x is the projected x-coordinate
  
      z is the negative depth (z increases as we move closer to the viewer) 

      Assumes X0 <= X1
  
 */

#define ROUNDED(X) ((float) (X))

SEXP depthrel(SEXP X0, 
	      SEXP Z0, 
	      SEXP X1, 
	      SEXP Z1,
	      SEXP Verb) 
{
  double *x0, *z0, *x1, *z1;
  int i, j, k, m, n, kmax, kmaxnew;
  int *front, *back;
  double xleft, xright, dxi, dxj, z0i, z0j, z1i, z1j;
  char verbose;
  int status;

  SEXP out, fout, bout, sout;
  int *fp, *bp;

  PROTECT(X0 = AS_NUMERIC(X0));
  PROTECT(Z0 = AS_NUMERIC(Z0));
  PROTECT(X1 = AS_NUMERIC(X1));
  PROTECT(Z1 = AS_NUMERIC(Z1));
  PROTECT(Verb = AS_INTEGER(Verb));

  x0 = NUMERIC_POINTER(X0);
  z0 = NUMERIC_POINTER(Z0);
  x1 = NUMERIC_POINTER(X1);
  z1 = NUMERIC_POINTER(Z1);
  
  verbose = (*(INTEGER_POINTER(Verb)) == 1);

  n = LENGTH(X0);

  status = 0;

  /* initial guess for max number of related pairs */
  kmax = 4 * (n+1);

  /* allocate space for list of related pairs */
  k = 0;
  front = (int *) R_alloc(kmax, sizeof(int));
  back  = (int *) R_alloc(kmax, sizeof(int));

  if(n >= 2) {
    for(i = 1; i < n; i++) {
      for(j = 0; j < i; j++) {
        if(x1[i] > x0[j] && x1[j] > x0[i]) {
          /* overlap occurs */
          /* consider left side */
          z0i = z0[i];
          z0j = z0[j];
          if(x0[j] < x0[i]) {
            xleft = x0[i];
	    dxj = x1[j]-x0[j];
	    if(dxj != 0.0) 
	      z0j += (z1[j] - z0j) * ((xleft - x0[j])/dxj);
          } else {
            xleft = x0[j];
	    dxi = x1[i]-x0[i];
	    if(dxi != 0.0) 
	      z0i += (z1[i] - z0i) * ((xleft - x0[i])/dxi);
          }
          /* consider right side */
          z1i = z1[i];
          z1j = z1[j];
          if(x1[j] > x1[i]) {
            xright = x1[i];
	    dxj = x1[j]-xleft;
	    if(dxj != 0.0)
	      z1j = z0j + (z1j - z0j) * ((xright - xleft)/dxj);
          } else {
            xright = x1[j];
	    dxi = x1[i]-xleft;
	    if(dxi != 0.0) 
	      z1i = z0i + (z1i - z0i) * ((xright - xleft)/dxi);
          }
          /* now determine which is in front */
          if(ROUNDED(z0i) >= ROUNDED(z0j) && ROUNDED(z1i) >= ROUNDED(z1j)) {
            /* 'i' is in front */
            front[k] = i + 1;
            back[k]  = j + 1;
          } else if(ROUNDED(z0i) <= ROUNDED(z0j) && ROUNDED(z1i) <= ROUNDED(z1j)){
            /* 'j' is in front */
            front[k] = j + 1;
            back[k]  = i + 1;
          } else {
	    if(verbose) 
	      warning("segments %d and %d cross over", i+1, j+1);
	    status = 1;
	  }
	  k++;
	  if(k >= kmax) {
	    /* storage overflow */
	    kmaxnew = 2 * kmax;
	    front = (int *) S_realloc((char *) front,
				      kmaxnew,  kmax, sizeof(int));
	    back  = (int *) S_realloc((char *) back,
				      kmaxnew,  kmax, sizeof(int));
	    kmax =  kmaxnew;
	  }
        }
      }
    }
  }
  /* copy to output */
  PROTECT(out = NEW_LIST(3));
  PROTECT(fout = NEW_INTEGER(k));
  PROTECT(bout = NEW_INTEGER(k));
  PROTECT(sout = NEW_INTEGER(1));
  if(k > 0) {
    fp = INTEGER_POINTER(fout);
    bp = INTEGER_POINTER(bout);
    for(m = 0; m < k; m++) {
      fp[m] = front[m];
      bp[m] = back[m];
    }
  }
  *(INTEGER_POINTER(sout)) = status;

  SET_VECTOR_ELT(out, 0, fout);
  SET_VECTOR_ELT(out, 1, bout);
  SET_VECTOR_ELT(out, 2, sout);

  UNPROTECT(9);
  return(out);
}
