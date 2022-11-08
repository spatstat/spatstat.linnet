/* 

  linequad.h

  Template code, #included several times in linequad.c

  Macros used:
  FUNNAME    function name (unmarked version)
  FMKNAME    function name (marked version)
  ALEA       #defined if grid location should be randomised 
  HUH        #defined if debugging is on
  SWAP       swap macro

  $Revision: 1.5 $ $Date: 2022/11/08 00:44:41 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

void FUNNAME(     /* 
	A linear network with *ns segments and *nv vertices
	is specified by the vectors from, to, xv, yv.

	Data points on the network are specified by *ndat, sdat, tdat.
	*** Assumed to be sorted in order of 'sdat' **

	Dummy points will be placed every 'eps' units along each segment.

	Output vectors:
	     wdat   quadrature weights for the data points
	     wdum   quadrature weights for the dummy points

	     xdum, |
	     ydum, | coordinates of dummy points
	     sdum, |
	     tdum  |
	    
	 Space must be allocated for sum(ceiling(lengths/eps)) dummy points. 
	 
        
      */
  /* ............ INPUT ..............................*/
  /* Network Segments */
  int *ns,    /* number of segments */
  int *from,
  int *to,    /* endpoints of each segment */
  /* Network Vertices */
  int *nv,    /* number of vertices */
  double *xv,
  double *yv, /* cartesian coords of vertices */
  /* Spacing parameter */
  double *eps, /* desired spacing of dummy points */
  /* Data point locations */
  int *ndat, /* number of data points (INPUT) */
  int *sdat, /* segment id (local coordinate) (INPUT) */
  double *tdat, /* location (local coordinate tp) (INPUT) */
  /* ................ OUTPUT............................... */
  double *wdat, /* quadrature weights for data points (OUTPUT) */
  /* Dummy points */
  int *ndum, /* number of dummy points (OUTPUT) */
  double *xdum,
  double *ydum, /* spatial coordinates of dummy points (OUTPUT) */
  int *sdum,    /* segment id (local coordinate tp) (OUTPUT)  */
  double *tdum, /* location (local coordinate) (OUTPUT) */
  double *wdum, /* quadrature weights (OUTPUT) */
  /* SCRATCH SPACE */
  int *maxscratch
) {
  int Nseg, Ndat, Ndum, Lmax, i, j, k, ll, m, fromi, toi;
#ifdef HUH
  int Nvert;
#endif
  int SegmentForData, nwhole, nentries, npieces, npieces1;
  double x0, y0, x1, y1, dx, dy;
  double seglength, ratio, epsilon, rump, epsfrac, rumpfrac, gridstart;
  double tfirst, tlast, tcurrent, plen, w;

  int *serial, *count, *pieceid;
  char *isdata;
  double *tvalue, *pieceweight;

  Nseg  = *ns;
  Ndat  = *ndat;

  Ndum = 0;
  Lmax = *maxscratch;

  epsilon = *eps;

#ifdef HUH
  Nvert = *nv;
  Rprintf("Nseg=%d, Nvert=%d, Ndat=d, Lmax = %d\n\n", Nseg, Nvert, Ndat, Lmax);
#endif

  /* allocate scratch space, one for each data/dummy point in current segment */
  serial = (int *) R_alloc(Lmax, sizeof(int));
  isdata = (char *) R_alloc(Lmax, sizeof(char));
  tvalue = (double *) R_alloc(Lmax, sizeof(double));
  pieceid = (int *) R_alloc(Lmax, sizeof(int));

  /* allocate scratch space, one for each piece of current segment */
  count = (int *) R_alloc(Lmax, sizeof(int));
  pieceweight = (double *) R_alloc(Lmax, sizeof(double));

  /* 
     initialise pointer at start of point pattern
     Determine which segment contains first point
  */
  k = 0;
  SegmentForData = (Ndat > 0) ? sdat[0] : -1;

#ifdef ALEA
  GetRNGstate();
#endif

  /* loop over line segments */
  for(i = 0; i < Nseg; i++) {

#ifdef HUH
    Rprintf("Segment %d\n", i);
#endif

    /* endpoints of segment */
    fromi = from[i];
    toi   = to[i];

    x0 = xv[fromi];
    y0 = yv[fromi];
    x1 = xv[toi];
    y1 = yv[toi];

    dx = x1 - x0;
    dy = y1 - y0;
    seglength = sqrt(dx * dx + dy * dy);

    /* divide segment into pieces of length eps 
       with shorter bits at each end */
    ratio = seglength/epsilon;
    nwhole = (int) floor(ratio);
    if(nwhole > 2 && ratio - nwhole < 0.5) --nwhole;
    rump = (seglength - nwhole * epsilon)/2.0;
    epsfrac = epsilon/seglength;
    rumpfrac = rump/seglength;
    /* 
       There are nwhole+2 pieces, with endpoints
       0, rumpfrac, rumpfrac+epsfrac, rumpfrac+2*epsfrac, ..., 1-rumpfrac, 1
    */

    /* Now place dummy points in these pieces */
#ifdef ALEA
    tfirst = rumpfrac * unif_rand();
#else
    tfirst = rumpfrac/2.0;
#endif
#ifdef HUH
    Rprintf("\tnwhole=%d, epsfrac=%lf, rumpfrac=%lf, tfirst=%lf\n", 
	    nwhole, epsfrac, rumpfrac, tfirst);
    Rprintf("\tsegment length %lf divided into %d pieces\n",
	    seglength, nwhole+2);
#endif

    /* create a new dummy point in each piece */
#ifdef HUH
    Rprintf("\tMaking left dummy point %d\n", Ndum);
#endif
    tvalue[0] = tfirst;
    serial[0] = Ndum;  
    isdata[0] = NO;
    count[0] = 1;
    pieceid[0] = 0;
    xdum[Ndum] = x0 + dx * tfirst;
    ydum[Ndum] = y0 + dy * tfirst;
    sdum[Ndum] = i;
    tdum[Ndum] = tfirst;
    ++Ndum;
    if(nwhole > 0) {
#ifdef HUH
      Rprintf("\tMaking %d middle dummy points\n", nwhole);
#endif
#ifdef ALEA
      gridstart = rumpfrac - unif_rand() * epsfrac;
#else
      gridstart = rumpfrac - epsfrac/2.0;
#endif
      for(j = 1; j <= nwhole; j++) {
	serial[j] = Ndum;
	tvalue[j] = tcurrent = gridstart + ((double) j) * epsfrac;
	isdata[j] = NO;
	count[j] = 1;
	pieceid[j] = j;
	xdum[Ndum] = x0 + dx * tcurrent;
	ydum[Ndum] = y0 + dy * tcurrent;
	sdum[Ndum] = i;
	tdum[Ndum] = tcurrent;
	++Ndum;
      }
    } 
    j = nwhole + 1;
#ifdef HUH
    Rprintf("\tMaking right dummy point %d\n", Ndum);
#endif
    serial[j] = Ndum;
    isdata[j] = NO;
    tvalue[j] = tlast = 1.0 - tfirst;
    count[j] = 1;
    pieceid[j] = j;
    xdum[Ndum] = x0 + dx * tlast;
    ydum[Ndum] = y0 + dy * tlast;
    sdum[Ndum] = i;
    tdum[Ndum] = tlast;
    ++Ndum;

    nentries = npieces = nwhole + 2;
    npieces1 = npieces-1;

    /* add any data points lying on current segment i */
    while(SegmentForData == i) {
#ifdef HUH
      Rprintf("\tData point %d lies on segment %d\n", k, i);
#endif
      serial[nentries] = k;
      tvalue[nentries] = tcurrent = tdat[k];
      isdata[nentries] = YES;
      /* determine which piece contains the data point */
      ll = (int) ceil((tcurrent - rumpfrac)/epsfrac);
      if(ll < 0) ll = 0; else if(ll >= npieces) ll = npieces1;
#ifdef HUH
      Rprintf("\tData point %d mapped to piece %d\n", k, ll);
#endif
      count[ll]++;
      pieceid[nentries] = ll;
      ++nentries;
      ++k;
      SegmentForData = (k < Ndat) ? sdat[k] : -1;
    }

    /* compute counting weights for each piece of segment */
#ifdef HUH
    Rprintf("\tcounting weights..\n");
#endif
    for(ll = 0; ll < npieces; ll++) {
      plen = (ll == 0 || ll == npieces1)? rump : epsilon;
      pieceweight[ll] = plen/count[ll];
    }
    
    /* apply weights to data/dummy points */
#ifdef HUH
    Rprintf("\tdistributing weights..\n");
#endif
    for(j = 0; j < nentries; j++) {
      m = serial[j];
      ll = pieceid[j];
      if(ll >= 0 && ll < npieces) {
	w = pieceweight[ll];
	if(isdata[j]) {
#ifdef HUH
	  Rprintf("\t\tEntry %d: data point %d, piece %d\n", j, m, ll);
#endif
	  wdat[m] = w;
	} else {
#ifdef HUH
	  Rprintf("\t\tEntry %d: dummy point %d, piece %d\n", j, m, ll);
#endif
	  wdum[m] = w;
	}
      }
    }
  }

  *ndum = Ndum;

#ifdef ALEA
  PutRNGstate();
#endif

}

void FMKNAME(
  /* NETWORK */	     
  int *ns,       /* number of segments */	     
  int *from,
  int *to,       /* endpoints of each segment */
  int *nv,       /* number of vertices */
  double *xv,
  double *yv,    /* cartesian coords of vertices */
  /* PARAMETERS */
  double *eps,   /* desired spacing of dummy points */
  /* DATA POINTS */
  int *ntypes,   /* number of types (INPUT) */
  int *ndat,     /* number of data points (INPUT)  */
  double *xdat,
  double *ydat,  /* spatial coordinates of data points (INPUT) */
  int *mdat,     /* mark values of data points (INPUT) */
  int *sdat,     /* segment id (local coordinate) of data points (INPUT) */
  double *tdat,  /* location (local coordinate) of data points (INPUT) */
  double *wdat,  /* quadrature weights of data points (>>OUTPUT<<)  */
  /* DUMMY POINTS */
  int *ndum,     /* number of dummy points (OUTPUT) */
  double *xdum,
  double *ydum,  /* spatial coordinates of dummy points (OUTPUT) */
  int *mdum,     /* mark values of dummy points (OUTPUT) */
  int *sdum,     /* segment id (local coordinate) of dummy points (OUTPUT) */
  double *tdum,  /* location (local coordinate) of dummy points (OUTPUT) */
  double *wdum,  /* quadrature weights of dummy points (OUTPUT) */
  /* SCRATCH SPACE */
  int *maxscratch
){
     /* 
	A linear network with *ns segments and *nv vertices
	is specified by the vectors from, to, xv, yv.

	Data points on the network are specified by 
	*ndat, xdat, ydat, mdat, sdat, tdat.
	*** Assumed to be sorted in order of 'sdat' **

	Dummy points will be placed every 'eps' units along each segment
	and replicated for each possible mark.
	Each data point location is also replicated by dummy points
	with each possible mark except the mark of the data point.

	Output vectors:
	     wdat   quadrature weights for the data points
	     wdum   quadrature weights for the dummy points

	     xdum, |
	     ydum, | coordinates of dummy points
	     sdum, |
	     tdum  |
	    
             mdum    marks for dummy points

	 Space must be allocated for 
	 ntypes * sum(ceiling(lengths/eps)) dummy points. 
	 
        
      */
  int Nseg, Ndat, Ndum, Ntypes, Lmax, i, k, ll, m, fromi, toi;
#ifdef HUH
  int Nvert;
#endif
  int SegmentForData, nwhole, nentries, npieces, npieces1, nMpieces;
  int jpiece, jentry, jpdata, type, mcurrent;
  double x0, y0, x1, y1, dx, dy, xcurrent, ycurrent;
  double seglength, ratio, epsilon, rump, epsfrac, rumpfrac, gridstart;
  double tfirst, tlast, tcurrent, plen, w;

  int *serial, *count, *mkpieceid;
  char *isdata;
  double *tvalue, *countingweight;

  Nseg  = *ns;
  Ndat  = *ndat;
  Ntypes = *ntypes;

  Ndum = 0;
  Lmax = *maxscratch;

  epsilon = *eps;

#ifdef HUH
  Nvert = *nv;
  Rprintf("Nseg=%d, Nvert=%d, Ndat=d, Lmax = %d\n\n", Nseg, Nvert, Ndat, Lmax);
#endif

  /* allocate scratch space, one for each data/dummy point in current segment */
  serial = (int *) R_alloc(Lmax, sizeof(int));
  isdata = (char *) R_alloc(Lmax, sizeof(char));
  tvalue = (double *) R_alloc(Lmax, sizeof(double));
  mkpieceid = (int *) R_alloc(Lmax, sizeof(int));

  /* allocate scratch space, one for each piece of current segment */
  count = (int *) R_alloc(Lmax, sizeof(int));
  countingweight = (double *) R_alloc(Lmax, sizeof(double));

  /* 
     initialise pointer at start of point pattern
     Determine which segment contains first point
  */
  k = 0;
  SegmentForData = (Ndat > 0) ? sdat[0] : -1;

#ifdef ALEA
  GetRNGstate();
#endif

  /* loop over line segments */
  for(i = 0; i < Nseg; i++) {

#ifdef HUH
    Rprintf("Segment %d\n", i);
#endif

    /* endpoints of segment */
    fromi = from[i];
    toi   = to[i];

    x0 = xv[fromi];
    y0 = yv[fromi];
    x1 = xv[toi];
    y1 = yv[toi];

    dx = x1 - x0;
    dy = y1 - y0;
    seglength = sqrt(dx * dx + dy * dy);

    /* divide segment into pieces of length eps 
       with shorter bits at each end */
    ratio = seglength/epsilon;
    nwhole = (int) floor(ratio);
    if(nwhole > 2 && ratio - nwhole < 0.5) --nwhole;
    npieces = nwhole + 2;
    rump = (seglength - nwhole * epsilon)/2.0;
    epsfrac = epsilon/seglength;
    rumpfrac = rump/seglength;
    /* 
       There are nwhole+2 pieces, with endpoints
       0, rumpfrac, rumpfrac+epsfrac, rumpfrac+2*epsfrac, ..., 1-rumpfrac, 1
    */

    /* Now place dummy points in these pieces */
#ifdef ALEA
    tfirst = rumpfrac * unif_rand();
    gridstart = rumpfrac - epsfrac * unif_rand();
#else
    tfirst = rumpfrac/2.0;
    gridstart = rumpfrac - epsfrac/2.0;
#endif
    tlast = 1.0 - tfirst;
#ifdef HUH
    Rprintf("\tnwhole=%d, epsfrac=%lf, rumpfrac=%lf, tfirst=%lf\n", 
	    nwhole, epsfrac, rumpfrac, tfirst);
    Rprintf("\tsegment length %lf divided into %d pieces\n",
	    seglength, npieces);
#endif

    /* 
       'Marked pieces' of segment are numbered in order
       (piece 0, mark 0), (piece 0, mark 1), ..., (piece 0, mark Ntypes-1),
       (piece 1, mark 0), .....
       
       mpieceid = type + pieceid * Ntypes

    */


#ifdef HUH
      Rprintf("\tMaking %d x %d = %d dummy points\n", 
	      npieces, Ntypes, npieces * Ntypes);
#endif

    /* create a new dummy point in each piece */
    npieces1 = npieces-1;
    for(jpiece = 0; jpiece < npieces; jpiece++) {
      tcurrent = (jpiece == 0) ? tfirst :
	(jpiece == npieces1) ? tlast : 
	(gridstart + ((double) jpiece) * epsfrac);
      xcurrent = x0 + dx * tcurrent;
      ycurrent = y0 + dy * tcurrent;
      
      for(type = 0; type < Ntypes; type++) {
	/* position in list of relevant data/dummy points */
	jentry = type + jpiece * Ntypes; 
	/* serial number of marked piece */
	ll = jentry; 

	tvalue[jentry] = tcurrent;
	serial[jentry] = Ndum;  
	isdata[jentry] = NO;
	mkpieceid[jentry] = ll;

	count[ll] = 1;

	xdum[Ndum] = xcurrent;
	ydum[Ndum] = ycurrent;
	mdum[Ndum] = type;
	sdum[Ndum] = i;
	tdum[Ndum] = tcurrent;
	++Ndum;
      }
    }
    nentries = npieces * Ntypes;

    /* handle any data points lying on current segment i */
    while(SegmentForData == i) {
#ifdef HUH
      Rprintf("\tData point %d lies on segment %d\n", k, i);
#endif
      xcurrent = xdat[k];
      ycurrent = ydat[k];
      tcurrent = tdat[k];
      mcurrent = mdat[k];
      /* determine which piece contains the data point */
      jpdata = (int) ceil((tcurrent - rumpfrac)/epsfrac);
      if(jpdata < 0) jpdata = 0; else if(jpdata >= npieces) jpdata = npieces1;
#ifdef HUH
      Rprintf("\tData point %d falls in piece %d\n", k, jpdata);
#endif
      /* 
	 copy data point, 
	 and create dummy points at same location with different marks  
      */
      for(type = 0; type < Ntypes; type++) {
	tvalue[nentries] = tcurrent;
	ll = type + jpdata * Ntypes;
	mkpieceid[nentries] = ll; 
	count[ll]++;
	if(type == mcurrent) {
	  /* data point */
	  isdata[nentries] = YES;
	  serial[nentries] = k;
	} else {
	  /* create dummy point */
	  isdata[nentries] = NO;
	  serial[nentries] = Ndum;
	  xdum[Ndum] = xcurrent;
	  ydum[Ndum] = ycurrent;
	  mdum[Ndum] = type;
	  sdum[Ndum] = i;
	  tdum[Ndum] = tcurrent;
	  ++Ndum;
	}
	++nentries;
      }

      ++k;
      SegmentForData = (k < Ndat) ? sdat[k] : -1;
    }

    /* compute counting weights for each piece of segment */
#ifdef HUH
    Rprintf("\tcounting weights..\n");
#endif
    for(jpiece = 0; jpiece < npieces; jpiece++) {
      plen = (jpiece == 0 || jpiece == npieces1)? rump : epsilon;
      for(type = 0; type < Ntypes; type++) {
	ll = type + jpiece * Ntypes;
	countingweight[ll] = plen/count[ll];
      }
    }
    
    /* apply weights to data/dummy points */
#ifdef HUH
    Rprintf("\tdistributing weights..\n");
#endif
    nMpieces = npieces * Ntypes;
    for(jentry = 0; jentry < nentries; jentry++) {
      m = serial[jentry];
      ll = mkpieceid[jentry];
      if(ll >= 0 && ll < nMpieces) {
	w = countingweight[ll];
	if(isdata[jentry]) {
#ifdef HUH
	  Rprintf("\t\tEntry %d: data point %d, piece %d\n", jentry, m, ll);
#endif
	  wdat[m] = w;
	} else {
#ifdef HUH
	  Rprintf("\t\tEntry %d: dummy point %d, piece %d\n", jentry, m, ll);
#endif
	  wdum[m] = w;
	}
      }
    }
  }

  *ndum = Ndum;

#ifdef ALEA
  PutRNGstate();
#endif

}
