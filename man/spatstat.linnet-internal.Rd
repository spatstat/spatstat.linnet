\name{spatstat.linnet-internal} 
\title{Internal spatstat.linnet functions}

%% Linear networks code is corralled here
\alias{ApplyConnected}
\alias{DoCountCrossEnds}
\alias{DoCountEnds}
\alias{FDMKERNEL}
\alias{as.linfun.linfun}
\alias{as.owin.lintess}
\alias{default.linnet.tolerance}
\alias{makeLinnetTolerance}
\alias{print.lintess}
\alias{print.summary.linim}
\alias{print.summary.linnet}
\alias{print.summary.lintess}
\alias{resolve.heat.steps}
\alias{summary.lintess}
\alias{nobjects.lintess}
\alias{Window.lintess}
\alias{Window<-.linnet}
\alias{Window<-.lpp}
\alias{densitypointsLPP}
\alias{evalCovar.lppm}
\alias{flatdensityfunlpp}
\alias{flatdensityatpointslpp}
\alias{getlambda.lpp}
\alias{local2lpp}
\alias{looHeatLPP}
\alias{looVoronoiLPP}
\alias{validate.lpp.coords}
\alias{as.ppm.lppm}
\alias{pointsAlongNetwork}
\alias{linearKengine}
\alias{linearKmulti}
\alias{linearKmulti.inhom}
\alias{linearKmultiEngine}
\alias{linearpcfengine}
\alias{linearpcfmulti}
\alias{linearpcfmulti.inhom}
\alias{linearPCFmultiEngine}
\alias{resampleNetworkDataFrame}
\alias{sortalongsegment}
\alias{vnnFind}
\alias{ldtEngine}
\alias{qkdeEngine}
\alias{Math.linimlist}
\alias{Ops.linimlist}
\alias{Summary.linimlist}
\alias{Complex.linimlist}
\alias{LinimOp}
\alias{LinimListOp}


%%%%%%%
\description{
  Internal spatstat.linnet functions.
}
\usage{

%% Linear networks code is corralled here
ApplyConnected(X, Engine, r, \dots, rule, auxdata)
DoCountEnds(X, D, toler)
DoCountCrossEnds(X, I, J, DIJ, toler)
FDMKERNEL(lppobj, dtt, dtx, M, nsave, weights,
          stepnames, setuponly, verbose)
\method{as.linfun}{linfun}(X, \dots)
\method{as.owin}{lintess}(W, \dots)
default.linnet.tolerance(L)
makeLinnetTolerance(toler)
\method{print}{lintess}(x, \dots)
\method{print}{summary.linim}(x, \dots)
\method{print}{summary.linnet}(x, \dots)
\method{print}{summary.lintess}(x, \dots)
\method{summary}{lintess}(object, \dots)
\method{nobjects}{lintess}(x)
\method{Window}{lintess}(X, \dots)
\method{Window}{linnet}(X, \dots, check=TRUE) <- value
\method{Window}{lpp}(X, \dots, check=TRUE) <- value
densitypointsLPP(x, sigma, \dots,
                 weights, nsigma, leaveoneout, fast,
                 fastmethod, floored,
                 dx, dt, iterMax, verbose, debug)
\method{evalCovar}{lppm}(model, covariate, \dots, lambdatype, 
          eps, nd, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset)
flatdensityfunlpp(X, \dots, disconnect, weights, what)
flatdensityatpointslpp(X, \dots, leaveoneout, disconnect, weights, what)
getlambda.lpp(lambda, X, subset, \dots,
              update, leaveoneout, loo.given, lambdaname)
local2lpp(L, seg, tp, X, df.only)
looHeatLPP(U0, Amatrix, npts, niter, nsave,
           lixelweight, lixelmap, verbose) 
looVoronoiLPP(X)
validate.lpp.coords(X, fatal, context)
\method{as.ppm}{lppm}(object)
pointsAlongNetwork(L, delta)
linearKengine(X, \dots, r, reweight, denom, correction,
              ratio, showworking)
linearKmulti(X, I, J, r, \dots, correction)
linearKmulti.inhom(X, I, J, lambdaI, lambdaJ, r, \dots, correction,
             normalise)
linearpcfengine(X, \dots, r, reweight, denom, correction, ratio)
linearpcfmulti(X, I, J, r, \dots, correction)
linearpcfmulti.inhom(X, I, J, lambdaI, lambdaJ, r, \dots,
                     correction, normalise)
linearKmultiEngine(X, I, J, \dots,
                   r, reweight, denom, correction, showworking)
linearPCFmultiEngine(X, I, J, \dots,
                   r, reweight, denom, correction, showworking)
resampleNetworkDataFrame(df, template)
sortalongsegment(df)
vnnFind(seg, tp, ns, nv, from, to, seglen, huge, tol, kmax)
ldtEngine(nv, ns, from, to, seglen, huge,
          coUXord, vnndist, vnnwhich, vnnlab)
resolve.heat.steps(sigma, \dots, dx, dt,
                   niter, iterMax, nsave,
                   seglengths, maxdegree, AMbound, L,
                   finespacing, fineNsplit, fineNlixels,
                   W, eps, dimyx, xy, 
                   allow.adjust, warn.adjust,
                   verbose, stepnames)
qkdeEngine(X, sigma, \dots, at, what,
           leaveoneout, diggle, raw, edge2D, edge,
           weights, varcov, positive, shortcut,
           precomputed, savecomputed)
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm = FALSE)}
%NAMESPACE S3method("Math", "linimlist")
%NAMESPACE S3method("Ops", "linimlist")
%NAMESPACE S3method("Complex", "linimlist")
%NAMESPACE S3method("Summary", "linimlist")
LinimOp(e1, e2, op)
LinimListOp(e1, e2, op)

}
\details{
  These internal \pkg{spatstat.linnet} functions should not be called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.linnet} to the next.
}
\value{
The return values of these functions are not documented,
and may change without warning.
}
\keyword{internal}
