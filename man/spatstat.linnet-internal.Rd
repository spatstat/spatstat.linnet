\name{spatstat.linnet-internal} 
\title{Internal spatstat.linnet functions}
\alias{spatstat.linnet-internal} %DoNotExport


%% Linear networks code is corralled here
\alias{ApplyConnected}
\alias{DoCountCrossEnds}
\alias{DoCountEnds}
\alias{FDMKERNEL}
\alias{as.linfun.linfun}
\alias{as.lintess}
\alias{as.lintess.lintess}
\alias{as.lintess.linnet}
\alias{as.lintess.linearquadratcount}
\alias{affine.lintess}
\alias{as.owin.lintess}
\alias{as.owin.linearquadratcount}
\alias{Window.linearquadratcount}
\alias{diagLppmEngine}
\alias{domain.linearquadratcount}
\alias{default.linnet.tolerance}
\alias{evaluateNetCovariate}
\alias{evaluateNetCovariateAlongNetwork}
\alias{evaluateNetCovariateAtPoints}
\alias{exactlppm}
\alias{flipxy.lintess}
\alias{getglmdata.lppm}
\alias{getglmfit.lppm}
\alias{getglmsubset.lppm}
\alias{hasglmfit.lppm}
\alias{intensity.linearquadratcount}
\alias{is.poisson.exactlppm}
\alias{is.stationary.exactlppm}
\alias{makeLinnetTolerance}
\alias{marginalIntegralOfLinim}
\alias{plot.linearquadratcount}
\alias{plot.qqlppm}
\alias{predict.exactlppm}
\alias{print.diaglppm}
\alias{print.exactlppm}
\alias{print.lintess}
\alias{print.qqlppm}
\alias{print.summary.linim}
\alias{print.summary.linnet}
\alias{print.summary.lintess}
\alias{rescale.lintess}
\alias{resolve.heat.steps}
\alias{response.exactlppm}
\alias{rmaxEuclidean}
\alias{rotate.lintess}
\alias{scalardilate.lintess}
\alias{shift.lintess}
\alias{shift.linearquadratcount}
\alias{summary.lintess}
\alias{nobjects.lintess}
\alias{Window.lintess}
\alias{Window<-.linnet}
\alias{Window<-.lpp}
\alias{densitypointsLPP}
\alias{flatdensityfunlpp}
\alias{flatdensityatpointslpp}
\alias{local2lpp}
\alias{looHeatLPP}
\alias{looVoronoiLPP}
\alias{validate.lpp.coords}
\alias{as.ppm.lppm}
\alias{pointsAlongNetwork}
\alias{lineardiscEngine}
\alias{linearEuclidEngine}
\alias{linearKengine}
\alias{linearKmulti}
\alias{linearKmulti.inhom}
\alias{linearKmultiEngine}
\alias{linearpcfengine}
\alias{linearpcfmulti}
\alias{linearpcfmulti.inhom}
\alias{linearPCFmultiEngine}
\alias{resampleNetworkDataFrame}
\alias{resolve.lambda.lpp}
\alias{sortalongsegment}
\alias{spatialCovariateEvidence.lppm}
\alias{spatialCovariateEvidence.exactlppm}
\alias{vnnFind}
\alias{ldtEngine}
\alias{qkdeEngine}
\alias{updateData.lppm}
\alias{Math.linimlist}
\alias{Ops.linimlist}
\alias{Summary.linimlist}
\alias{Complex.linimlist}
\alias{LinimOp}
\alias{LinimListOp}
\alias{traceTessLinnet}


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
as.lintess(X, \dots)
\method{as.lintess}{lintess}(X, \dots)
\method{as.lintess}{linnet}(X, \dots)
\method{as.lintess}{linearquadratcount}(X, \dots)
\method{affine}{lintess}(X, \dots)
\method{as.owin}{lintess}(W, \dots)
\method{as.owin}{linearquadratcount}(W, \dots, fatal)
\method{Window}{linearquadratcount}(X, \dots)
diagLppmEngine(object, \dots, type, typename, opt,
               sigma, compute.sd, compute.cts,
               envelope, nsim, nrank,
               rv, oldstyle, splineargs, verbose)
\method{domain}{linearquadratcount}(X, \dots)
\method{flipxy}{lintess}(X, \dots)
\method{getglmdata}{lppm}(object, \dots)
\method{getglmfit}{lppm}(object, \dots)
\method{getglmsubset}{lppm}(object, \dots)
\method{hasglmfit}{lppm}(object)
default.linnet.tolerance(L)
evaluateNetCovariate(covariate, locations, \dots)
evaluateNetCovariateAlongNetwork(covariate, locations, \dots, types)
evaluateNetCovariateAtPoints(covariate, locations, \dots, allow.column)
exactlppm(X, baseline, \dots, subset)
\method{intensity}{linearquadratcount}(X, \dots, fun)
\method{is.poisson}{exactlppm}(x)
\method{is.stationary}{exactlppm}(x)
makeLinnetTolerance(toler)
marginalIntegralOfLinim(f, margin)
\method{plot}{linearquadratcount}(x, \dots, main, add, do.plot, textargs)
\method{plot}{qqlppm}(x, \dots, limits, monochrome, limcol)
\method{predict}{exactlppm}(object, locations,\dots)
\method{print}{diaglppm}(x, \dots)
\method{print}{exactlppm}(x, \dots)
\method{print}{lintess}(x, \dots)
\method{print}{qqlppm}(x, \dots)
\method{print}{summary.linim}(x, \dots)
\method{print}{summary.linnet}(x, \dots)
\method{print}{summary.lintess}(x, \dots)
\method{rescale}{lintess}(X, s, unitname)
\method{response}{exactlppm}(object)
\method{rotate}{lintess}(X, \dots)
\method{shift}{lintess}(X, \dots)
\method{shift}{linearquadratcount}(X, \dots)
\method{scalardilate}{lintess}(X, f,\dots)
\method{summary}{lintess}(object, \dots)
\method{nobjects}{lintess}(x)
\method{Window}{lintess}(X, \dots)
\method{Window}{linnet}(X, \dots, check=TRUE) <- value
\method{Window}{lpp}(X, \dots, check=TRUE) <- value
densitypointsLPP(x, sigma, \dots,
                 weights, nsigma, leaveoneout, fast,
                 fastmethod, floored,
                 dx, dt, iterMax, verbose, debug)
flatdensityfunlpp(X, \dots, disconnect, weights, what)
flatdensityatpointslpp(X, \dots, leaveoneout, disconnect, weights, what)
local2lpp(L, seg, tp, X, df.only)
looHeatLPP(U0, Amatrix, npts, niter, nsave,
           lixelweight, lixelmap, verbose) 
looVoronoiLPP(X)
validate.lpp.coords(X, fatal, context)
\method{as.ppm}{lppm}(object)
pointsAlongNetwork(L, delta)
lineardiscEngine(L, x, r, want)
linearEuclidEngine(X, fun, \dots, r, reweight, denom,
                   samplesize, showworking, correction)
linearKengine(X, \dots, r, reweight, denom, samplesize,
              correction, ratio, showworking)
linearKmulti(X, I, J, r, \dots, correction)
linearKmulti.inhom(X, I, J, lambdaI, lambdaJ, r, \dots, correction,
             normalise, sigma)
linearpcfengine(X, \dots, r, reweight, denom, samplesize, correction, ratio)
linearpcfmulti(X, I, J, r, \dots, correction)
linearpcfmulti.inhom(X, I, J, lambdaI, lambdaJ, r, \dots,
                     correction, normalise,
                     sigma, adjust.sigma, bw, adjust.bw)
linearKmultiEngine(X, I, J, \dots,
                   r, reweight, denom, samplesize, correction, showworking)
linearPCFmultiEngine(X, I, J, \dots,
                   r, reweight, denom, samplesize, correction, showworking)
resampleNetworkDataFrame(df, template)
\method{resolve.lambda}{lpp}(X, lambda, subset, \dots,
       update, leaveoneout, everywhere, loo.given, sigma, lambdaname)
sortalongsegment(df)
\method{spatialCovariateEvidence}{lppm}(model, covariate, \dots, lambdatype, 
          eps, dimyx, xy, rule.eps,
          delta, nd, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict, raster.action)
\method{spatialCovariateEvidence}{exactlppm}(model, covariate, \dots,
          lambdatype, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict, raster.action)
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
rmaxEuclidean(L, verbose, show)
qkdeEngine(x, sigma, \dots, at, what,
           leaveoneout, diggle, raw, edge2D, edge,
           weights, varcov, positive, shortcut,
           precomputed, savecomputed)
\method{updateData}{lppm}(model, X, \dots)
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
traceTessLinnet(A, L)

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
