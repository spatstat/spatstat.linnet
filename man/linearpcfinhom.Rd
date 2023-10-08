\name{linearpcfinhom}
\alias{linearpcfinhom}
\title{
  Inhomogeneous Linear Pair Correlation Function
}
\description{
  Computes an estimate of the inhomogeneous linear pair correlation function
  for a point pattern on a linear network.
}
\usage{
linearpcfinhom(X, lambda=NULL, r=NULL, ..., correction="Ang",
               normalise=TRUE, normpower=1,
	       update = TRUE, leaveoneout = TRUE,
               sigma=NULL, adjust.sigma=1,
               bw="nrd0", adjust.bw=1,
	       ratio = FALSE)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Intensity values for the point pattern. Either a numeric vector,
    a \code{function}, a pixel image (object of class \code{"im"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.default}}
    to control the smoothing of the estimates of pair correlation.
  }
  \item{correction}{
    Geometry correction.
    Either \code{"none"} or \code{"Ang"}. See Details.
  }
  \item{normalise}{
    Logical. If \code{TRUE} (the default), the denominator of the estimator is 
    data-dependent (equal to the sum of the reciprocal intensities at the data
    points, raised to \code{normpower}), which reduces the sampling variability.
    If \code{FALSE}, the denominator is the length of the network.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See explanation in \code{\link{linearKinhom}}.
  }
  \item{update}{
    Logical value indicating what to do when \code{lambda} is a fitted model
    (class \code{"lppm"} or \code{"ppm"}).
    If \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link{update.lppm}} or \code{\link{update.ppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
  \item{leaveoneout}{
    Logical value specifying whether to use a
    leave-one-out rule when calculating the intensity.
    See Details.
  }
  \item{sigma}{
    Smoothing bandwidth (passed to \code{\link{density.lpp}})
    for kernel density estimation of the intensity when
    \code{lambda=NULL}.
  }
  \item{adjust.sigma}{
    Numeric value. \code{sigma} will be multiplied by this value.
  }
  \item{bw}{
    Smoothing bandwidth (passed to \code{\link[stats]{density.default}})
    for one-dimensional kernel smoothing of the pair correlation function.
    Either a numeric value, or a character string recognised
    by \code{\link[stats]{density.default}}.
  }
  \item{adjust.bw}{
    Numeric value. \code{bw} will be multiplied by this value.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\details{
  This command computes the inhomogeneous version of the 
  linear pair correlation function from point pattern data on a linear network.

  The argument \code{lambda} should provide estimated values
  of the intensity of the point process at each point of \code{X}.

  If \code{lambda=NULL}, the intensity will be estimated by kernel
  smoothing by calling \code{\link{density.lpp}} with the smoothing
  bandwidth \code{sigma}, and with any other relevant arguments
  that might be present in \code{\dots}.  A leave-one-out kernel estimate
  will be computed if \code{leaveoneout=TRUE}.

  If \code{lambda} is given, 
  it may be a numeric vector (of length equal to
  the number of points in \code{X}), or a \code{function(x,y)} that will be
  evaluated at the points of \code{X} to yield numeric values, 
  or a pixel image (object of class \code{"im"}) or a fitted point 
  process model (object of class \code{"ppm"} or \code{"lppm"}).

  If \code{lambda} is a fitted point process model,
  the default behaviour is to update the model by re-fitting it to
  the data, before computing the fitted intensity.
  This can be disabled by setting \code{update=FALSE}.
  The intensity at data points will be computed
  by \code{\link{fitted.lppm}} or \code{\link{fitted.ppm}}.
  A leave-one-out estimate will be computed if \code{leaveoneout=TRUE}
  and \code{update=TRUE}.
  
  If \code{correction="none"}, the calculations do not include
  any correction for the geometry of the linear network.
  If \code{correction="Ang"}, the pair counts are weighted using
  Ang's correction (Ang, 2010). 

  The bandwidth for smoothing the pairwise distances
  is determined by arguments \code{\dots}
  passed to \code{\link{density.default}}, mainly the arguments
  \code{bw} and \code{adjust}. The default is
  to choose the bandwidth by Silverman's rule of thumb 
  \code{bw="nrd0"} explained in \code{\link{density.default}}.
}
\section{Warning}{
  Older versions of \code{\link{linearpcfinhom}} interpreted
  \code{lambda=NULL} to mean that the homogeneous function
  \code{\link{linearpcf}} should be computed. This was changed to the
  current behaviour in version \code{3.1-0} of \pkg{spatstat.linnet}.
}
\value{
  Function value table (object of class \code{"fv"}).

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{g(r)}. 
}
\author{
  \wei and
  \adrian.
}
\references{
  Ang, Q.W. (2010) Statistical methodology for spatial point patterns
  on a linear network. MSc thesis, University of Western Australia.

  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  Okabe, A. and Yamada, I. (2001) The K-function method on a network and
  its computational implementation. \emph{Geographical Analysis}
  \bold{33}, 271-290.
}

\seealso{
  \code{\link{linearpcf}},
  \code{\link{linearKinhom}},
  \code{\link{lpp}}
}
\examples{
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X ~x)
  g <- linearpcfinhom(X, lambda=fit, update=FALSE)
  plot(g)
  ge <- linearpcfinhom(X, sigma=bw.lppl)
}
\keyword{spatial}
\keyword{nonparametric}


\concept{Linear network}
