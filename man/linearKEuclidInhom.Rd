\name{linearKEuclidInhom}
\alias{linearKEuclidInhom}
\title{
    Inhomogeneous Linear K Function Based on Euclidean Distances
}
\description{
  Computes an estimate of the inhomogeneous linear \eqn{K} function
  based on Euclidean distances, for a point pattern on a linear network.
}
\usage{
   linearKEuclidInhom(X, lambda = NULL, r = NULL, \dots,
      normalise = TRUE, normpower = 2, update = TRUE,
      leaveoneout = TRUE, sigma=NULL)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Intensity values for the point pattern. Either a numeric vector,
    a \code{function}, a pixel image
    (object of class \code{"im"} or \code{"linim"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}) or \code{NULL}.
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Ignored.
  }
  \item{normalise}{
    Logical. If \code{TRUE} (the default), the denominator of the estimator is 
    data-dependent (equal to the sum of the reciprocal intensities at the data
    points, raised to \code{normpower}), which reduces the sampling variability.
    If \code{FALSE}, the denominator is the length of the network.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See Details.
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
}
\details{
  This command computes the inhomogeneous version of the 
  linear \eqn{K} function based on \emph{Euclidean} distances,
  for a point pattern on a linear network.

  This is different from the inhomogeneous \eqn{K} function based on
  \emph{shortest-path} distances,
  which is computed by \code{\link{linearKinhom}}.

  The inhomogeneous \eqn{K} function based on \emph{Euclidean} distances
  is defined in equation (23) of Rakshit, Nair and Baddeley (2017).
  Estimation is performed as described in equation (28).

  The argument \code{lambda} should provide estimated values
  of the intensity of the point process at each point of \code{X}.

  If \code{lambda=NULL}, the intensity will be estimated by kernel
  smoothing by calling \code{\link{density.lpp}} with the smoothing
  bandwidth \code{sigma}, and with any other relevant arguments
  that might be present in \code{\dots}. A leave-one-out kernel estimate
  will be computed if \code{leaveoneout=TRUE}.

  If \code{lambda} is given, then it is expected to provide estimated values
  of the intensity of the point process at each point of \code{X}. 
  The argument \code{lambda} may be a numeric vector (of length equal to
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
  
  If \code{normalise=TRUE} (the default), then the estimate
  is multiplied by \eqn{c^{\mbox{normpower}}}{c^normpower} where 
  \eqn{
    c = \mbox{length}(L)/\sum (1/\lambda(x_i)).
  }{
    c = length(L)/sum[i] (1/lambda(x[i])).
  }
  This rescaling reduces the variability and bias of the estimate
  in small samples and in cases of very strong inhomogeneity.
  The default value of \code{normpower} is 1 (for consistency with
  previous versions of \pkg{spatstat})
  but the most sensible value is 2, which would correspond to rescaling
  the \code{lambda} values so that
  \eqn{
    \sum (1/\lambda(x_i)) = \mbox{area}(W).
  }{
    sum[i] (1/lambda(x[i])) = area(W).
  }
}
\section{Warning}{
  Older versions of \code{\link{linearKEuclidInhom}} interpreted
  \code{lambda=NULL} to mean that the homogeneous function
  \code{\link{linearKEuclid}} should be computed. This was changed to the
  current behaviour in version \code{3.1-0} of \pkg{spatstat.linnet}.
}
\value{
  Function value table (object of class \code{"fv"}).
}
\references{
  Rakshit. S., Nair, G. and Baddeley, A. (2017)
  Second-order analysis of point patterns on a network
  using any distance metric. \emph{Spatial Statistics} \bold{22} (1) 129--154.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{linearpcfEuclidInhom}}, \code{\link{linearKEuclid}}.

  See \code{\link{linearKinhom}} for the corresponding function
  based on shortest-path distances.
}
\examples{
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X ~x)
  K <- linearKEuclidInhom(X, lambda=fit)
  plot(K)
}
\keyword{spatial}
\keyword{nonparametric}
\concept{Linear network}
