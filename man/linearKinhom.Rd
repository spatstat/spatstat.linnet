\name{linearKinhom}
\alias{linearKinhom}
\title{
  Inhomogeneous Linear K Function
}
\description{
  Computes an estimate of the inhomogeneous linear \eqn{K} function
  for a point pattern on a linear network.
}
\usage{
linearKinhom(X, lambda=NULL, r=NULL, ..., correction="Ang",
             normalise=TRUE, normpower=1,
	     update=TRUE, leaveoneout=TRUE, sigma=NULL, ratio=FALSE)
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
    There is a sensible default. Users are advised not to specify
    \code{r} in normal usage.
  }
  \item{\dots}{
    Ignored.
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
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    the estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\details{
  This command computes the inhomogeneous version of the 
  linear \eqn{K} function from point pattern data on a linear network.

  The argument \code{lambda} should provide estimated values
  of the intensity of the point process at each point of \code{X}.

  If \code{lambda=NULL}, the intensity will be estimated by kernel
  smoothing by calling \code{\link{density.lpp}} with the smoothing
  bandwidth \code{sigma}, and with any other relevant arguments
  that might be present in \code{\dots}. A leave-one-out kernel estimate
  will be computed if \code{leaveoneout=TRUE}.

  If \code{lambda} is given, it may be a numeric vector (of length equal to
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

  Each estimate is initially computed as 
  \deqn{
    \widehat K_{\rm inhom}(r) = \frac{1}{\mbox{length}(L)}
    \sum_i \sum_j \frac{1\{d_{ij} \le r\}
      e(x_i,x_j)}{\lambda(x_i)\lambda(x_j)}
  }{
    K^inhom(r)= (1/length(L)) sum[i] sum[j] 1(d[i,j] <= r) * 
    e(x[i],x[j])/(lambda(x[i]) * lambda(x[j]))
  }
  where \code{L} is the linear network,
  \eqn{d_{ij}}{d[i,j]} is the distance between points
  \eqn{x_i}{x[i]} and \eqn{x_j}{x[j]}, and
  \eqn{e(x_i,x_j)}{e(x[i],x[j])} is a weight.
  If \code{correction="none"} then this weight is equal to 1,
  while if  \code{correction="Ang"} the weight is
  \eqn{e(x_i,x_j,r) = 1/m(x_i, d_{ij})}{e(x[i],x[j],r) = 1/m(x[i],d[i,j])}
  where \eqn{m(u,t)} is the number of locations on the network that lie
  exactly \eqn{t} units distant from location \eqn{u} by the shortest
  path.

  If \code{normalise=TRUE} (the default), then the estimates
  described above
  are multiplied by \eqn{c^{\mbox{normpower}}}{c^normpower} where 
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
  Older versions of \code{\link{linearKinhom}} interpreted
  \code{lambda=NULL} to mean that the homogeneous function
  \code{\link{linearK}} should be computed. This was changed to the
  current behaviour in version \code{3.1-0} of \pkg{spatstat.linnet}.
}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  \wei and
  \adrian
}
\references{
  Ang, Q.W. (2010) Statistical methodology for spatial point patterns
  on a linear network. MSc thesis, University of Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.
}
\seealso{
  \code{\link{lpp}}
}
\examples{
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X ~x)
  K <- linearKinhom(X, lambda=fit)
  plot(K)
  Ke <- linearKinhom(X, sigma=bw.lppl)
  plot(Ke)
}
\keyword{spatial}
\keyword{nonparametric}


\concept{Linear network}
