\name{linearpcfEuclidInhom}
\alias{linearpcfEuclidInhom}
\title{
    Inhomogeneous Linear Pair Correlation Function Based on Euclidean Distances
}
\description{
  Computes an estimate of the inhomogeneous pair correlation function
  based on Euclidean distances, for a point pattern on a linear network.
}
\usage{
   linearpcfEuclidInhom(X, lambda = NULL, r = NULL, \dots,
      normalise = TRUE, normpower = 2, update = TRUE, leaveoneout = TRUE)
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
    or \code{"lppm"}).
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
    Logical value (passed to \code{\link{fitted.lppm}} or
    \code{\link{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity,
    when \code{lambda} is a fitted model.
    Supported only when \code{update=TRUE}.
  }
}
\details{
  This command computes the inhomogeneous version of the 
  pair correlation function based on \emph{Euclidean} distances,
  for a point pattern on a linear network.

  This is different from the inhomogeneous pair correlation function based on
  \emph{shortest-path} distances, which is computed by \code{\link{linearpcfinhom}}.

  If \code{lambda = NULL} the result is equivalent to the
  homogeneous pair correlation function \code{\link{linearpcfEuclid}}.
  If \code{lambda} is given, then it is expected to provide estimated values
  of the intensity of the point process at each point of \code{X}. 
  The argument \code{lambda} may be a numeric vector (of length equal to
  the number of points in \code{X}), or a \code{function(x,y)} that will be
  evaluated at the points of \code{X} to yield numeric values, 
  or a pixel image (object of class \code{"im"}) or a fitted point 
  process model (object of class \code{"ppm"} or \code{"lppm"}).

  The inhomogeneous pair correlation function based on \emph{Euclidean} distances
  is defined in equation (30) of Rakshit, Nair and Baddeley (2017).
  Estimation is performed as described in equation (34) of 
  Rakshit, Nair and Baddeley (2017).
  
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
  \code{\link{linearKEuclidInhom}}, \code{\link{linearpcfEuclid}}.

  See \code{\link{linearpcfinhom}} for the corresponding function
  based on shortest-path distances.
}
\examples{
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X ~x)
  g <- linearpcfEuclidInhom(X, lambda=fit)
  plot(g)
}
\keyword{spatial}
\keyword{nonparametric}
