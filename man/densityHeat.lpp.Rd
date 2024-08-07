\name{densityHeat.lpp}
\alias{densityHeat.lpp}
\title{
  Kernel Density on a Network using Heat Equation
}
\description{
  Given a point pattern on a linear network,
  compute a kernel estimate of intensity,
  by solving the heat equation.
}
\usage{
  \method{densityHeat}{lpp}(x, sigma=NULL, \dots,
              at=c("pixels", "points"), leaveoneout=TRUE,
              weights = NULL,
              dx = NULL, dt = NULL, iterMax = 1e+06,
              finespacing = TRUE, verbose=FALSE)
}
\arguments{
  \item{x}{
    Point pattern on a linear network (object of class \code{"lpp"})
    to be smoothed.
  }
  \item{sigma}{
    Smoothing bandwidth (standard deviation of the kernel).
    A numeric value in the same units as the spatial coordinates of \code{x}.
    Alternatively \code{sigma} may be a function which selects a
    bandwidth when applied to \code{X},
    for example, \code{\link[spatstat.explore]{bw.scott.iso}} or \code{\link{bw.lppl}}.
    There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} determining the
    resolution of the result. (Any other arguments are ignored.)
  }
  \item{at}{
    String specifying whether to compute the intensity values
    at a fine grid of pixel locations on the network
    (\code{at="pixels"}, the default) or
    only at the data points of \code{x} (\code{at="points"}).
  }
  \item{leaveoneout}{
    Logical value indicating whether to compute a leave-one-out
    estimator. Applicable only when \code{at="points"}.
  }
  \item{weights}{
    Optional. Numeric vector of weights associated with the
    points of \code{x}. Weights may be positive, negative or zero.
  }
  \item{dx}{
    Optional. Spacing of the sampling points along the network.
    A single number giving a distance value in the same units as \code{x}.
  }
  \item{dt}{
    Optional. Time step in the heat equation solver.
    A single number. 
  }
  \item{iterMax}{
    Maximum number of iterations.
  }
  \item{finespacing}{
    Logical value specifying whether the discrete approximation is
    required to be accurate along every segment of the network,
    no matter how short the segment is. See the section on Discretisation.
  }
  \item{verbose}{
    Logical value specifying whether to print progress reports.
  }
}
\details{
  The function \code{\link[spatstat.explore]{densityHeat}} is generic.
  This is the method for the class \code{"lpp"} of points on a linear
  network.
  
  Kernel smoothing is applied to the points of \code{x}
  using a kernel based on path distances in the network.
  If \code{at="pixels"} (the default),
  the result is a pixel image on the linear network (class
  \code{"linim"}) which can be plotted.
  If \code{at="points"} the result is a numeric vector giving the
  density estimates at the data points of \code{x}. 

  The smoothing operation is equivalent to the
  \dQuote{equal-split continuous} rule described in
  Section 9.2.3 of Okabe and Sugihara (2012).
  However, the actual computation is performed rapidly, by solving the classical
  time-dependent heat equation on the network,
  as described in McSwiggan et al (2016).
  Computational time is short, but increases quadratically with
  \code{sigma}.

  If \code{at="points"} and \code{leaveoneout=TRUE}, 
  a leave-one-out estimate is computed at each data point
  (that is, the estimate at each data point \code{x[i]} is based
  on all of the points except \code{x[i]})
  using the truncated series approximation
  of McSwiggan et al (2019).

  The argument \code{sigma} specifies the smoothing bandwidth.
  If \code{sigma} is missing or \code{NULL},
  the default is one-eighth of the length of the shortest side
  of the bounding box of \code{x}.
  If \code{sigma} is a function in the \R language, it is assumed
  to be a bandwidth selection rule, and it will be applied to \code{x}
  to compute the bandwidth value.
}
\section{Infinite bandwidth}{
  If \code{sigma=Inf}, the resulting density estimate is 
  constant over all locations,
  and is equal to the average density of points per unit length.
  (If the network is not connected, then this rule
  is applied separately to each connected component of the network).
}
\section{Discretisation and Error Messages}{
  The arguments \code{dx}, \code{dt} and \code{iterMax}
  determine the discretisation of the network, according to
  a set of rules. The argument \code{finespacing} determines
  which rule will be applied. 

  The arguments \code{dx}, \code{dt}, \code{iterMax}
  are connected by several constraints; specifying one of these
  arguments will affect the default values of the other two arguments.
  
  The argument \code{finespacing} specifies whether
  a very fine spacing of sample points is required, in order to attain
  hihg accuracy.
  \itemize{
    \item
    If \code{finespacing=TRUE} (the default), then the
    sample point spacing \code{dx} must not exceed
    one-third of the length of the shortest
    segment of the network. This ensures that the discrete approximation
    is accurate along every segment, no matter how short the segment
    is. However, this may not be feasible
    if it implies a very large number of sample points, or a large number
    of iterations: in such cases, the code may terminate with an error
    about illegal values of \code{dx}, \code{dt} or \code{iterMax}.
    \item
    If \code{finespacing=FALSE}, then
    the sample point spacing \code{dx} will be about one-half the width
    of a pixel in the default pixellation of the window of \code{x}.
    This is usually a much coarser resolution than the one
    selected by \code{finespacing=TRUE}. If it is too coarse,
    the pixel resolution can be refined using
    the arguments \code{dimyx}, \code{eps} or \code{xy}
    passed to \code{\link[spatstat.geom]{as.mask}}. For example,
    \code{dimyx=512} would specify a 512 x 512 pixel grid.
    The default pixel resolution can be changed for the remainder of
    the \R session by \code{\link[spatstat.geom]{spatstat.options}('npixel')}.
  }
}
\value{
  If \code{at="pixels"} (the default),
  a pixel image on the linear network (object of class \code{"linim"}).

  If \code{at="points"}, a numeric vector with one entry for each point
  of \code{x}.
}
\references{
  McSwiggan, G., Baddeley, A. and Nair, G. (2016)
  Kernel density estimation on a linear network.
  \emph{Scandinavian Journal of Statistics} \bold{44}, 324--345.

  McSwiggan, G., Baddeley, A. and Nair, G. (2019)
  Estimation of relative risk for events on a linear network.
  \emph{Statistics and Computing} \bold{30}, 469--484.

  Okabe, A. and Sugihara, K. (2012)
  \emph{Spatial analysis along networks}.
  Wiley.
}
\author{
  \adrian and Greg McSwiggan.
}
\seealso{
  \code{\link{density.lpp}}
}
\examples{
  X <- runiflpp(3, simplenet)
  D <- densityHeat(X, 0.2)
  plot(D, style="w", main="", adjust=2)
  densityHeat.lpp(X, 0.2, at="points")
  Dw <- densityHeat(X, 0.2, weights=c(1,2,-1))
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Linear network}
