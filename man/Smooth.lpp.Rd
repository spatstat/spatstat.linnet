\name{Smooth.lpp}
\alias{Smooth.lpp}
\title{Spatial Smoothing of Observations on a Network}
\description{
  Performs spatial smoothing of numeric values observed
  at a set of locations on a network. Uses kernel smoothing.
}
\usage{
\method{Smooth}{lpp}(X, sigma,
                     ...,
                     at=c("pixels", "points"),
                     weights=rep(1, npoints(X)),
                     leaveoneout=TRUE) 
}
\arguments{
  \item{X}{
    A marked point pattern on a linear network
    (object of class \code{"lpp"}).
  }
  \item{sigma}{
    Smoothing bandwidth.
    A single positive number.
    See \code{\link{density.lpp}}.
  }
  \item{\dots}{
    Further arguments passed to
    \code{\link{density.lpp}}
    to control the kernel smoothing and
    the pixel resolution of the result.
  }
  \item{at}{
    String specifying whether to compute the smoothed values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
  \item{weights}{
    Optional weights attached to the observations.
    \pointweightsVector.
  }
  \item{leaveoneout}{
    Logical value indicating whether to compute a leave-one-out
    estimator. Applicable only when \code{at="points"}.
  }
}
\details{
  The function \code{Smooth.lpp}
  performs spatial smoothing of numeric values
  observed at a set of irregular locations on a linear network.
  
  \code{Smooth.lpp} is a method for the generic function
  \code{\link[spatstat.explore]{Smooth}} for the class \code{"lpp"} of point patterns.
  Thus you can type simply \code{Smooth(X)}.
  
  Smoothing is performed by kernel weighting, using the Gaussian kernel
  by default. If the observed values are \eqn{v_1,\ldots,v_n}{v[1],...,v[n]}
  at locations \eqn{x_1,\ldots,x_n}{x[1],...,x[n]} respectively,
  then the smoothed value at a location \eqn{u} is
  \deqn{
    g(u) = \frac{\sum_i k(u, x_i) v_i}{\sum_i k(u, x_i)}
  }{
    g(u) = (sum of k(u, x[i]) v[i])/(sum of k(u, x[i]))
  }
  where \eqn{k} is the kernel.
  This is known as the Nadaraya-Watson smoother
  (Nadaraya, 1964, 1989; Watson, 1964).
  The type of kernel is determined by further arguments \code{...}
  which are passed to \code{\link{density.lpp}}

  The argument \code{X} must be a marked point pattern on a linear
  network (object of class \code{"lpp"}).
  The points of the pattern are taken to be the
  observation locations \eqn{x_i}{x[i]}, and the marks of the pattern
  are taken to be the numeric values \eqn{v_i}{v[i]} observed at these
  locations.

  The marks are allowed to be a data frame.
  Then the smoothing procedure is applied to each column of marks. 
  
  The numerator and denominator are computed by \code{\link{density.lpp}}.
  The arguments \code{...} control the smoothing kernel parameters.

  The optional argument \code{weights} allows numerical weights to
  be applied to the data. If a weight \eqn{w_i}{w[i]}
  is associated with location \eqn{x_i}{x[i]}, then the smoothed
  function is 
  (ignoring edge corrections)
  \deqn{
    g(u) = \frac{\sum_i k(u, x_i) v_i w_i}{\sum_i k(u, x_i) w_i}
  }{
    g(u) = (sum of k(u, x[i]) v[i] w[i])/(sum of k(u, x[i]) w[i])
  }
}
\section{Very small bandwidth}{
  If the chosen bandwidth \code{sigma} is very small,
  kernel smoothing is mathematically equivalent
  to nearest-neighbour interpolation.
}
\value{
  \emph{If \code{X} has a single column of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is
    a pixel image on the network (object of class \code{"linim"}). 
    Pixel values are values of the interpolated function.
    \item
    If \code{at="points"}, the result is a numeric vector
    of length equal to the number of points in \code{X}.
    Entries are values of the interpolated function at the points of \code{X}.
  }
  \emph{If \code{X} has a data frame of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is a named list of 
    pixel images on the network (objects of class \code{"linim"}). There is one
    image for each column of marks. This list also belongs to
    the class \code{"solist"}, for which there is a plot method.
    \item
    If \code{at="points"}, the result is a data frame
    with one row for each point of \code{X},
    and one column for each column of marks. 
    Entries are values of the interpolated function at the points of \code{X}.
  }
  The return value has attribute
  \code{"sigma"} which reports the smoothing
  bandwidth that was used.
}
\seealso{
  \code{\link[spatstat.explore]{Smooth}},
  \code{\link{density.lpp}}.
}
\examples{
  X <- spiders
  if(!interactive()) X <- X[owin(c(0,1100), c(0, 500))]
  marks(X) <- coords(X)$x
  plot(Smooth(X, 50))
  Smooth(X, 50, at="points")
}
\author{
  \spatstatAuthors.
}
\references{
  Nadaraya, E.A. (1964) On estimating regression.
  \emph{Theory of Probability and its Applications}
  \bold{9}, 141--142.

  Nadaraya, E.A. (1989) 
  \emph{Nonparametric estimation of probability densities
  and regression curves}.
  Kluwer, Dordrecht.

  Watson, G.S. (1964)
  Smooth regression analysis.
  \emph{Sankhya A} \bold{26}, 359--372.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
\concept{Linear network}
