\name{integral.linim}
\alias{integral.linim}
\alias{integral.linfun}
\title{
  Integral on a Linear Network
}
\description{
  Computes the integral (total value) of a function or pixel image
  over a linear network.
}
\usage{
\method{integral}{linim}(f, domain=NULL, weight=NULL, ...)

\method{integral}{linfun}(f, domain=NULL, weight=NULL, ..., exact=FALSE, delta, nd)
}
\arguments{
  \item{f}{
    A pixel image on a linear network (class \code{"linim"})
    or a function on a linear network (class \code{"linfun"}).
  }
  \item{domain}{
    Optional window specifying the domain of integration.
    Alternatively a tessellation (class \code{"tess"} or \code{"lintess"}).
  }
  \item{weight}{
    Optional numerical weight function for the integration.
    A pixel image (object of class \code{"linim"} or \code{"im"}),
    a function (object of class \code{"linfun"}, \code{"funxy"} or a
    a \code{function(x,y)}) or anything acceptable to \code{\link{as.linim}}.
  }
  \item{\dots}{
    Ignored.
  }
  \item{exact}{
    Logical value specifying whether to use a more accurate
    (and slower) calculation method.
    See Details.
  }
  \item{delta}{
    Optional. 
    The step length (in coordinate units)
    for computing the approximate integral (if \code{exact=FALSE}).
    A single positive number.
  }
  \item{nd}{
    Optional. 
    Integer giving the approximate number of sample points on the
    network (if \code{exact=FALSE}).
  }
}
\details{
  The integral (total value of the function over the network) is calculated.

  If \code{domain} is a window (class \code{"owin"}) then the integration
  will be restricted to this window. If \code{domain} is a tessellation
  (class \code{"tess"} or \code{"lintess"})
  then the integral of \code{f} in each
  tile of \code{domain} will be computed.

  If \code{weight} is given, effectively the integral of \code{weight * f}
  is computed.

  For objects of class \code{"linfun"} there is the option of using
  a more accurate calculation method in which the integral along each
  segment of the network is computed separately using the utility 
  \code{\link[stats]{integrate}} from the \pkg{stats} package.
  If \code{exact=TRUE}, additional arguments \code{\dots} are
  passed to \code{\link[stats]{integrate}} to control the computation.
}
\value{
  A single numeric or complex value
  (or a vector of such values if \code{domain} is a tessellation).
}
\seealso{
  \code{\link{linim}},
  \code{\link[spatstat.geom]{integral.im}},
  \code{\link[stats]{integrate}}
}
\examples{
  # make a function and image data
  xcoord <- linfun(function(x,y,seg,tp) { x }, simplenet)
  integral(xcoord)
  integral(xcoord, exact=TRUE)
  X <- as.linim(xcoord)
  integral(X)

  # integrals inside each tile of a tessellation
  A <- quadrats(Frame(simplenet), 3)
  integral(X, A)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
\concept{Linear network}
