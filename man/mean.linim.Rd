\name{mean.linim}
\alias{mean.linim}
\alias{median.linim}
\alias{quantile.linim}
\alias{quantilefun.linim}
\title{Mean, Median, Quantiles of Pixel Values on a Linear Network}
\description{
  Calculates the mean, median, or quantiles
  of the pixel values in a pixel image on a linear network.
}
\usage{
  \method{mean}{linim}(x, \dots)

  \method{median}{linim}(x, \dots)

  \method{quantile}{linim}(x, probs=seq(0,1,0.25), \dots)

  \method{quantilefun}{linim}(x, \dots, type=1)
}
\arguments{
  \item{x}{
    A pixel image on a linear network (object of class
    \code{"linim"}).
  }
  \item{probs}{
    Vector of probabilities for which quantiles should be
    calculated.
  } 
  \item{\dots}{Arguments passed to other methods.}
  \item{type}{
    Integer specifying the type of quantiles,
    as explained in \code{\link[stats]{quantile.default}}.
    Only types 1 and 2 are currently implemented.
  }
}
\details{
  These functions calculate the mean, median and quantiles
  of the pixel values in the image
  \code{x} on a linear network.

  An object of class \code{"linim"}
  describes a pixel image on a linear network. See \code{\link{linim}}.

  The functions described here are methods for the 
  generic \code{\link{mean}}, \code{\link[stats]{median}}
  and \code{\link[stats]{quantile}} for the class \code{"linim"}.
}
\value{
  For \code{mean} and \code{median}, a single number.
  For \code{quantile}, a numeric vector of the same length as \code{probs}.
  For \code{quantilefun}, a function.
}
\seealso{
  \code{\link{mean}},
  \code{\link[stats]{median}},
  \code{\link[stats]{quantile}},

  \code{\link{mean.im}},
  \code{\link{quantile.im}},
  \code{\link{quantilefun}}
}
\examples{
  M <- psp2mask(as.psp(simplenet))
  Z <- as.im(function(x,y) {x-y}, W=M)
  X <- linim(simplenet, Z)
  X
  mean(X)
  median(X)
  quantile(X)
  f <- quantilefun(X)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
