\name{density.linnet}
\alias{density.linnet}
\title{Kernel Smoothing of Linear Network}
\description{
  Compute a kernel smoothed intensity function for
  the line segments of a linear network.
}
\usage{
  \method{density}{linnet}(x, \dots)
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"})
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.core]{density.psp}}
    to control the amount of smoothing and the spatial resolution
    of the result.
  }
}
\value{
  A pixel image in two dimensions (object of class \code{"im"})
  or a numeric vector.
}
\details{
  This is the method for the generic function \code{\link[stats]{density}}
  for the class \code{"linnet"} (linear networks).

  The network \code{x} is first converted to a line segment pattern
  (object of class \code{"psp"}). Then the method
  \code{\link[spatstat.core]{density.psp}} is applied to the
  segment pattern.

  A kernel estimate of the intensity of the line segment pattern
  is computed. The result is 
  the convolution of the isotropic Gaussian kernel, of
  standard deviation \code{sigma}, with the line segments.

  The intensity of a line segment pattern
  is the (spatially-varying) amount of segment length per unit area,
  expressed in the same units as the coordinates of \code{x}.
  If the units of \code{x} are in metres, then an intensity value of 3
  means that there are 3 metres of segment length per square metre of
  spatial domain.
  
  See \code{\link[spatstat.core]{density.psp}}
  for more details.
}
\seealso{
  \code{\link[spatstat.core]{density.psp}},
  \code{\link[spatstat.core]{im.object}},
  \code{\link[stats]{density}}.
}
\examples{
  D <- density(simplenet, 0.1)
  plot(D)
  plot(simplenet, add=TRUE, col="white")
  ## compare with average intensity
  volume(simplenet)/area(Window(simplenet))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
