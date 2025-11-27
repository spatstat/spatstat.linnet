\name{as.linfun} 
\alias{as.linfun}
\alias{as.linfun.linim}
\alias{as.linfun.linnet}
\alias{as.linfun.lintess}
\alias{as.function.linim}
\alias{as.function.linnet}
\alias{as.function.lintess}
\title{
  Convert Data to a Function on a Linear Network
}
\description{
  Convert some kind of data to an object of class \code{"linfun"}
  representing a function on a linear network.
}
\usage{
  as.linfun(X, \dots)

  \method{as.linfun}{linim}(X, \dots)

  \method{as.linfun}{linnet}(X, \dots, values=marks(X))

  \method{as.linfun}{lintess}(X, \dots, values=marks(X), navalue=NA)

  \method{as.function}{linim}(x, \dots)

  \method{as.function}{linnet}(x, \dots)

  \method{as.function}{lintess}(x, \dots)
}
\arguments{
  \item{X,x}{
    Some kind of data to be converted.
  }
  \item{\dots}{
    Other arguments passed to methods.
  }
  \item{values}{
    Optional. Vector of function values,
    one entry associated with each tile of the tessellation.
  }
  \item{navalue}{
    Optional. Function value associated with locations that
    do not belong to a tile of the tessellation.
  }
}
\details{
  An object of class \code{"linfun"} represents
  a function defined on a linear network.
  This page documents methods for converting other kinds of objects
  to a \code{"linfun"} object.

  The function \code{as.linfun} is generic. There are methods for
  converting linear networks (class \code{"linnet"}),
  pixel images on a linear network (class \code{"linnet"}) and
  tessellations on a linear network (class \code{"lintess"})
  to functions on a network.
  Equivalent methods are also provided for the
  generic \code{\link[base]{as.function}}.
  
  The methods \code{as.linfun.linim} and \code{as.function.linim}
  convert objects of class \code{"linim"} (pixel images on a linear
  network) to functions on the network.

  The methods \code{as.linfun.linnet} and \code{as.function.linnet}
  converts a linear network (object
  of class \code{"linnet"}) to a function on the network. The function
  values are specified by the argument \code{values}.
  It should be a vector with one entry for each segment of the network;
  any point lying on segment number \code{i}
  will return the value \code{values[i]}.
  If \code{values} is missing or \code{NULL}, the function values are
  taken to be the marks attached to the segments (\code{values=marks(X)});
  if there are no marks attached to the segments, the function value is
  the integer index of the segment (\code{values=seq_len(nsegments(X))}).

  The methods \code{as.linfun.lintess} and \code{as.function.lintess}
  convert a tessellation
  on a linear network into a function with a different value
  on each tile of the tessellation.
  The function values are specified by the argument \code{values}.
  It should be a vector with one entry for each tile of the tessellation;
  any point lying in tile number \code{i} will return the value \code{values[i]}.
  If \code{values} is missing, the marks of the tessellation
  are taken as the function values. 
  If \code{values} is missing and the tessellation has no marks,
  or if \code{values} is given as \code{NULL},
  then the function returns factor values identifying which tile
  contains each given point.
}
\value{
  Object of class \code{"linfun"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{linfun}}
}
\examples{
   X <- runiflpp(2, simplenet)
   Y <- runiflpp(5, simplenet)

   # image on network
   D <- density(Y, 0.1)

   f <- as.linfun(D)
   f
   f(X)

   h <- as.linfun(simplenet)

   # tessellation on network
   Z <- lineardirichlet(Y)
   g <- as.linfun(Z)
   g(X)
   h <- as.linfun(Z, values = runif(5))
   h(X)
}
\keyword{spatial}
\keyword{manip}
\concept{Linear network}
