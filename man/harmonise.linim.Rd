\name{harmonise.linim}
\alias{harmonise.linim}
\alias{harmonize.linim}
\title{Make Pixel Images on a Network Compatible}
\description{
  Convert several pixel images to a common pixel raster.
}
\usage{
\method{harmonise}{linim}(\dots)

\method{harmonize}{linim}(\dots)
}
\arguments{
  \item{\dots}{
    Any number of pixel images on a network
    (objects of class \code{"linim"})
    or data which can be converted to pixel images on a network
    by \code{\link{as.linim}}.
  } 
}
\details{
  This function makes any number of pixel images on a network compatible,
  by converting them all to a common pixel grid.

  The command \code{\link[spatstat.geom]{harmonise}} is generic. This is the
  method for objects of class \code{"linim"}.
  
  At least one of the arguments \code{\dots} must be a pixel image
  on a network (object of class \code{"linim"}) or a network
  (class \code{"linnet"}) so that the network is defined.
  
  If several arguments contain network information then they must
  specify the same network.

  Other arguments may be two-dimensional images
  (class \code{"im"}), windows (class \code{"owin"}),
  functions (\code{function(x,y)}) or numerical constants. These will be
  converted to images using \code{\link{as.linim}}.

  The return value is a list, with entries corresponding to the input
  arguments, in which each entry is a pixel image on the same network.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.
}
\value{
  A list, of length equal to the number of arguments \code{\dots},
  whose entries are pixel images on a network.
}
\author{
  \spatstatAuthors.
}
\examples{
   g <- linfun(function(x,y, seg, tp) { seg }, simplenet)
   Image1 <- as.linim(g)
   Distfun <- distfun(runiflpp(3, simplenet))
   Image2 <- as.im(function(x,y) { x }, Window(simplenet))
   harmonise(Image1, Distfun, Image2)
}
\seealso{
  \code{\link[spatstat.geom]{harmonise}},
  \code{\link[spatstat.geom]{harmonise.im}},
  \code{\link{as.linim}}
}
\keyword{spatial}
\keyword{manip}
