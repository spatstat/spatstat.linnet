\name{mergeTiles.lintess}
\alias{mergeTiles.lintess}
\title{
  Merge Groups of Tiles in a Tessellation on a Network
}
\description{
  Given a tessellation on a linear network,
  and a factor which classifies the tiles of the
  tessellation into groups, form the set union of each group of tiles,
  and construct the tessellation based on these merged tiles.
}
\usage{
 \method{mergeTiles}{lintess}(x, group, \dots)
}
\arguments{
  \item{x}{
    A tessellation on a network (object of class \code{"lintess"}).
  }
  \item{group}{
    A \code{factor}, with one entry for each tile of \code{x},
    which classifies the tiles into groups.
  }
  \item{\dots}{Other arguments, currently ignored.}
}
\value{
  A tessellation on a network (object of class \code{"lintess"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.geom]{mergeTiles}}
}
\examples{
  A <- as.lintess(simplenet)
  A
  g <- sample(letters[1:3], nobjects(A), replace=TRUE)
  g
  mergeTiles(A, g) 
}
\keyword{spatial}
\keyword{manip}
