\name{tilenames.lintess}
\alias{tilenames.lintess}
\alias{tilenames<-.lintess}
\title{Names of Tiles in a Tessellation on a Network}
\description{
  Extract or Change the Names of the Tiles in a Tessellation on a Network.
}
\usage{
\method{tilenames}{lintess}(x)

\method{tilenames}{lintess}(x) <- value
}
\arguments{
  \item{x}{
    A tessellation 
    on a linear network (object of class \code{"lintess"}).
  }
  \item{value}{Character vector giving new names for the tiles.}
}
\details{
  These functions extract or change the names of the 
  tiles that make up the tessellation \code{x}.

  If the tessellation is a regular grid, the tile names
  cannot be changed.
}
\value{
  \code{tilenames} returns a character vector.
}
\seealso{
  \code{\link{lintess}},
  \code{\link[spatstat.geom]{tiles}}
}
\examples{
  B <- lineardirichlet(runiflpp(5, simplenet))
  tilenames(B)
  tilenames(B) <- letters[1:5]
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
