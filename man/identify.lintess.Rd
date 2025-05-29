\name{identify.lintess}
\alias{identify.lintess}
\title{Interactively Identify Tiles of a Tessellation on a Network}
\description{
  If a tessellation on a linear network is plotted in the graphics window,
  then each time the left mouse button is pressed,
  this function will find the tile which contains
  the mouse position, print the serial number of the tile, and draw the tile
  in a different colour. 
}
\usage{
  \method{identify}{lintess}(x, \dots, labels=seq_len(nobjects(x)),
                             n=nobjects(x), plot=TRUE)
}
\arguments{
  \item{x}{
    A tessellation on a linear network
    (object of class \code{"lintess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{identify.default}}.
  }
  \item{labels}{
    Labels associated with the tiles of the tessellation, to be plotted when the
    tiles are identified. A character vector or numeric vector
    of length equal to the number of tiles of \code{x}.
  }
  \item{n}{
    Maximum number of tiles to be identified.
  }
  \item{plot}{
    Logical. Whether to plot the \code{labels} when a tile is identified.
  }
}
\value{
  A vector containing the serial numbers of the tiles of
  \code{x} that were identified.
}
\details{
  This is a method for the generic function \code{\link[graphics]{identify}}
  for tessellations on a linear network.

  The tessellation \code{x}, or the underlying network,
  should first be plotted using \code{\link{plot.lintess}}
  or \code{\link{plot.linnet}}.
  
  Then \code{identify(x)}
  reads the position of the graphics pointer each time the
  left mouse button is pressed.  It then determines which tile of
  \code{x} contains the mouse position.
  The index of this tile 
  will be returned as part of the value of the call.

  Each time a tile is identified,
  text will be displayed alongside the tile
  showing its serial number, and the tile will be re-drawn in a
  different colour.
}
\seealso{
  \code{\link[graphics]{identify}},
  \code{\link{plot.lintess}}
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{iplot}
