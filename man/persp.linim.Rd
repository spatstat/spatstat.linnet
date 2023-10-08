\name{persp.linim}
\alias{persp.linim}
\title{
  Perspective View of Pixel Image on a Linear Network
}
\description{
  Given a pixel image on a linear network, generate a perspective view.
}
\usage{
 \method{persp}{linim}(x, \dots, main,
                       grid = TRUE, ngrid = 10,
                       col.grid = "grey", col.base = "white",
                       neg.args=list(), warncross=FALSE)
}
\arguments{
  \item{x}{
    Pixel image on a linear network (object of class \code{"linim"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{persp.default}}
    to control the perspective view, or passed to 
    \code{\link[graphics]{segments}} or \code{\link[graphics]{polygon}}
    to control the appearance of the vertical panes.
  }
  \item{main}{
    Main title for the plot.
  }
  \item{grid}{
    Logical value indicating whether to draw a rectangular grid
    at height zero, to assist the perception of perspective.
  }
  \item{ngrid}{
    Number of grid lines to draw, if \code{grid=TRUE}.
  }
  \item{col.grid}{
    Colour of grid lines, if \code{grid=TRUE}.
  }
  \item{col.base}{
    Colour of base plane, if \code{grid=TRUE}.
  }
  \item{neg.args}{
    Optional list of arguments passed to \code{\link[graphics]{polygon}}
    when displaying negative values of the function.
  }
  \item{warncross}{
    Logical value indicating whether to issue a warning if
    two segments of the network cross each other (which causes
    difficulty for the algorithm).
  }
}
\details{
  The pixel values are interpreted as the spatially-varying height of a
  vertical surface erected on each segment of the linear network. These
  surfaces are drawn in perspective view.

  This style of plot is often attributed to Okabe and Sugihara (2012).
}
\value{
  (Invisibly) the perspective transformation matrix,
  as described in the help for \code{\link[graphics]{persp.default}}.
}
\author{
  \adrian and Greg McSwiggan.
}
\seealso{
  \code{\link{persp.linfun}}
}
\references{
  Okabe, A. and Sugihara, K. (2012)
  \emph{Spatial Analysis Along Networks}. John Wiley and Sons, New York.
}
\examples{
  if(interactive()) {
    Z <- density(chicago, 100)
  } else {
    X <- runiflpp(10, simplenet)
    Z <- density(X, 0.1)
  }
  persp(Z, theta=30, phi=20)
}
\keyword{spatial}
\keyword{hplot}
\concept{Linear network}
