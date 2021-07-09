\name{persp.linfun}
\alias{persp.linfun}
\title{
  Perspective View of Function on a Linear Network
}
\description{
  Given a function on a linear network, generate a
  perspective view.
}
\usage{
 \method{persp}{linfun}(x, \dots, main, eps = NULL, dimyx = NULL, xy = NULL)
}
\arguments{
  \item{x}{
    The function to be plotted. An object of class \code{"linfun"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{persp.linim}}
    controlling the appearance of the plot.
  }
  \item{main}{
    Main title for the plot.
  }
  \item{eps,dimyx,xy}{
    Arguments passed to \code{\link[spatstat.linnet]{as.linim}} determining the
    spatial resolution when the function is converted to an image.
  }
}
\details{
  The function \code{x} is converted to a pixel image on the linear
  network using \code{\link[spatstat.linnet]{as.linim}}. Then
  \code{\link{persp.linim}} is invoked to generate the perspective plot.

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
  \code{\link{persp.linim}}
}
\examples{
f <- linfun(function(x,y,seg,tp) { abs(sin(25*x)) + abs(sin(15*y)) }, simplenet)
persp(f, phi=20)
}
\references{
  Okabe, A. and Sugihara, K. (2012)
  \emph{Spatial Analysis Along Networks}. John Wiley and Sons, New York.
}
\keyword{spatial}
\keyword{hplot}

