\name{transect.linim}
\alias{transect.linim}
\title{
  Pixel Values Along a Transect of a Linear Network
}
\description{
  Extract values of a pixel image on a linear network,
  along a path in the network.
}
\usage{
transect.linim(X, \dots, path=NULL, click=FALSE, add=FALSE,
               nsample=512, Xname=NULL)
}
\arguments{
  \item{X}{
    A pixel image on a linear network (object of class \code{"linim"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{path}{
    Integer vector (of length at least 2)
    specifying the path as a sequence of vertices on the network.
  }
  \item{click}{
    Logical value. 
    If \code{TRUE}, the path is determined interactively
    by the user, who selects a sequence vertices of the network
    by clicking on the current plot.
  }
  \item{add}{
    Logical. If \code{click=TRUE}, this argument determines
    whether to perform interactive tasks on the current plot (\code{add=TRUE})
    or to start by plotting \code{X} (\code{add=FALSE}).
  }
  \item{nsample}{
    Integer. Total number of sample locations along the path.
  }
  \item{Xname}{
    Optional character string giving the name of \code{X} to be used
    when plotting the transect.
  }
}
\details{
  The pixel values of the image \code{X} along a path on the network
  will be extracted. The result is a function table (\code{"fv"} object)
  which can be plotted directly.

  The argument \code{path} specifies the path in the network.
  Alternatively, if \code{path=NULL} and \code{click=TRUE}, the user can
  interactively specify the path.
  \itemize{
    \item
    If \code{path} is given, it specifies the path
    as a sequence of vertices of the network, where each vertex is identified
    by its integer index. For example, \code{path=c(1,3,7)} specifies
    a path going from \code{V[1]} to \code{V[3]} to \code{V[7]}, where
    \code{V = vertices(as.linnet(X))}.
    \item
    If \code{path=NULL} and \code{click=TRUE}, the user will be able to
    select successive vertices along the path by clicking on a plot of
    the network.
    (Left-click to select a vertex; right-click or ESC to finish). 
  }
  For each leg of the path,
  if the successive vertices are not joined directly by a segment of the
  network, the shortest path between them will be taken.

  A sequence of \code{nsample} sample points along the path are
  generated, and the pixel value of the image at each sample point is
  extracted. The pixel values are assembled into a function table
  (object of class \code{"fv"}). The default plot of this object is
  a graph of pixel value against distance along the path.
  
  The result has attributes \code{"path"} giving the indices of
  the vertices along the (expanded) path,
  and \code{"lines"} giving the line segments making up the path.
  See the examples.
}
\value{
  An object of class \code{"fv"} which can be plotted.
  The result has attributes \code{"path"} (integer vector)
  and \code{"lines"} (object of class \code{"psp"})
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.explore]{fv.object}},
  \code{\link[spatstat.explore]{transect.im}}
}
\examples{
  if(interactive()) {
    D <- density(chicago, 150)
    A <- transect.linim(D, path=c(1, 248))
  } else {
    Z <- as.im(function(x,y) {x-y}, Frame(simplenet))
    D <- linim(simplenet, Z)
    A <- transect.linim(D, path=c(1, 7))
  }
  opa <- par(mfrow=c(1,2))
  plot(A, main="transect")
  plot(D, main="image")
  plot(attr(A, "lines"), add=TRUE, col=8, lwd=7)
  par(opa)
  attr(A, "path")
}
\keyword{spatial}
\keyword{manip}
\keyword{iplot}
