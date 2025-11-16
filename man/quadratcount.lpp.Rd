\name{quadratcount}
\alias{quadratcount.lpp}
\title{Quadrat counting for a point pattern on a linear network}
\description{
  Given a point pattern on a linear network,
  divide the network into tiles,
  and count the numbers of points in each file.
}
\usage{
  \method{quadratcount}{lpp}(X, \dots, nx=5, ny=nx,
               xbreaks=NULL, ybreaks=NULL, left.open=TRUE,
               tess=NULL)
}
\arguments{
  \item{X}{
    A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{nx,ny}{
    Numbers of rectangular quadrats in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{xbreaks} and \code{ybreaks}.
  }
  \item{\dots}{Additional arguments are ignored.}
  \item{xbreaks}{
    Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the rectangular quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the rectangular quadrats. Incompatible with \code{ny}.
  }
  \item{tess}{
    Tessellation (object of class \code{"tess"} or \code{"lintess"})
    determining the quadrats. Incompatible
    with \code{nx,ny,xbreaks,ybreaks}.
  }
  \item{left.open}{
    Logical value specifying whether rectangular quadrats are left-open and
    right-closed (\code{left.open=TRUE}, the default) or
    left-closed and right-open (\code{left.open=FALSE}).
  }
}
\value{
  Contingency table containing the number of points counted in each
  tile. The table is also an object of the
  special class \code{"linearquadratcount"}
  and there is a plot method for this class.
}
\details{
  Quadrat counting is an elementary technique for analysing spatial
  point patterns. See Diggle (2003).

  The function \code{\link[spatstat.geom]{quadratcount}} is generic.
  This page documents the method \code{quadratcount.lpp}
  for the class \code{"lpp"} of point patterns on a linear network.

  First the network is divided into pieces, as described below.
  Then the number of points of \code{X} falling in each piece
  of the network is counted. These numbers are returned as a contingency table.

  \itemize{
    \item
    By default, if none of the arguments
    \code{nx, ny, xbreaks, ybreaks, tess} is given,
    every segment of the network is taken as a separate piece.
    The number of points in each segment of the network is counted.
    \item
    If \code{nx, ny} are given, the window containing
    the point pattern \code{X} is divided into
    an \code{nx * ny} grid of rectangular tiles or `quadrats'.
    These tiles are then intersected with the network on which \code{X}
    is defined. The number of points falling in each rectangle
    is counted.
    \item
    If \code{xbreaks} is given, 
    the window containing the point pattern \code{X} will be divided
    into rectangles, with \code{xbreaks} and \code{ybreaks}
    giving the \eqn{x} and \eqn{y} coordinates of the rectangle
    boundaries, respectively.
    The lengths of \code{xbreaks} and \code{ybreaks} may be different.
    \item 
    The argument \code{tess} can be a tessellation on the network
    (object of class \code{"lintess"}) whose tiles will serve as the
    quadrats.
    \item
    Alternatively \code{tess} can be a two-dimensional tessellation
    (object of class \code{"tess"}) which will be intersected with the
    network to determine the tessellation of the network.
  }
  
  The algorithm counts the number of points of \code{X}
  falling in each tile of the tessellation, and returns these counts as a
  contingency table. 

  The return value is a \code{table} which can be printed neatly.
  The return value is also a member of the special class
  \code{"linearquadratcount"}. Plotting the object will display the
  quadrats, annotated by their counts. See the examples.

  To calculate an estimate of intensity based on the quadrat counts,
  use \code{\link{intensity.linearquadratcount}}.

  To extract the quadrats used in a \code{linearquadratcount} object,
  use \code{\link{as.lintess}}.

  Marks attached to the points are ignored by \code{quadratcount.lpp}.
  To obtain a separate contingency table for each type of point
  in a multitype point pattern,
  first separate the different points using \code{\link[spatstat.geom]{split.ppx}},
  then apply \code{quadratcount.lpp} to each pattern.
}
\section{Treament of data points on the boundary}{
  The treatment of points which lie
  on the boundary of two quadrats is undefined, and may depend on the
  hardware.
}
\seealso{
  \code{\link{lintess}}.
}
\references{
  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Academic Press, 2003.

  Stoyan, D. and Stoyan, H. (1994)
  \emph{Fractals, random shapes and point fields:
  methods of geometrical statistics.}
  John Wiley and Sons.
} 
\examples{
  X <- runiflpp(40, simplenet)
  A <- quadratcount(X)
  A
  plot(A)
  B <- quadratcount(X, nx=2)
  B
  plot(B)
}
\author{\spatstatAuthors.}
\keyword{spatial}
\keyword{math}
 
 
