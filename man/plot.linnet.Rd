\name{plot.linnet}
\alias{plot.linnet}
\title{
  Plot a linear network
}
\description{
  Plots a linear network
}
\usage{
 \method{plot}{linnet}(x, ..., main=NULL, add=FALSE, 
                               vertices=FALSE, window=FALSE,
                               do.plot=TRUE,
                               args.vertices=list(), args.segments=list())
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Graphics arguments passed to \code{\link{plot.psp}}
    and \code{\link{plot.ppp}}.
  }
  \item{main}{
    Main title for plot. Use \code{main=""} to suppress it.
  }
  \item{add}{
    Logical. If code{TRUE}, superimpose the graphics
    over the current plot. If \code{FALSE}, generate a new plot.
  }
  \item{vertices}{
    Logical. Whether to plot the vertices as well.
  }
  \item{window}{
    Logical. Whether to plot the window containing the linear network.
  }
  \item{do.plot}{
    Logical. Whether to actually perform the plot.
  }
  \item{args.segments}{
    Optional list of arguments passed to \code{\link[spatstat.geom]{plot.psp}}
    when plotting the line segments of the network.
    These arguments override any arguments in \code{\dots}.
  }
  \item{args.vertices}{
    Optional list of arguments passed to \code{\link[spatstat.geom]{plot.ppp}}
    when plotting the vertices of the network (only when \code{vertices=TRUE}).
    These arguments override any arguments in \code{\dots}.
  }
}
\details{
  This is the plot method for class \code{"linnet"}.

  The line segments of the network \code{x} are plotted
  using \code{\link[spatstat.geom]{plot.psp}}. 
  If \code{vertices=TRUE}, the vertices of the network will also be
  plotted, using \code{\link[spatstat.geom]{plot.ppp}}.
  If \code{window=TRUE}, the window surrounding the network will also
  be plotted.

  By default, marks attached to the vertices or the line segments
  are not displayed. Use \code{args.segments} or \code{args.vertices}
  to change the default.  
}
\value{
  An (invisible) object of class \code{"owin"}
  giving the bounding box of the network.
}
\author{
  \wei and
  \adrian 
}
\seealso{
   \code{\link{linnet}}
}
\examples{
   plot(simplenet)
}
\keyword{spatial}
