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
                               do.plot=TRUE,
                               show.vertices=FALSE, show.window=FALSE,
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
  \item{do.plot}{
    Logical value specifying whether to actually perform the plot.
  }
  \item{show.vertices}{
    Logical value specifying whether to plot the vertices as well.
  }
  \item{show.window}{
    Logical value specifying whether to plot the window containing the
    linear network. 
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
  If \code{show.vertices=TRUE}, the vertices of the network will also be
  plotted, using \code{\link[spatstat.geom]{plot.ppp}}.
  If \code{show.window=TRUE}, the window surrounding the network will also
  be plotted.

  If the vertices or line segments of \code{x} are marked,
  the marks are not displayed by default. To plot the marks,
  set \code{use.marks=TRUE}. To plot the marks and plot the associated legends,
  set \code{use.marks=TRUE, legend=TRUE}. To plot only the marks of the
  segments and not the marks of the vertices, set
  \code{args.segments=list(use.marks=TRUE)} and so on.
}
\value{
  An (invisible) list with two elements, \code{segments} and
  \code{vertices} describing the representation of the marks.
  The element \code{segments} contains
  the result of \code{\link[spatstat.geom]{plot.psp}}
  (either a \code{colourmap}, a numeric value or an \code{owin}).
  The element \code{vertices} contains the result
  of \code{\link[spatstat.geom]{plot.ppp}} (a \code{symbolmap})
  or \code{NULL}.

  The result also has attribute \code{"bbox"}
  giving the bounding box for the plot.
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

   L <- simplenet
   marks(L, "vertices") <- letters[1:nvertices(L)]
   marks(L, "segments") <- runif(nsegments(L))
   plot(L, show.vertices=TRUE, use.marks=TRUE, legend=TRUE)
}
\keyword{spatial}
