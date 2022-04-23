\name{plot.lpp}
\alias{plot.lpp}
\title{
  Plot Point Pattern on Linear Network
}
\description{
  Plots a point pattern on a linear network.
  Plot method for the class \code{"lpp"} of point patterns on a linear network.
}
\usage{
\method{plot}{lpp}(x, \dots, main, add = FALSE,
     type = c("p", "n"),
     use.marks=TRUE, which.marks=NULL,
     legend=TRUE,
     leg.side=c("left", "bottom", "top", "right"),
     leg.args=list(),
     show.all = !add, show.window=FALSE, show.network=TRUE, 
     do.plot = TRUE, multiplot=TRUE) 
}
\arguments{
  \item{x}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{plot.linnet}}
    or \code{\link{plot.ppp}}.
  }
  \item{main}{
    Main title for plot.
  }
  \item{add}{
    Logical value indicating whether the plot is to be added to the
    existing plot (\code{add=TRUE}) or whether a new plot should be
    initialised (\code{add=FALSE}, the default).
  }
  \item{type}{
    Type of plot: either \code{"p"} or \code{"n"}.
    If \code{type="p"} (the default), both the points and the observation window
    are plotted. If \code{type="n"}, only the window is plotted.
  }
  \item{use.marks}{
    logical flag; if \code{TRUE}, plot points using a different
    plotting symbol for each mark;
    if \code{FALSE}, only the locations of the points will be plotted,
    using \code{\link{points}()}.
  }
  \item{which.marks}{
    Index determining which column of marks to use,
    if the marks of \code{x} are a data frame.
    A character or integer vector identifying one or more
    columns of marks. 
    If \code{add=FALSE} then
    the default is to plot all columns of marks, in a series of separate
    plots.
    If \code{add=TRUE} then only one column of marks can be plotted,
    and the default is \code{which.marks=1}
    indicating the first column of marks.
  }
  \item{legend}{
    Logical value indicating whether to add a legend showing the mapping
    between mark values and graphical symbols (for a marked point pattern).
  }
  \item{leg.side}{
    Position of legend relative to main plot.
  }
  \item{leg.args}{
    List of additional arguments passed to \code{\link{plot.symbolmap}}
    or \code{\link{symbolmap}} to control the legend.
    In addition to arguments documented under
    \code{\link{plot.symbolmap}}, and graphical arguments recognised
    by \code{\link{symbolmap}}, the list may also include the argument
    \code{sep} giving the separation between the main plot and the
    legend, or \code{sep.frac} giving the separation as a fraction
    of the relevant dimension (width or height) of the main plot.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title and the window containing the network.
  }
  \item{show.window}{
    Logical value indicating whether to plot the
    window containing the network. Overrides \code{show.all}.
  }
  \item{show.network}{
    Logical value indicating whether to plot the network.
  }
  \item{do.plot}{
    Logical value determining whether to actually perform the plotting.
  }
  \item{multiplot}{
    Logical value giving permission to display multiple plots.
  }
}
\details{
  The linear network is plotted by \code{\link{plot.linnet}},
  then the points are plotted using code equivalent to \code{\link{plot.ppp}}.

  Commonly-used arguments include:
  \itemize{
    \item \code{col} and \code{lwd} for the colour and width of lines
    in the linear network
    \item \code{cols} for the colour or colours of the points
    \item \code{chars} for the plot characters representing different
    types of points
    \item \code{shape} to control the shape of the symbol
    (this argument takes precedence over \code{chars}).
  }
  These are documented in the help file for
  \code{\link{plot.ppp}}.

  If \code{shape="crossticks"}, the points are
  drawn as short line segments perpendicular to the network.

  Note that the linear network will be plotted
  even when \code{add=TRUE}, unless \code{show.network=FALSE}.
}
\value{
  (Invisible) object of class \code{"symbolmap"}
  giving the correspondence between 
  mark values and plotting characters.
}
\seealso{
  \code{\link{lpp}}.

  See \code{\link{plot.ppp}} for options for representing the points.

  See also \code{\link{points.lpp}}, \code{\link{text.lpp}}.
}
\examples{
  plot(chicago, cols=1:7)
  plot(dendrite, shape="crossticks", cols=2:4, size=8,
                 leg.side="bottom", leg.args=list(lwd=3))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
