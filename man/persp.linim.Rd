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
                       neg.args=list(), warncross=FALSE,
                       zadjust=1,
                       extrapolate=c("linear", "constant"))
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
    on the base plane (at height zero), to assist the perception of
    perspective. 
  }
  \item{ngrid}{
    Number of grid lines to draw, if \code{grid=TRUE}.
  }
  \item{col.grid}{
    Colour of grid lines, if \code{grid=TRUE}.
  }
  \item{col.base}{
    Colour of base plane.
    A single colour value, or a pixel image.
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
  \item{zadjust}{
    Adjustment factor for vertical scale, relative to the default scale.
  }
  \item{extrapolate}{
    Character string (partially matched) specifying how to extrapolate
    the value at the endpoint of each segment.
  }
}
\details{
  The pixel values are interpreted as the spatially-varying height of a
  vertical surface erected on each segment of the linear network. These
  surfaces are drawn in perspective view.
  This style of plot is often attributed to Okabe and Sugihara (2012).

  \enumerate{
    \item
    The horizontal plane at height zero is drawn, in perspective view,
    in the colour specified by \code{col.base}.
    
    If \code{col.base} is a pixel image, it will be rendered as a 
    colour image shown in perspective view on the
    horizontal plane. The argument \code{colmap} controls the mapping
    from pixel values of \code{col.base} to physical colours.

    \item
    A grid of lines on the horizontal plane is drawn
    if \code{grid=TRUE} (the default).
    \item 
    For each segment of the network, a vertical polygon 
    is constructed, with a straight lower edge aligned with the network segment,
    and a crooked upper edge whose height is proportional to 
    values of \code{x}. The polygon linearly interpolates between
    the values of pixels that lie along the segment.
    At each end of the segment,
    \itemize{
      \item If \code{extrapolate="linear"} (the default), the polygon
      height at the end of the segment is determined by linearly
      extrapolating from the two nearest pixel values.
      \item
      If \code{extrapolate="constant"}, the polygon height at the end of
      the segment is defined to be equal to the nearest pixel value.
    }
    The vertical polygons are drawn in the colour and style specified by
    the additional arguments \code{\dots}, for example, \code{col} for
    colour.
  }

  If \code{x} contains negative values, they will be represented as
  polygons extending downwards below the horizontal plane.
  These would be obscured if \code{col.base} is an opaque colour
  other than white, or if \code{col.base} is a pixel image.
  A transparent colour for \code{col.base} can be used
  if it is supported by the graphics device.

  Like all spatial plots in the \pkg{spatstat} family,
  \code{persp.linim} does not independently rescale
  the \eqn{x} and \eqn{y} coordinates. A long narrow window will be
  represented as a long narrow window in the perspective view.
  To override this and allow the coordinates to be independently
  rescaled, use the argument \code{scale=TRUE} which will be passed
  to \code{\link[graphics]{persp.default}}.
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
