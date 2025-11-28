\name{linim.apply}
\alias{linim.apply}
\title{
  Apply Function Pixelwise to List of Images on a Network
}
\description{
  Returns a pixel image obtained by applying a given function
  to corresponding pixels in a list of several pixel images
  on a linear network.
}
\usage{
linim.apply(X, FUN, \dots, fun.handles.na=FALSE, check=TRUE, verbose=TRUE)
}
\arguments{
  \item{X}{
    A list of pixel images on the same network
    (objects of class \code{"linim"}).
  }
  \item{FUN}{
    A function that can be applied to vectors,
    or a character string giving the name of such a function.
  }
  \item{\dots}{
    Additional arguments to \code{FUN}.
  }
  \item{fun.handles.na}{
    Logical value specifying what to do when the data include
    \code{NA} values. See Details.
  }
  \item{check}{
    Logical value specifying whether to check that the images in
    \code{X} are compatible (for example that they have the same
    grid of pixel locations) and to convert them to compatible images
    if necessary.
  }
  \item{verbose}{
    Logical value specifying whether to print informative messages.
  }
}
\details{
  The argument \code{X} should be a list of pixel images on a network
  (objects of class \code{"linim"}). They should all be defined
  on the same network. If the images do not have
  identical pixel grids, they will be converted to a common
  grid using \code{\link{harmonise.linim}}.

  At each pixel location, the values of the images in \code{X}
  at that pixel will be extracted as a vector; the function
  \code{FUN} will be applied to this vector; and the return value of
  \code{FUN} will become the pixel value of the resulting image.
  For example \code{linim.apply(X, mean)} will return a pixel image
  in which the value of each pixel is the average of the corresponding
  pixel values in the images in \code{X}.

  If the result of \code{FUN} is a vector, then the result of
  \code{linim.apply} will be a list of images.
  For example \code{linim.apply(X, range)} will
  return a list of two images containing the pixelwise minimum
  and pixelwise maximum, respectively, of the input images in \code{X}.

  The argument \code{fun.handles.na} specifies what to do
  when some of the pixel values are \code{NA}.
  \itemize{
    \item If \code{fun.handles.na=FALSE} (the default),
    the function \code{FUN} is never applied to data that include
    \code{NA} values; the result is defined to be \code{NA} whenever
    the data contain \code{NA}.
    \item
    If \code{fun.handles.na=TRUE}, the function \code{FUN} will be applied to
    all pixel data, including those which contain \code{NA} values.
  }
}
\value{
  A pixel image on a network (object of class \code{"linim"})
  or a list of pixel images on the same network.
}
\seealso{
  \code{\link{eval.linim}} for algebraic operations with images.
}
\examples{
  # list of two pixel images
  X <- runiflpp(5, simplenet, nsim=2)
  Y <- solapply(lapply(X, distfun), as.linim)
  plot(Y)
  linim.apply(Y, max)
  linim.apply(Y, sum)
  linim.apply(Y, range)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
