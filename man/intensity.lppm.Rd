\name{intensity.lppm}  
\alias{intensity.lppm}
\title{
  Intensity of Fitted Point Process Model on a Linear Network
}
\description{
  Computes the intensity of a fitted point process model
  on a linear network.
}
\usage{
 \method{intensity}{lppm}(X, \dots)
}
\arguments{
  \item{X}{
    A fitted point process model on a linear network
    (object of class \code{"lppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.lppm}}.
  }
}
\details{
  This is a method for the generic function \code{\link[spatstat.geom]{intensity}}
  for fitted point process models on a linear network (objects of
  class \code{"lppm"}) created by
  the model-fitting function \code{\link{lppm}}).

  The intensity of a point process model on a linear network is the expected
  number of random points per unit length (McSwiggan, 2019;
  Baddeley, Rubak and Turner, 2015, Chapter 17).

  The result of \code{intensity.lppm(X)} is a numerical value if \code{X}
  is a homogeneous Poisson point process, and a pixel image if \code{X}
  is inhomogeneous. (In the latter case, the resolution of the pixel
  image is controlled by the arguments \code{\dots} which are passed
  to \code{\link{predict.lppm}}.)
}
\value{
  A numeric value or numeric vector (if the model is homogeneous)
  or a pixel image (object of class \code{"linim"})
}
\references{
  \baddrubaturnbook.
  
  McSwiggan, G. (2019)  
  Spatial point process methods for linear networks
  with applications to road accident analysis.
  PhD thesis, University of Western Australia.

}
\seealso{
  \code{\link[spatstat.geom]{intensity}},
  \code{\link[spatstat.linnet]{intensity.lpp}},
  \code{\link[spatstat.model]{intensity.ppm}}.
  
  Type \code{methods(intensity)} to see methods for other classes.
}
\examples{
  m1 <- lppm(spiders ~ 1)
  intensity(m1)
  mx <- lppm(spiders ~ x)
  intensity(mx)
}
\author{
  \adrian and Greg McSwiggan.
}
\keyword{spatial}
\keyword{models}
