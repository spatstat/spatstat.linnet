\name{predict.lppm}
\alias{predict.lppm}
\title{
  Predict Point Process Model on Linear Network
}
\description{
  Given a fitted point process model on a linear network,
  compute the fitted intensity or conditional intensity of the model.
}
\usage{
\method{predict}{lppm}(object, ..., type = "trend",
                       locations = NULL, covariates = NULL,
                       se = FALSE,
                       new.coef=NULL)
}
\arguments{
  \item{object}{
    The fitted model. An object of class \code{"lppm"},
    see \code{\link{lppm}}.
  }
  \item{type}{
    Type of values to be computed. Either \code{"trend"} or
    \code{"cif"}. Currently ignored.
  }
  \item{locations}{
    Optional. Locations at which predictions should be computed.
    Either a data frame with two columns of coordinates,
    or a binary image mask.
  }
  \item{covariates}{
    Values of external covariates required by the model.
    Either a data frame, or a list of images and/or functions.
  }
  \item{se}{
    Logical value indicating whether to calculate standard errors as well.
  }
  \item{new.coef}{
    Optional. Numeric vector of model coefficients,
    to be used instead of the fitted coefficients
    \code{coef(object)} when calculating
    the prediction.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    to determine the
    pixel resolution (if \code{locations} is missing).
  }
}
\details{
  This function computes the fitted point process intensity,
  optionally with standard errors,
  for a point process model on a linear network. 
  It is a method for the generic \code{\link[stats]{predict}}
  for the class \code{"lppm"}.

  The argument \code{object} should be an object of class \code{"lppm"}
  (produced by \code{\link{lppm}}) representing a point process model
  on a linear network.

  Currently the argument \code{type} has no effect.
  The fitted intensity is computed in all cases.
  This occurs because currently all fitted models of class \code{"lppm"}
  are Poisson point processes, where the trend, intensity, and
  conditional intensity are the same.
  
  Predicted values are computed at the locations given by the
  argument \code{locations}. If this argument is missing,
  then predicted values are computed at a fine grid of points
  on the linear network.

  \itemize{
    \item 
    If \code{locations} is missing or \code{NULL} (the default),
    the return value is a pixel image (object of class \code{"linim"}
    and \code{"im"}) corresponding to a discretisation
    of the linear network, with numeric pixel values giving the
    predicted values at each location on the linear network.
    (If the model is multitype, the result is a list of such pixel
    images, one for each possible type of point.)
    \item
    If \code{locations} is a data frame, the result is a 
    numeric vector of predicted values at the locations specified by
    the data frame.
    \item
    If \code{locations} is a binary mask, the result is a pixel image
    with predicted values computed at the pixels of the mask.
    (If the model is multitype, the result is a list of such pixel
    images, one for each possible type of point.)
  }

  If \code{se=TRUE}, standard errors are also computed.
  The result is a list of two elements, each following the format
  described above; the first element contains the fitted estimates,
  and the second element contains the standard errors.
}
\value{
  If \code{se=FALSE} (the default), the result is 
  a pixel image (object of class \code{"linim"} and \code{"im"}) or
  a list of pixel images,
  or a numeric vector, depending on the argument \code{locations}.
  See Details.

  If \code{se=TRUE}, the result is a list of two elements,
  each with the format described above.
}
\author{
  \adrian
}
\seealso{
  \code{\link{lpp}},
  \code{\link{linim}}
}
\examples{
  X <- runiflpp(12, simplenet)
  fit <- lppm(X ~ x)
  v <- predict(fit, type="trend")
  plot(v)
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  Rakshit, S., McSwiggan, G., Nair, G. and Baddeley, A. (2021)
  Variable selection using penalised likelihoods for
  point patterns on a linear network.
 \emph{Australian and New Zealand Journal of Statistics}
 \bold{63}. DOI 10.1111/anzs.12341.

  Baddeley, A., Nair, G., Rakshit, S., McSwiggan, G. and Davies, T.M. (2021)
  Analysing point patterns on networks --- a review.
  \emph{Spatial Statistics} \bold{42}, 100435.
}
\keyword{spatial}
\keyword{models}
