\name{auc.lpp}
\alias{auc.lpp}
\alias{auc.lppm}
\title{
  Area Under ROC Curve for Network Data
}
\description{
  Compute the AUC (area under the Receiver Operating Characteristic
  curve) for a point pattern on a network,
  or a fitted point process model on a network.
}
\usage{
\method{auc}{lpp}(X, covariate, \dots, high = TRUE,
                     subset=NULL)

\method{auc}{lppm}(X, \dots, subset=NULL)
}
\arguments{
  \item{X}{
    Point pattern on a network (object of class \code{"lpp"})
    or fitted point process model on a network
    (object of class \code{"lppm"}).
  }
  \item{covariate}{
    Spatial covariate. Either a \code{function(x,y)},
    a pixel image (object of class \code{"im"} or \code{"linim"}), or
    one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates.    
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.explore]{roc}},
    and arguments passed to \code{\link[spatstat.geom]{as.mask}}
    controlling the pixel resolution for calculations.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, for which the AUC should be
    calculated.
  }
}
\details{
  The generic \code{\link[spatstat.explore]{auc}} computes the AUC,
  the area under the curve of the Receiver Operating Characteristic.
  The ROC curve itself is computed by the generic \code{\link[spatstat.explore]{roc}}.

  The functions \code{auc.lpp} and \code{auc.lppm}
  are methods for \code{\link[spatstat.explore]{auc}} for
  point patterns on a linear network (class \code{"lpp"}) and
  fitted point process models on a linear network (class \code{"lppm"}).

  For a point pattern \code{X} and a covariate \code{Z}, the
  AUC is a numerical index that measures the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  Let \eqn{x_i}{x[i]} be a randomly-chosen data point from \code{X}
  and \eqn{U} a randomly-selected location in the study region.
  The AUC is the probability that
  \eqn{Z(x_i) > Z(U)}{Z(x[i]) > Z(U)}
  assuming \code{high=TRUE}.
  That is, AUC is the probability that a randomly-selected data point
  has a higher value of the covariate \code{Z} than does a
  randomly-selected spatial location. The AUC is a number between 0 and 1. 
  A value of 0.5 indicates a complete lack of discriminatory power.
  
  For a fitted point process model \code{X},
  the AUC measures the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  Suppose \eqn{\lambda(u)}{\lambda(u)} is the intensity function of the model.
  The AUC is the probability that
  \eqn{\lambda(x_i) > \lambda(U)}{\lambda(x[i]) > \lambda(U)}.
  That is, AUC is the probability that a randomly-selected data point
  has higher predicted intensity than does a randomly-selected spatial
  location.
  The AUC is \bold{not} a measure of the goodness-of-fit of the model
  (Lobo et al, 2007).
}
\value{
  Numeric.
  For \code{auc.lpp}, the result is a single number
  giving the AUC value.
  
  For \code{auc.lppm}, the result is a 
  numeric vector of length 2 giving the AUC value 
  and the theoretically expected AUC value for this model.
}
\references{
  \rocpaper.
  
  Lobo, J.M.,
  \ifelse{latex}{\out{Jim{\'e}nez}}{Jimenez}-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in 
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S. 
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \ifelse{latex}{\out{Birkh{\"a}user}}{Birkhauser}, Basel.
}
\author{
  \rocketAuthors.
}
\seealso{
  \code{\link[spatstat.explore]{roc}}.
  
  \code{\link[spatstat.explore]{auc}},
  \code{\link[spatstat.model]{auc.ppm}}.

  \code{\link[spatstat.explore]{youden}}.

}
\examples{
  Crimes <- unmark(chicago)
  fit <- lppm(Crimes ~ x)
  auc(fit)
  auc(Crimes, "x")
}
\keyword{spatial}
