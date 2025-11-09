\name{roc.lpp}
\alias{roc.lpp}
\alias{roc.lppm}
\title{
  Receiver Operating Characteristic for Data on a Network
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a point pattern on a linear network
  or a fitted point process model on a linear network.
}
\usage{
\method{roc}{lpp}(X, covariate, 
                  \dots,
                  baseline = NULL, high = TRUE, weights = NULL,
                  method = "raw",
                  CI = "none", alpha=0.05,
                  subset=NULL)

\method{roc}{lppm}(X, covariate=NULL,
                   \dots, baseline=NULL, high=TRUE,
                   method = "raw",
                   CI = "none", alpha=0.05,
                   leaveoneout=FALSE, subset=NULL)
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
    Traditionally omitted when \code{X} is a fitted model.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} controlling the
    pixel resolution for calculations.
  }
  \item{baseline}{
    Optional. A spatial object giving a baseline intensity.
    Usually a \code{function(x,y)} or
    a pixel image (object of class \code{"im"} or \code{"linim"})
    giving the baseline intensity at any location on the network.
    Alternatively a point pattern on the network
    (object of class \code{"lpp"})
    giving the locations of the reference population.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
  \item{weights}{
    Optional. Numeric vector of weights attached to the data points.
  }
  \item{method}{
    The method or methods that should be used to estimate the ROC curve.
    A character vector: current choices are
    \code{"raw"}, \code{"monotonic"}, \code{"smooth"} and \code{"all"}.
    See Details.
  }
  \item{CI}{
    Character string (partially matched) specifying whether confidence
    intervals should be computed, and for which method.
    See Details.
  }
  \item{alpha}{
    Numeric value between 0 and 1. The confidence intervals will have
    confidence level \code{1-alpha}. The default gives 95\%
    confidence intervals.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, from which the ROC should be
    calculated.
  }
  \item{leaveoneout}{
    Logical value specifying (for \code{roc.lppm}) whether the fitted intensity
    of the model at each of the original data points should be computed
    by the leave-one-out procedure 
    (i.e. by removing the data point in question from the point pattern,
    re-fitting the model to the reduced point pattern, and computing the
    intensity of this modified model at the point in question)
    as described in \rocpapercite.
    It is also possible to specify \code{leaveoneout=c(TRUE,FALSE)}
    so that both versions are calculated.
  }

}
\details{
  The command \code{roc}
  computes the Receiver Operating
  Characteristic curve.
  The area under the ROC is computed by \code{\link[spatstat.explore]{auc}}.

  The function \code{\link[spatstat.explore]{roc}} is generic,
  with methods for point patterns, fitted point process models
  and other kinds of data.
  
  This help file describes the methods for classes \code{"lpp"} and
  \code{"lppm"}.

  For a point pattern \code{X} and a covariate \code{Z}, the
  ROC is a plot showing the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  For each possible threshold \eqn{z}, the algorithm calculates
  the fraction \eqn{a(z)} of area in the study region where the
  covariate takes a value greater than \eqn{z}, and the
  fraction \eqn{b(z)} of data points for which the covariate value
  is greater than \eqn{z}. The ROC is a plot of \eqn{b(z)} against
  \eqn{a(z)} for all thresholds \eqn{z}. 
  
  For a fitted point process model, 
  the ROC shows the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  The ROC is \bold{not} a diagnostic for the goodness-of-fit of the model
  (Lobo et al, 2007).
}
\value{
  Function value table (object of class \code{"fv"})
  which can be plotted to show the ROC curve.
  Also belongs to class \code{"roc"}.
}
\references{
  \rocpaper.
  
  Lobo, J.M.,
  \Jimenez-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in 
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S. 
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \Birkhauser, Basel.
}
\author{
\spatstatAuthors.
}
\seealso{
\code{\link[spatstat.explore]{auc}}
}
\examples{
  plot(roc(spiders, "x"))
  fit <- lppm(spiders ~ x)
  plot(roc(fit))
}
\keyword{spatial}
\concept{Linear network}
