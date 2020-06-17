\name{model.images.lppm}  
\alias{model.images.lppm}
\title{Compute Images of Constructed Covariates}
\description{
  For a point process model fitted to spatial point pattern data
  on a linear network,
  this function computes pixel images of the covariates
  in the design matrix.
}
\usage{
  \method{model.images}{lppm}(object, L = as.linnet(object), ...)
}
\arguments{
  \item{object}{
    Fitted point process model on a linear network.
    An object of class \code{"lppm"}.
  }
  \item{L}{
    A linear network (object of class \code{"linnet"}) in which the
    images should be computed. Defaults to the network
    in which the model was fitted.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link[stats:model.matrix]{model.matrix.lm}}.
  }
}
\details{
  This command is similar to \code{\link{model.matrix.lppm}} except
  that it computes pixel images of the covariates,
  instead of computing the covariate values at certain points only.

  The \code{object} must be a fitted spatial point process model
  on a linear network (object of class \code{"lppm"}
  produced by the model-fitting
  function \code{\link{lppm}}).

  The spatial covariates required by the model-fitting procedure
  are computed at every location on the network \code{L}.

  Note that the spatial covariates computed here
  are not necessarily the original covariates
  that were supplied when fitting the
  model. Rather, they are the canonical covariates,
  the covariates that appear in the
  loglinear representation of the (conditional) intensity
  and in the columns of the design matrix. For example, they might include
  dummy or indicator variables for different levels of a factor,
  depending on the contrasts that are in force.

  The format of the result depends on whether the original point pattern
  data were marked or unmarked.
  \itemize{
    \item
    If the original dataset was unmarked,
    the result is a named list of pixel images on the network (objects of class
    \code{"linim"}) containing the values of the spatial covariates.
    The names of the list elements are the names of the covariates
    determined by \code{\link[stats:model.matrix]{model.matrix.lm}}.
    The result is also of class \code{"solist"} so that it can
    be plotted immediately.
    \item 
    If the original dataset was a multitype point pattern,
    the result is a \code{\link{hyperframe}}
    with one column for each possible type of points.
    Each column is a named list of pixel images on the network (objects of class
    \code{"linim"}) containing the values of the spatial covariates.
    The row names of the hyperframe are the names of the covariates
    determined by \code{\link[stats:model.matrix]{model.matrix.lm}}.
  }
  The pixel resolution is determined by the arguments \code{\dots}
  and \code{\link[spatstat.core]{spatstat.options}}.

}
\value{
  A list (of class \code{"solist"}) or
  array (of class \code{"hyperframe"}) containing
  pixel images on the network (objects of class \code{"linim"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.core]{model.matrix.ppm}},
  \code{\link[stats]{model.matrix}},
  \code{\link{lppm}}.
}
\examples{
   fit <- lppm(spiders ~ x + polynom(y, 2))
   model.images(fit)
}
\keyword{spatial}
\keyword{models}
