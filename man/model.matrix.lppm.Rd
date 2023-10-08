\name{model.matrix.lppm}
\alias{model.matrix.lppm}
\title{Extract Design Matrix from Point Process Model on a Network}
\description{
  Given a point process model that has been
  fitted to spatial point pattern data on a linear network,
  this function extracts the design matrix of the model.
}
\usage{
   \method{model.matrix}{lppm}(object,
                              data=model.frame(object, na.action=NULL),
                              \dots, 
                              keepNA=TRUE)
}
\arguments{
  \item{object}{
    The fitted point process model. An object of class \code{"lppm"}.
  }
  \item{data}{
    A model frame, containing the data required for the Berman-Turner device.
  }
  \item{keepNA}{
    Logical. Determines whether rows containing NA values will be
    deleted or retained.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link[stats:model.matrix]{model.matrix.lm}}.
  }
}
\details{
  This is a method for the generic function
  \code{\link[stats]{model.matrix}}.
  It extracts the design matrix of a spatial point process model
  on a linear network (object of class \code{"lppm"}).

  More precisely, this command extracts
  the design matrix of the generalised linear model associated with
  a spatial point process model. 
  
  The \code{object} must be a fitted point process model
  on a network (object of class \code{"lppm"})
  produced by the model-fitting function \code{\link{lppm}}.
  The method \code{model.matrix.lppm}
  extracts the model matrix for the GLM.

  The result is a matrix, with one row for every quadrature point
  in the fitting procedure, and one column for every canonical
  covariate in the design matrix.

  If there are \code{NA} values in the covariates,
  the argument \code{keepNA} determines whether to retain or delete
  the corresponding rows of the model matrix. The default
  \code{keepNA=TRUE} is to retain them. Note that this differs from
  the default behaviour of many other methods for \code{model.matrix},
  which typically delete rows containing \code{NA}.
}
\value{
  A matrix. Columns of the matrix are canonical covariates in the model.
  Rows of the matrix correspond to quadrature points
  in the fitting procedure (provided \code{keepNA=TRUE}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[stats]{model.matrix}},
  \code{\link{model.images.lppm}},
  \code{\link{lppm}}
}

\examples{
   fit <- lppm(spiders ~ x + y)
   head(model.matrix(fit))
}
\keyword{spatial}
\keyword{models}
\concept{Linear network}
