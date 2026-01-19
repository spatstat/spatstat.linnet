\name{lppm}
\alias{lppm}
\alias{lppm.formula}
\alias{lppm.lpp}
\title{
  Fit Point Process Model to Point Pattern on Linear Network
}
\description{
  Fit a point process model to a point pattern dataset on a linear network
}
\usage{
lppm(X, ...)

\method{lppm}{formula}(X, interaction=NULL, ..., data=NULL)

\method{lppm}{lpp}(X, ..., eps=NULL, nd=1000, random=FALSE)
}
\arguments{
  \item{X}{
    Either an object of class \code{"lpp"} specifying a point pattern
    on a linear network, or a \code{formula} specifying the
    point process model.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.model]{ppm}}.
  }
  \item{interaction}{
    An object of class \code{"interact"}
    describing the point process interaction
    structure, or \code{NULL} indicating that a Poisson process (stationary
    or nonstationary) should be fitted.
  }
  \item{data}{
    Optional. The values of spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A list whose entries are images,
    functions, windows, tessellations or single numbers. 
  }
  \item{eps}{
    Optional. Spacing between dummy points along each segment of the
    network. 
  }
  \item{nd}{
    Optional. Total number of dummy points placed on 
    the network. Ignored if \code{eps} is given.
  }
  \item{random}{
    Logical value indicating whether the grid of dummy points should be
    placed at a randomised starting position.
  }
}
\details{
  This function fits a point process model to data that specify
  a point pattern on a linear network. It is a counterpart of
  the model-fitting function \code{\link[spatstat.model]{ppm}} designed
  to work with objects of class \code{"lpp"} instead of \code{"ppp"}.

  The function \code{lppm} is generic, with methods for
  the classes \code{formula} and \code{lppp}.

  In \code{lppm.lpp}
  the first argument \code{X} should be an object of class \code{"lpp"}
  (created by the command \code{\link{lpp}}) specifying a point pattern
  on a linear network.

  In \code{lppm.formula},
  the first argument is a \code{formula} in the \R language
  describing the spatial trend model to be fitted. It has the general form
  \code{pattern ~ trend} where the left hand side \code{pattern} is usually
  the name of a point pattern on a linear network
  (object of class \code{"lpp"})
  to which the model should be fitted, or an expression which evaluates
  to such a point pattern;
  and the right hand side \code{trend} is an expression specifying the
  spatial trend of the model.
  Variable names which appear in the \code{trend} can be
  \itemize{
    \item the name of an object in the current environment
    \item the name of an entry in the list \code{covariates}
    \item one of the reserved names \code{x}, \code{y}, \code{seg},
    \code{tp} representing respectively the spatial coordinates \eqn{x}
    and \eqn{y}, and the local coordinates \code{seg} (line segment
    index) and \code{tp} (relative position along the segment).
  }
  Covariates which are objects in the environment or entries in the
  list \code{covariates} may have any of the following formats:
  \describe{
    \item{a pixel image,}{
      giving the values of a spatial covariate at 
      a fine grid of locations. It should be an object of
      class \code{"im"}, see \code{\link[spatstat.geom]{im.object}},
      or class \code{"linim"}, see \code{\link[spatstat.linnet]{linim}}.
    }
    \item{a function,}{
      which can be evaluated
      at any location on the network to obtain the value of the spatial
      covariate.
      This may be a function of class \code{"linfun"} (function on a
      network) or \code{"funxy"} (function in two dimensional space).
      Alternatively it may be any function in the \R language:
      the first two arguments of the function should be the
      Cartesian coordinates \eqn{x} and \eqn{y}. The function may have
      additional arguments include \code{seg}, \code{tp} and
      \code{marks} and other arguments.
    }
    \item{a window,}{
      interpreted as a logical variable
      which is \code{TRUE} inside the window and \code{FALSE} outside
      it. This should be an object of class \code{"owin"}.
    }
    \item{a tessellation,}{
      interpreted as a factor covariate.
      For each spatial location, the factor value indicates
      which tile of the tessellation it belongs to.
      This should be an object of class \code{"tess"} or \code{"lintess"}.
    }
    \item{a single number,}{indicating a covariate that is
      constant in this dataset.
    }
  }

  Other arguments \code{...} are passed from \code{lppm.formula}
  to \code{lppm.lpp} and from \code{lppm.lpp} to \code{\link[spatstat.model]{ppm}}.
}
\value{
  An object of class \code{"lppm"} representing the fitted model.
  There are methods for \code{print}, \code{predict},
  \code{coef} and similar functions.
}
\author{
  \adrian
  and Greg McSwiggan.
}
\seealso{
  \code{\link{methods.lppm}},
  \code{\link{predict.lppm}},
  \code{\link[spatstat.model]{ppm}},
  \code{\link{lpp}}.
}
\examples{
  X <- runiflpp(15, simplenet)
  lppm(X ~1)
  lppm(X ~x)
  marks(X) <- factor(rep(letters[1:3], 5))
  lppm(X ~ marks)
  lppm(X ~ marks * x)
}
\references{
  McSwiggan, G. (2019)  
  Spatial point process methods for linear networks
  with applications to road accident analysis.
  PhD thesis, University of Western Australia.
}
\keyword{spatial}
\keyword{models}
\concept{Linear network}
