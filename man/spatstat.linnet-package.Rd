\name{spatstat.linnet-package}
\alias{spatstat.linnet-package} 
\alias{spatstat.linnet} 
\docType{package}
\title{The spatstat.linnet Package}
\description{
  The \pkg{spatstat.linnet} package
  belongs to the \pkg{spatstat} family of packages.
  It contains the functionality
  for analysing spatial data on a linear network.
}
\details{
  \pkg{spatstat} is
  a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  The original \pkg{spatstat} package
  has now been split into several
  sub-packages.


  This sub-package \pkg{spatstat.linnet} contains the
  user-level functions from \pkg{spatstat} 
  that are concerned with spatial data on a linear network.
}
\section{Structure of the spatstat family}{  
  The orginal \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.core} containing the main functionality
    for statistical analysis and modelling of spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.
  For an overview of all the functions available in these sub-packages,
  see the help file for \pkg{spatstat} in the \pkg{spatstat} package,

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{Overview of functionality}{

  Here is a list of the main functionality in \pkg{spatstat.linnet}.
  
  \bold{Point patterns on a linear network}

    An object of class \code{"linnet"} represents a linear network
    (for example, a road network).

    \tabular{ll}{
      \code{\link{linnet}} \tab create a linear network \cr
      \code{\link{clickjoin}} \tab interactively join vertices in network \cr
      \code{spatstat.gui::iplot.linnet} \tab interactively plot network \cr
      \code{\link[spatstat.data]{simplenet}} \tab simple example of network \cr
      \code{\link{lineardisc}} \tab disc in a linear network \cr
      \code{\link{delaunayNetwork}} \tab network of Delaunay triangulation \cr
      \code{\link{dirichletNetwork}} \tab network of Dirichlet edges \cr
      \code{\link{methods.linnet}} \tab methods for \code{linnet} objects\cr
      \code{\link{vertices.linnet}} \tab nodes of network \cr
      \code{\link{joinVertices}} \tab join existing vertices in a network \cr
      \code{\link{insertVertices}} \tab insert new vertices at positions
      along a network \cr
      \code{\link{addVertices}} \tab add new vertices, extending a network \cr
      \code{\link{thinNetwork}} \tab remove vertices or lines from a network \cr
      \code{\link{repairNetwork}} \tab repair internal format \cr
      \code{\link{pixellate.linnet}} \tab approximate by pixel image
    }
    
    An object of class \code{"lpp"} represents a 
    point pattern on a linear network (for example,
    road accidents on a road network).
    
    \tabular{ll}{
      \code{\link{lpp}} \tab create a point pattern on a linear network \cr
      \code{\link{methods.lpp}} \tab methods for \code{lpp} objects \cr
      \code{\link{subset.lpp}} \tab method for \code{subset} \cr
      \code{\link{rpoislpp}} \tab simulate Poisson points on linear network \cr
      \code{\link{runiflpp}} \tab simulate random points on a linear network \cr
      \code{\link[spatstat.data]{chicago}} \tab Chicago crime data \cr
      \code{\link[spatstat.data]{dendrite}} \tab Dendritic spines data \cr
      \code{\link[spatstat.data]{spiders}} \tab Spider webs on mortar lines of brick wall 
    }
  \bold{Summary statistics for a point pattern on a linear network:}

  These are for point patterns on a linear network (class \code{lpp}).
  For unmarked patterns:
  
  \tabular{ll}{
    \code{\link{linearK}} \tab
    \eqn{K} function on linear network \cr
    \code{\link{linearKinhom}} \tab
    inhomogeneous \eqn{K} function on linear network \cr
    \code{\link{linearpcf}} \tab
    pair correlation function on linear network \cr
    \code{\link{linearpcfinhom}} \tab
    inhomogeneous pair correlation on linear network
  }

  For multitype patterns:
  \tabular{ll}{
    \code{\link{linearKcross}} \tab
    \eqn{K} function between two types of points \cr
    \code{\link{linearKdot}} \tab
    \eqn{K} function from one type to any type \cr
    \code{\link{linearKcross.inhom}} \tab
    Inhomogeneous version of \code{\link{linearKcross}} \cr
    \code{\link{linearKdot.inhom}} \tab
    Inhomogeneous version of \code{\link{linearKdot}} \cr
    \code{\link{linearmarkconnect}} \tab
    Mark connection function  on linear network \cr
    \code{\link{linearmarkequal}} \tab
    Mark equality function on linear network \cr
    \code{\link{linearpcfcross}} \tab
    Pair correlation between two types of points \cr
    \code{\link{linearpcfdot}} \tab
    Pair correlation from one type to any type \cr
    \code{\link{linearpcfcross.inhom}} \tab
    Inhomogeneous version of \code{\link{linearpcfcross}} \cr
    \code{\link{linearpcfdot.inhom}} \tab
    Inhomogeneous version of \code{\link{linearpcfdot}} 
  }

  Related facilities:
  
  \tabular{ll}{
    \code{\link{pairdist.lpp}} \tab distances between pairs  \cr
    \code{\link{crossdist.lpp}} \tab distances between pairs \cr
    \code{\link{nndist.lpp}} \tab nearest neighbour distances  \cr
    \code{\link{nncross.lpp}} \tab nearest neighbour distances  \cr
    \code{\link{nnwhich.lpp}} \tab find nearest neighbours  \cr
    \code{\link{nnfun.lpp}} \tab find nearest data point  \cr
    \code{\link{density.lpp}} \tab kernel smoothing estimator of intensity  \cr
    \code{\link{distfun.lpp}} \tab distance transform  \cr
    \code{\link{envelope.lpp}} \tab simulation envelopes  \cr
    \code{\link{rpoislpp}} \tab simulate Poisson points on linear network \cr
    \code{\link{runiflpp}} \tab simulate random points on a linear network 
  }
  
  It is also possible to fit point process models to \code{lpp} objects.

  \bold{Point process models on a linear network:}

  An object of class \code{"lpp"} represents a pattern of points on
  a linear network. Point process models can also be fitted to these
  objects. Currently only Poisson models can be fitted.

  \tabular{ll}{
    \code{\link{lppm}} \tab point process model on linear network \cr
    \code{\link{anova.lppm}} \tab analysis of deviance for \cr
    \tab point process model on linear network \cr
    \code{\link{envelope.lppm}} \tab simulation envelopes for \cr
    \tab point process model on linear network \cr
    \code{\link{fitted.lppm}} \tab fitted intensity values \cr
    \code{\link{predict.lppm}} \tab model prediction on linear network \cr
    \code{\link{linim}} \tab pixel image on linear network \cr
    \code{\link{plot.linim}} \tab plot a pixel image on linear network \cr
    \code{\link{eval.linim}} \tab evaluate expression involving images \cr
    \code{\link{linfun}} \tab function defined on linear network \cr
    \code{\link{methods.linfun}} \tab conversion facilities
  }
  
}

\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Ottmar Cronie,
  Tilman Davies,
  Greg McSwiggan and
  Suman Rakshit
  made substantial contributions of code.
}
\keyword{spatial}
\keyword{package}

