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

