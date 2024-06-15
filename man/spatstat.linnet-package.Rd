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
    \item \pkg{spatstat.explore} containing the main functionality
    for exploratory and non-parametric analysis of spatial data
    \item \pkg{spatstat.model} containing the main functionality
    for statistical modelling and inference for spatial data
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
\section{Overview of \pkg{spatstat.linnet}}{
  A linear network is a subset of the two-dimensional plane
  composed of straight line segments. It could represent a road network, for
  example. Our code requires that, if two segments intersect each other,
  then the intersection is a single point, and the intersection point is
  treated as a vertex of the network.

  The \pkg{spatstat.linnet} package supports spatial data analysis on
  a linear network. The primary aim is to analyse spatial patterns of
  points on a network. The points could represent road accidents on a
  road network, for example.

  The \pkg{spatstat.linnet} package provides code for handling
  \itemize{
    \item \code{linear networks}
    \item \code{point patterns on a linear network}
    \item \code{pixel images on a linear network} (where the network is
  divided into small segments and a numerical value is assigned to each segment)
    \item \code{functions on a linear network} (i.e. functions that are
    defined at every location along the network)
    \item \code{tessellations of a linear network} (where the network is
    subdivided into disjoint subsets with different labels)
    \item \code{point process models on a linear network}
  }

  Here is a list of the main functionality
  provided in \pkg{spatstat.linnet}.

  \bold{Linear networks}
  
  An object of class \code{"linnet"} represents a linear network.
  Examples of such objects include the dataset
  \code{\link[spatstat.data]{simplenet}} provided in the package.
  
  Linear network objects can be created by the following functions:

  \tabular{ll}{
    \code{\link{linnet}} \tab create a linear network \cr
    \code{\link{as.linnet}} \tab convert other data to a network \cr
    \code{\link{delaunayNetwork}} \tab network of Delaunay triangulation \cr
    \code{\link{dirichletNetwork}} \tab network of Dirichlet edges \cr
  }

  Utilities for manipulating networks include:
  \tabular{ll}{
    \code{\link{[.linnet}}\tab  extract subset of linear network \cr
    \code{\link{clickjoin}} \tab interactively join vertices in network \cr         \code{\link{joinVertices}} \tab join existing vertices in a network \cr
    \code{\link{insertVertices}} \tab insert new vertices at positions
    along network \cr
    \code{\link{addVertices}} \tab add new vertices, extending a network \cr
    \code{\link{thinNetwork}} \tab remove vertices or lines from a network \cr
    \code{\link{repairNetwork}} \tab repair internal format \cr
    \code{\link{vertices.linnet}} \tab extract the vertices of network \cr
    \code{\link{terminalvertices}} \tab find terminal vertices of
    network \cr
    \code{\link{affine.linnet}} \tab apply affine transformation \cr
    \code{\link{shift.linnet}} \tab apply vector translation \cr
    \code{\link{rotate.linnet}} \tab apply rotation \cr
    \code{\link{rescale.linnet}} \tab rescale the unit of length \cr
    \code{\link{scalardilate.linnet}} \tab physically rescale the
    network \cr
    \code{\link{diameter.linnet}} \tab diameter of linear  network \cr
    \code{\link{is.connected.linnet}} \tab  determine whether network
    is connected \cr
    \code{\link{lineardisc}} \tab compute disc of given radius in
    network \cr
    \code{\link{marks.linnet}} \tab extract marks of a network \cr
    \code{\link{marks<-.linnet}} \tab assign marks to a network \cr
    \code{\link{plot.linnet}} \tab plot a network \cr
    \code{\link{as.owin.linnet}} \tab extract window containing network \cr
    \code{\link{as.psp.linnet}} \tab extract line segments comprising
    network \cr
    \code{\link{nsegments.linnet}} \tab number of segments in network\cr
    \code{\link{nvertices.linnet}} \tab number of vertices in network\cr
    \code{\link{pixellate.linnet}} \tab convert network to 2D pixel
    image \cr
    \code{\link{print.linnet}} \tab print basic information \cr
    \code{\link{summary.linnet}} \tab print summary information \cr
    \code{\link{unitname.linnet}} \tab extract name of unit of length \cr
    \code{\link{unitname<-.linnet}} \tab assign name of unit of length \cr
    \code{\link{vertexdegree}} \tab number of segments meeting
    each vertex \cr
    \code{\link{volume.linnet}} \tab total length of network \cr
    \code{\link{Window.linnet}} \tab extract window containing network \cr
    \code{\link{density.linnet}} \tab smoothed 2D spatial density of lines \cr
  }

  A network is called a tree if it has no closed loops.
  The following functions support the creation and manipulation of
  trees:
  \tabular{ll}{
    \code{\link{begins}}\tab check start of character string\cr
    \code{\link{branchlabelfun}}\tab tree branch membership labelling
    function \cr
    \code{\link{deletebranch}}\tab delete a branch of a tree \cr
    \code{\link{extractbranch}}\tab extract a branch of a tree \cr
    \code{\link{treebranchlabels}}\tab label vertices of a tree by
    branch membership \cr
    \code{\link{treeprune}}\tab prune tree to given level\cr
  }

  \bold{Point patterns on a linear network}

    An object of class \code{"lpp"} represents a 
    point pattern on a linear network (for example,
    road accidents on a road network). 

    Examples of such objects include the following datasets
    provided in the \pkg{spatstat.data} package:
    
    \tabular{ll}{
      \code{\link[spatstat.data]{chicago}} \tab Chicago crime data \cr
      \code{\link[spatstat.data]{dendrite}} \tab Dendritic spines data \cr
      \code{\link[spatstat.data]{spiders}} \tab Spider webs on mortar lines of brick wall 
    }
    There is also a dataset provided in the extension package
    \pkg{spatstat.Knet}:
    \tabular{ll}{
      \code{wacrashes} \tab Road accidents in Western Australia
    }
    
    Point patterns on a network can be created by the following
    functions:

    \tabular{ll}{
      \code{\link{lpp}} \tab create a point pattern on a linear network \cr
      \code{\link{as.lpp}} \tab convert other data to point pattern on network \cr
      \code{\link{clicklpp}}\tab interactively add points on a linear
      Network \cr
      \code{\link{crossing.linnet}}\tab crossing points between network and other lines
    }

    Point patterns on a network can be generated randomly
    using the following functions:
    
    \tabular{ll}{
      \code{\link{rpoislpp}} \tab Poisson points on linear network \cr
      \code{\link{runiflpp}} \tab uniform random points on a linear network \cr
      \code{\link{rlpp}}\tab random points on a linear network\cr
      \code{\link{rSwitzerlpp}}\tab simulate Switzer-type point process on linear
      network \cr
      \code{\link{rThomaslpp}}\tab simulate Thomas process on linear network \cr
      \code{\link{rcelllpp}}\tab  simulate cell process on linear
      network \cr
      \code{\link{rjitter.lpp}}\tab randomly perturb a point pattern on
      a network \cr
    }

    Functions for manipulating a point pattern on a network include
    the following. An object of class \code{"lpp"} also belongs to the
    class \code{"ppx"}, for which additional support is available.
    
    \tabular{ll}{
      \code{\link{as.ppp.lpp}} \tab convert to 2D point pattern \cr
      \code{\link{as.psp.lpp}} \tab extract line segments \cr
      \code{\link[spatstat.geom]{marks.ppx}} \tab extract marks associated with points \cr
      \code{\link[spatstat.geom]{marks<-.ppx}} \tab assign marks to points on network \cr
      \code{\link{nsegments.lpp}} \tab count number of segments \cr
      \code{\link{print.lpp}} \tab print basic information \cr
      \code{\link{summary.lpp}} \tab print summary information \cr
      \code{\link{unitname.lpp}} \tab extract name of unit of length \cr
      \code{\link{unitname<-.lpp}} \tab assign name of unit of length \cr
      \code{\link{unmark.lpp}} \tab remove marks \cr
      \code{\link{subset.lpp}} \tab subset of points satisfying a
      condition \cr
      \code{\link{[.lpp}} \tab extract subset of point pattern\cr
      \code{\link{Window.lpp}} \tab extract window containing network \cr
      \code{\link{as.owin.lpp}} \tab extract window containing network \cr
      \code{\link{affine.lpp}} \tab apply affine transformation \cr
      \code{\link{shift.lpp}} \tab apply vector translation \cr
      \code{\link{rotate.lpp}} \tab apply rotation \cr
      \code{\link{rescale.lpp}} \tab rescale the unit of length \cr
      \code{\link{scalardilate.lpp}} \tab physically rescale the
      network and points \cr
      \code{\link{connected.lpp}}\tab find connected components of point
      pattern on network \cr
      \code{\link{cut.lpp}}\tab classify points in a Point Pattern on a
      Network \cr
      \code{\link{distfun.lpp}}\tab distance map (function) \cr
      \code{\link{distmap.lpp}}\tab distance map (image) \cr
      \code{\link{domain.lpp}}\tab  extract the linear network \cr
      \code{\link{identify.lpp}}\tab interactively identify points \cr
      \code{\link{is.multitype.lpp}}\tab recognize whether point pattern is
      multitype\cr
      \code{\link{nncross.lpp}}\tab nearest neighbours\cr
      \code{\link{nndist.lpp}}\tab  nearest neighbour distances \cr
      \code{\link{nnfromvertex}}\tab nearest data point from each vertex\cr
      \code{\link{nnfun.lpp}}\tab  nearest neighbour map \cr
      \code{\link{nnwhich.lpp}}\tab  identify nearest neighbours \cr
      \code{\link{pairdist.lpp}}\tab pairwise shortest-path distances \cr
      \code{\link{plot.lpp}}\tab plot point pattern on linear Network \cr
      \code{\link{points.lpp}}\tab draw points on existing plot \cr
      \code{\link{superimpose.lpp}}\tab  superimpose several point
      patterns \cr
      \code{\link{text.lpp}} \tab add text labels  \cr
      \code{\link{unstack.lpp}}\tab separate multiple columns of marks \cr
    }

    \bold{Pixel images on a network}

    An object of class \code{"linim"} represents a pixel image
    on a linear network. Effectively, the network is divided into small
    segments (lixels) and each small segment is assigned a value,
    which could be numeric, factor, logical or complex values.

    Pixel images on a network can be created using the following
    functions:

    \tabular{ll}{
      \code{\link{linim}}\tab create pixel image on linear network\cr
      \code{\link{as.linim}}\tab convert other data to pixel image on network\cr
    }

    Functions for manipulating a pixel image on a network include:

    \tabular{ll}{
      \code{\link{[.linim}} \tab extract subset of pixel image on linear network\cr
      \code{\link{[<-.linim}} \tab reset values in subset of image on linear network\cr
      \code{\link{Math.linim}}\tab S3 group generic methods for images on a linear network\cr
      \code{\link{eval.linim}}\tab evaluate expression involving pixel images on
      linear network\cr
      \code{\link{as.linnet.linim}} \tab extract linear network \cr
      \code{\link{integral.linim}}\tab integral of pixel image on a linear network\cr
      \code{\link{mean.linim}} \tab mean of pixel values \cr
      \code{\link{median.linim}} \tab median of pixel values \cr
      \code{\link{quantile.linim}} \tab quantiles of pixel values \cr
      \code{\link{as.data.frame.linim}} \tab convert to data frame \cr
      \code{\link{print.linim}} \tab print basic information \cr
      \code{\link{summary.linim}} \tab print summary information \cr
      \code{\link{affine.linim}} \tab apply affine transformation \cr
      \code{\link{scalardilate.linim}} \tab apply scalar dilation \cr
      \code{\link{shift.linim}} \tab apply vector translation \cr
      \code{\link{pairs.linim}} \tab scatterplot matrix for images \cr
      \code{\link{persp.linim}}\tab perspective view of pixel image on network\cr
      \code{\link{plot.linim}}\tab plot pixel image on linear network\cr
    }

    \bold{Functions on a linear network}

    An object of class \code{"linfun"} represents a function defined
    at any location along the network. Objects of this class are created
    by the following functions:

    \tabular{ll}{
      \code{\link{linfun}}\tab create function on a linear network \cr
      \code{\link{as.linfun}}\tab convert other data to function on network \cr
    }

    The following supporting code is available:

    \tabular{ll}{
      \code{\link{print.linfun}} \tab print basic information \cr
      \code{\link{summary.linfun}} \tab print summary information \cr
      \code{\link{plot.linfun}} \tab plot function on network \cr
      \code{\link{persp.linfun}}\tab perspective view of function on network\cr
      \code{\link{as.data.frame.linfun}} \tab convert to data frame \cr
      \code{\link{as.owin.linfun}} \tab extract window containing
      network \cr
      \code{\link{as.function.linfun}} \tab convert to ordinary \R
      function \cr
    }

    \bold{Tessellations of a linear network}

    An object of class \code{"lintess"} represents a tessellation of the
    network, that is, a subdivision of the network into disjoint subsets
    called \sQuote{tiles}. Objects of this class are created
    by the following functions:

    \tabular{ll}{
      \code{\link{lintess}}\tab create tessellation of network \cr
      \code{\link{chop.linnet}} \tab divide a linear network into tiles using
      infinite lines \cr
      \code{\link{divide.linnet}}\tab  divide linear network at cut
      points \cr
      \code{\link{lineardirichlet}}\tab Dirichlet tessellation on a linear network\cr
    }

    The following functions are provided for manipulating a tessellation
    on a network:

    \tabular{ll}{
      \code{\link{as.data.frame.lintess}}\tab convert to data frame \cr
      \code{\link{intersect.lintess}}\tab intersection of two
      tessellations on network \cr
      \code{\link{lineartileindex}}\tab determine which tile contains each
      given point on network \cr
      \code{\link{marks.lintess}}\tab extract marks of each tile \cr
      \code{\link{marks<-.lintess}}\tab assign marks to each tile \cr
      \code{\link{plot.lintess}}\tab plot tessellation on network \cr
      \code{\link{tile.lengths}}\tab compute lengths of tiles \cr
      \code{\link{tilenames.lintess}}\tab  names of tiles \cr
      \code{\link{as.linfun.lintess}}\tab convert tessellation to a
      function \cr
    }

  \bold{Smoothing a point pattern on a linear network:}

  Given a point pattern dataset on a linear network, it is often
  desired to estimate the spatially-varying density or intensity
  of points along the network. For example if the points represent
  road accidents, then we may wish to estimate the spatially-varying
  density of accidents per unit length (over a given period of time).

  Related tasks include estimation of relative risk, and smoothing of
  of values observed at the data points.

  \tabular{ll}{
    \code{\link{density.lpp}}\tab kernel estimate of intensity\cr
    \code{\link{densityEqualSplit}}\tab kernel estimate of intensity
    using equal-split algorithm \cr
    \code{\link{densityHeat.lpp}}\tab kernel estimate of intensity using heat equation\cr
    \code{\link{densityQuick.lpp}}\tab kernel estimate of intensity using a 2D kernel\cr
    \code{\link{densityVoronoi.lpp}}\tab intensity estimate using
    Voronoi-Dirichlet Tessellation\cr
    \code{\link{densityfun.lpp}}\tab kernel estimate of intensity as a
    function \cr
    \code{\link{bw.lppl}}\tab Bandwidth selection for kernel estimate of
    intensity \cr
    \code{\link{bw.voronoi}}\tab bandwidth selection for Voronoi estimator \cr
    \code{\link{relrisk.lpp}}\tab kernel estimate of relative risk\cr
    \code{\link{bw.relrisk.lpp}}\tab  Bandwidth selection for relative
    risk \cr
    \code{\link{Smooth.lpp}}\tab spatial smoothing of observations at
    points \cr

  }
  
  \bold{Exploration of dependence on a covariate:}

  Another task is to investigate how the spatially-varying intensity
  of points depends on an explanatory variable (covariate). The
  covariate may be given as a pixel image on the network
  (class \code{"linim"}) or
  as a function on the network (class \code{"linfun"}).

  \tabular{ll}{
    \code{\link{rhohat.lpp}}\tab nonparametric estimate of intensity as function
    of a covariate\cr
    \code{\link{roc.lpp}}\tab Receiver Operating Characteristic for data on a
    network\cr
    \code{\link{auc.lpp}}\tab  Area Under ROC Curve for data on a network\cr
    \code{\link{cdf.test.lpp}}\tab spatial distribution test for points on a
    linear network\cr
    \code{\link{berman.test.lpp}}\tab Berman's tests for point pattern
    on a network \cr
    \code{\link{sdr.lpp}}\tab Sufficient Dimension Reduction for a point
    pattern on a linear network\cr
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
    inhomogeneous pair correlation on linear network\cr
    \code{\link{linearJinhom}} \tab
    inhomogeneous \eqn{J} function on linear network \cr
    \code{\link{linearKEuclid}} \tab
    \eqn{K} function on linear network using Euclidean distance \cr
    \code{\link{linearKEuclidInhom}} \tab
    inhomogeneous \eqn{K} function on linear network using Euclidean distance\cr
    \code{\link{linearpcfEuclid}} \tab
    pair correlation function on linear network using Euclidean distance \cr
    \code{\link{linearpcfEuclidInhom}} \tab
    inhomogeneous pair correlation on linear network using Euclidean
    distance \cr
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
    \code{\link{data.lppm}} \tab extract original data \cr
    \code{\link{berman.test.lppm}} \tab Berman's tests of
    goodness-of-fit \cr
    \code{\link{is.marked.lppm}}\tab Recognise whether model is marked\cr
    \code{\link{is.multitype.lppm}}\tab Recognise whether model is multitype\cr
    \code{\link{is.stationary.lppm}}\tab Recognise whether model is
    stationary \cr
    \code{\link{model.frame.lppm}}\tab Extract the variables in model \cr
    \code{\link{model.images.lppm}}\tab Compute images of constructed covariates \cr
    \code{\link{model.matrix.lppm}}\tab Extract design matrix \cr
    \code{\link{plot.lppm}}\tab Plot fitted point process model\cr
    \code{\link{pseudoR2.lppm}}\tab Calculate Pseudo-R-Squared for model \cr
    \code{\link{simulate.lppm}}\tab simulate fitted point process model \cr
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
\concept{Linear network}
