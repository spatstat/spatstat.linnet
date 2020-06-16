##  spatstat.linnet/R/First.R

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.linnet"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatLinnetVersion", vs)
  packageStartupMessage(paste("spatstat.linnet", vs))
  return(invisible(NULL))
}

  
