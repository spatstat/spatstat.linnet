#'  eval.linim.R
#'
#'  Analogues of eval.im, harmonise.im and im.apply for 'linim'
#'
#'  $Revision: 1.5 $ $Date: 2025/11/28 06:19:31 $

harmonize.linim <- harmonise.linim <- function(...) {
  argz <- list(...)
  n <- length(argz)
  if(n < 2) return(argz)
  result <- vector(mode="list", length=n)
  isim <- sapply(argz, is.im)
  isfun <- sapply(argz, is.function)
  islinim <- sapply(argz, is.linim)
  islinnet <- sapply(argz, is.linnet)
  islinfun <- sapply(argz, inherits, what="linfun")
  hasnet <- islinim | islinnet | islinfun
  if(!any(hasnet))
    stop("No network information was supplied", call.=FALSE)
  ## Extract the network
  nets <- unique(solapply(argz[hasnet], as.linnet))
  if(length(nets) > 1) 
    stop("Not all data refer to the same linear network", call.=FALSE)
  L <- nets[[1L]]
  ## Select argument which determines the common resolution
  imax <- NULL
  if(any(islinim)) {
    ## Use linim objects to determine resolution
    if(sum(islinim) == 1) {
      imax <- which(islinim) # the unique linim object
    } else {
      ## Extract data frames of linim objects
      dframes <- lapply(argz[islinim], attr, which="df")
      dfexists <- !sapply(dframes, is.null)
      if(any(dfexists)) {
        if(sum(dfexists) == 1) {
          ## Use the unique linim object which has a data frame
          imax <- which(islinim)[which(dfexists)]
        } else {
          ## find the data frame with finest spacing
          nr <- sapply(dframes[dfexists], nrow)
          imax <- which(islinim)[which(dfexists)[which.max(nr)]]
        }
      } # else there are no data frames to compare
    }
  }
  if(is.null(imax) && any(isim)) {
    ## Find image object with finest pixel resolution
    xstep <- sapply(argz[[isim]], attr, which="xstep")
    ystep <- sapply(argz[[isim]], attr, which="ystep")
    imax <- which.min(xstep * ystep)
  }
  if(is.null(imax)) {
    ## There is no resolution information.
    ## Convert everything at the same default resolution.
    argz <- solapply(argz, as.linim, L=L)
  } else {
    ## Use argument number 'imax' to determine resolution
    template <- argz[[imax]]
    G <- as.im(template)
    dftemplate <- if(is.linim(template)) attr(template, which="df") else NULL
    if(is.null(dftemplate)) {
      ## Convert everything using the pixel grid of 'template'
      argz <- solapply(argz, as.linim, L=L, xy=G)
    } else {
      ## Convert everything using the grid and data frame of template
      argz <- solapply(argz, as.linim, L=template, xy=G)
    }
  }
  return(argz)
}

#' analogue of eval.im

eval.linim <- function(expr, envir, harmonize=TRUE, warn=TRUE) {
  sc <- sys.call()
  # Get names of all variables in the expression
  e <- as.expression(substitute(expr))
  varnames <- all.vars(e)
  allnames <- all.names(e, unique=TRUE)
  funnames <- allnames[!(allnames %in% varnames)]
  if(length(varnames) == 0)
    stop("No variables in this expression")
  # get the values of the variables
  if(missing(envir)) {
    envir <- parent.frame() # WAS: sys.parent()
  } else if(is.list(envir)) {
    envir <- list2env(envir, parent=parent.frame())
  }
  vars <- mget(varnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
  funs <- mget(funnames, envir=envir, inherits=TRUE, ifnotfound=list(NULL))
  # Find out which variables are (linear) images
  islinim <- unlist(lapply(vars, inherits, what="linim"))
  if(!any(islinim))
    stop("There are no linear images (class linim) in this expression")
  # ....................................
  # Evaluate the pixel values using eval.im
  # ....................................
  sc[[1L]] <- as.name('eval.im')
  sc$envir <- envir
  Y <- eval(sc)
  # .........................................
  # Then evaluate data frame entries if feasible
  # .........................................
  dfY <- NULL
  linims <- vars[islinim]
  nlinims <- length(linims)
  dframes <- lapply(linims, attr, which="df")
  nets <- lapply(linims, attr, which="L")
  isim <- unlist(lapply(vars, is.im))
  if(!any(isim & !islinim)) {
    # all images are 'linim' objects
    # Check that the images refer to the same linear network
    if(nlinims > 1) {
      agree <- unlist(lapply(nets[-1L], identical, y=nets[[1L]]))
      if(!all(agree))
        stop(paste("Images do not refer to the same linear network"))
    }
    dfempty <- unlist(lapply(dframes, is.null))
    if(!any(dfempty)) {
      # ensure data frames are compatible
      if(length(dframes) > 1 && (
          length(unique(nr <- sapply(dframes, nrow))) > 1   ||
           !allElementsIdentical(dframes, "seg")            ||
   	   !allElementsIdentical(dframes, "tp")
	)) {
        # find the one with finest spacing
	imax <- which.max(nr)
	# resample the others
	dframes[-imax] <- lapply(dframes[-imax],
	                         resampleNetworkDataFrame,
	                         template=dframes[[imax]])
      }
      # replace each image variable by its data frame column of values
      vars[islinim] <- lapply(dframes, getElement, "values")
      # now evaluate expression
      Yvalues <- eval(e, append(vars, funs))
      # pack up
      dfY <- dframes[[1L]]
      dfY$values <- Yvalues
    }
  }
  result <- linim(nets[[1L]], Y, df=dfY, restrict=FALSE)
  return(result)
}

resampleNetworkDataFrame <- function(df, template) {
  # resample 'df' at the points of 'template'
  invalues  <- df$values
  insegment <- df$mapXY
  inteepee  <- df$tp
  out <- template
  n <- nrow(out)
  outvalues <- vector(mode = typeof(invalues), length=n)
  outsegment <- out$mapXY
  outteepee  <- out$tp
  for(i in seq_len(n)) {
    relevant <- which(insegment == outsegment[i])
    if(length(relevant) > 0) {
      j <- which.min(abs(inteepee[relevant] - outteepee[i]))
      outvalues[i] <- invalues[relevant[j]]
    }
  }
  out$values <- outvalues
  return(out)
}

linim.apply <- function(X, FUN, ...,
                     fun.handles.na=FALSE, check=TRUE, verbose=TRUE) {
  ## determine function to be applied
  fun <- if(is.character(FUN)) get(FUN, mode="function") else
         if(is.function(FUN)) FUN else stop("Unrecognised format for FUN")
  funcode <- match(list(fun),
                   list(base::sum,
                        base::mean,
                        base::mean.default,
                        stats::var,
                        stats::sd),
                   nomatch=0L)
  funtype <- c("general", "sum", "mean", "mean", "var", "sd")[funcode+1L]
  ## ensure images are compatible
  if(check) 
    X <- do.call(harmonise.linim, X)
  if(!all(sapply(X, is.linim)))
    stop("All elements of X must be pixel images on a network (class linim)",
         call.=FALSE)
  ## extract linear network
  L <- as.linim(X[[1]])
  ## First apply function to pixel values
  result <- im.apply(X, FUN, ..., fun.handles.na=fun.handles.na,
                     check=check, verbose=verbose)
  single <- is.im(result)
  ## Consider data for sample points on network
  dframes <- lapply(X, attr, which="df") 
  if(any(sapply(dframes, is.null))) {
    ## Give up and convert result to linim
    if(single) {
      result <- linim(L=L, Z=result)
    } else {
      result <- solapply(result, linim, L=L)
    }
  } else {
    ## Extract values at sample points
    vals <- sapply(dframes, getElement, name="values")
    ## apply function to these values
    y <- ImApplyEngine(vals, fun, funtype, fun.handles.na, ...)
    ## attach to results
    df <- dframes[[1L]]
    if(single) {
      df$values <- y
      result <- linim(L=L, Z=result, df=df)
    } else {
      n <- length(result)
      class(result) <- "list"
      for(i in seq_len(n)) {
        df$values <- y[,i]
        result[[i]] <- linim(result[[i]], L=L, df=df)
      }
      class(result) <- "solist"
    }
  }
  return(result)
}

