#' shortestpath.R
#'
#'    Shortest path between two points on a network
#' 
#' Copyright (c) 2024 Adrian Baddeley 
#' GNU Public Licence (>= 2.0)
#'
#' $Revision: 1.4 $ $Date: 2024/11/08 00:54:27 $

shortestpath <- function(X, i=1, j=2) {
  stopifnot(is.lpp(X))
  check.1.integer(i)
  check.1.integer(j)
  nX <- npoints(X)
  if(nX < 2 || min(i,j) < 1 || max(i,j) > nX) {
    stopifnot(npoints(X) >= 2)
    stopifnot(i >= 1 && i <= npoints(X))
    stopifnot(j >= 1 && j <= npoints(X))
  }
  #' extract only the two endpoints; henceforth they are points 1 and 2
  X <- X[c(i,j)]
  #' network information
  W <- Window(X)
  L <- domain(X)
  from <- L$from
  to   <- L$to  
  len <- lengths_psp(L$lines)
  nedges <- length(len)
  XP <- as.ppp(unmark(X))
  #' coordinates of endpoints
  coX <- coords(X)
  START <- 1L
  END   <- 2L
  segSTART <- coX$seg[START]
  segEND <- coX$seg[END]
  tpSTART <- coX$tp[START]
  tpEND <- coX$tp[END]
  #' check whether start and finish are on the same segment
  if(segSTART == segEND) {
    Path <- as.psp(from=XP[START], to=XP[END])
    attr(Path, "steps") <- integer(0)
    return(Path)
  }
  #' distance from each vertex to end point
  dv <- nnfromvertex(X[END], "dist")
  #' vertices of segment containing end point
  vENDleft <- from[segEND]
  vENDright <- to[segEND]
  #' initialise
  steps <- integer(0)
  #' consider segment containing start point
  vleft <- from[segSTART]
  vright <- to[segSTART]
  dleft <- tpSTART * len[segSTART] + dv[vleft]
  dright <- (1-tpSTART) * len[segSTART] + dv[vright]
  steps[1] <- vcurrent <- if(dleft < dright) vleft else vright
  #'
  while(length(steps) < nedges + 2) {
    #' find edges adjacent to current vertex 
    edgesfromcurrent <- which(from == vcurrent)
    edgestocurrent <- which(to == vcurrent)
    #' vertices adjacent to current vertex
    verticesfromcurrent <- to[edgesfromcurrent]
    verticestocurrent <- from[edgestocurrent]
    vadjacent <- c(verticesfromcurrent, verticestocurrent)
    #' calculate distance via each possible route
    dfrom <- len[edgesfromcurrent] + dv[verticesfromcurrent]
    dto   <- len[edgestocurrent]   + dv[verticestocurrent]
    stepsize <- c(dfrom, dto)
    #' detect case where an edge contains the target point
    if(vcurrent == vENDleft && any(vadjacent == vENDright)) {
      i <- which(vadjacent == vENDright)
      stepsize[i] <- len[segEND] * tpEND
      vadjacent[i] <- 0 
    } else if(vcurrent == vENDright && any(vadjacent == vENDleft)) {
      i <- which(vadjacent == vENDleft)
      stepsize[i] <- len[segEND] * (1-tpEND)
      vadjacent[i] <- 0
    }
    #' find step giving shortest path
    ibest <- which.min(stepsize)
    vbest <- vadjacent[ibest]
    if(vbest == 0) break
    steps <- c(steps, vbest)
    vcurrent <- vbest
  }
  #' form result
  V <- vertices(L)
  PathPoints <- superimpose(XP[START], V[steps], XP[END], W=W)
  Path <- as.psp(from=PathPoints[-npoints(PathPoints)],
                 to=PathPoints[-1L],
                 window=W)
  attr(Path, "steps") <- steps
  return(Path)
}

  
  
