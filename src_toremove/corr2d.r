#+ Define correlations
corr2d <- function( values, par, label, dist2_is_global=F) {
# correlations between pair of points
#------------------------------------------------------------------------------
  if ( length( par == 1)) {
    dh <- par
    if (dist2_is_global) {
      dist2 <- envtmp$dist2
    } else {
      x  <- values[,1]
      y  <- values[,2]
      dist2 <- outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2
    }
    if (label == "gaussian") {
      res <- exp(-0.5 * dist2 / (dh*dh))
    } else if (label == "soar")  {
      distnorm <- sqrt(dist2) / dh 
      res <- (1+distnorm)*exp(-distnorm)
    } else if (label == "powerlaw")  {
      res <- 1 / (1 + 0.5*dist2 / (dh*dh))
    } else if (label == "toar")  {
      dist <- sqrt(dist2)
      res <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
    }
  }
  res
}
