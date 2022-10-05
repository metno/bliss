#+ Define correlations 
corr1d <- function( values, par, label) {
# correlations when input vectors are referred to a specific point (e.g. distances from a point)
#-----------------------------------------------------------------------------
  if ( length( par == 1)) {
    dist <- values
    dh   <- par
    if (label == "gaussian") {
      res <- exp( -0.5* (dist*dist) / (dh*dh) )
    } else if (label == "soar")  {
      res <- (1+dist/dh)*exp(-dist/dh)
    } else if (label == "powerlaw")  {
    res <- 1 / (1 + 0.5*(dist*dist)/(dh*dh))
    } else if (label == "toar")  {
      res <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
    } else {
      res <- rep( NA, length(dist))
    }
  }
  res
}
