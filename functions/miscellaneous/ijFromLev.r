#+ position of coefficients of a given spatial resolution level into the vector of all wavelet coefficients
ijFromLev <- function( n, lev, fw=F) { 
# n = max number of levels (i.e. 2**n is the length of one dimension of the dyadic domain)
# lev = spatial resolution level (1=finer; ... to coarser)
# fw = are these the father wavelet coefficients?
#-------------------------------------------------------------------
  i <- c( NA, NA, NA)
  j <- c( NA, NA, NA)
  if ( lev == 1) { 
    i[1] <- 1
  } else {
    i[1] <- sum( 3 * ( 2**n / 2**(1:(lev-1)))**2) + 1
  }
  i[2] <- i[1] + ( 2**n / 2**lev)**2
  i[3] <- i[2] + ( 2**n / 2**lev)**2
  j[3] <- i[1] + 3 * ( 2**n / 2**lev)**2 - 1
  j[2] <- j[3] - ( 2**n / 2**lev)**2
  j[1] <- j[2] - ( 2**n / 2**lev)**2 
  if (fw) {
    i <- j[3] + 1
    j <- i + ( 2**n / 2**lev)**2 - 1
    return( c(i,j))
  }
  return( cbind(i,j))
}
  

