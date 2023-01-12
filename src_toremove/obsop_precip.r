#+ observation operator, precipitation
obsop_precip <- function( r,
                          xcoord=NA,
                          ycoord=NA,
                          method="bilinear",
                          inf=0) {
# method. "simple" nearest neighbour. "bilinear" bilinear
# inf. lower limit for the plausible values
# - Global variables
# rmaster. raster with grid parameters
# n0. number of observations
# VecX, VecY. observation coordinates
# nens. number of ensemble members
# aix. indexes. the ensemble members (Xb) are defined only over aix cells
# Xb. matrix. ensemble members (aix cells, nens ensembles)
# - Returned values
# yb. matrix. ensemble at observation locations (ix observation loc, nens)
# ix. vector. indeces over the observation locations.
#------------------------------------------------------------------------------
  n0 <- length( xcoord)
  yb0 <- array( data=NA, dim=c( n0, nens))
  for (e in 1:nens) {
    r[] <- NA
    r[aix]  <- Xb[,e]
    yb0[,e] <- extract( r, cbind(xcoord,ycoord), method="bilinear")
    # if NAs, then extra-shot with nearest neighbour
    isna <- which( yb0[,e] < 0 | is.na( yb0[,e]) | is.nan( yb0[,e]))
    if ( length( isna)>0) 
      yb0[isna,e] <- extract( r, cbind( xcoord[isna], ycoord[isna]))
    rm(isna)
  }
  # index over yb with all the ensemble members finite and not NAs
  ix <- which( apply( yb0, MAR=1, FUN = function(x) 
                     {length( which( !is.na(x) & !is.nan(x)))}) == nens )
  n1 <- length(ix)
  yb <- array( data=NA, dim=c( n1, nens))
  for (e in 1:nens) yb[,e] <- yb0[ix,e]
  rm( yb0)
  yb[yb<inf] <- inf
  return( list( yb=yb, ix=ix))
}

