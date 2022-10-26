#+
clean_and_smooth_the_field <- function( r, threshold=NA, s=NULL, 
                                        val_obs=NULL, x_obs=NULL, y_obs=NULL, 
                                        lambda=2,
                                        area_small_clumps=100) {
#------------------------------------------------------------------------------
# Input:
# r, is a raster
# threshold, is a scalar value with the threshold for rain/no-rain
#------------------------------------------------------------------------------

  r_area_one_cell <- res(r)[1] * res(r)[2]
  ncells_small_clumps <- ceiling( area_small_clumps / r_area_one_cell)
  # 1st step: remove small precipitation values and replace NAs with 0s
  values <- getValues(r)
  if ( !is.na( threshold)) values[values<threshold] <- 0
  if ( any( is.na( values)))   values[is.na(values)]  <- 0

  # 2nd step: smooth 
  r[] <- values
  t   <- r - s
  t[] <- gauss2dsmooth( x=as.matrix(t), lambda=lambda, nx=ncol(r), ny=nrow(r))
  r   <- s + t
  values <- getValues(r)
  if ( !is.na( threshold)) values[values<threshold] <- 0
  r[] <- values

  # 3rd step: remove small clumps of precipitation or those that do not include wet obs (and they are not in s, the background field)
  values_or <- getValues(r)
  values_bk <- getValues(s)
  values <- abs( values_or - values_bk)
  values[which( !is.na(values) & values<threshold)] <- NA
  r[] <- values
  rclump <- clump(r)
  fr <- freq(rclump)
  oclump <- extract( rclump, cbind( x_obs, y_obs))
  ix <- which( !is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=ncells_small_clumps) | !(fr[,1] %in% oclump)) )
  if ( length( ixx <- which( getValues(rclump) %in% fr[ix,1])) > 0)
    values_or[ixx] <- values_bk[ixx]
  r[] <- values_or

  # return
  r
}
