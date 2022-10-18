#+
clean_and_smooth_the_field <- function( r, y_env, s=NULL, 
                                        val_obs=NULL, x_obs=NULL, y_obs=NULL, 
                                        lambda=2,
                                        area_small_clumps=100) {
#------------------------------------------------------------------------------

  # 1st step: remove small precipitation values and replace NAs with 0s
  values <- getValues(r)
  if ( !is.na( y_env$rain)) values[values<y_env$rain] <- 0
  if ( any( is.na( values)))   values[is.na(values)]  <- 0

  # 2nd step: smooth 
  r[] <- values
  t   <- r - s
  t[] <- gauss2dsmooth( x=as.matrix(t), lambda=lambda, nx=ncol(r), ny=nrow(r))
  r   <- s + t
  values <- getValues(r)
  if ( !is.na( y_env$rain)) values[values<y_env$rain] <- 0
  r[] <- values

  # 3rd step: remove small clumps of precipitation or those that do not include wet obs (and they are not in s, the background field)
  values_or <- getValues(r)
  values <- abs( getValues(r) - getValues(s))
  values[which( !is.na(values) & values<y_env$rain)] <- NA
  r[] <- values
  rclump <- clump(r)
  fr <- freq(rclump)
#  ixwet <- which(val_obs<y_env$rain)
  oclump <- extract( rclump, cbind( x_obs, y_obs))
  ix <- which( !is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=area_small_clumps) | !(fr[,1] %in% oclump)) )
  values_or[which( getValues(rclump) %in% fr[ix,1])] <- 0
  r[] <- values_or

  # return
  r
}
