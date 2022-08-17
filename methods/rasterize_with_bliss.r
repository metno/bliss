#+
rasterize_with_bliss <- function( argv, y_env, env) {
#==============================================================================
  cat("-~- Rasterize -~-\n")

  # default, rasterize using mean values
  r <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=env$rmaster, field=y_env$yo$value, fun=mean, na.rm=T )
  env$Xa <- getValues(r)
  y_env$yo$value_a <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
  if (env$cv_mode | env$cv_mode_random)
    y_env$yov$value_a <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
  
  # diagnostic info, standard deviation at each cell
  r <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=env$rmaster, field=y_env$yo$value, fun=sd, na.rm=T )
  env$Xa_sd <- getValues(r)
  y_env$yo$value_a_sd <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
  if (env$cv_mode | env$cv_mode_random)
    y_env$yov$value_a_sd <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")

  # diagnostic info, number of observations in each cell
  r <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=env$rmaster, field=y_env$yo$value, fun=function(x,na.rm){length(na.omit(x))} )
  env$Xa_n <- getValues(r)
  y_env$yo$value_a_n <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
  if (env$cv_mode | env$cv_mode_random)
    y_env$yov$value_a_n <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")

  # mask out cell with less than a predefined number in observations in them
  ix_raster <- integer(0)
  if (!is.na(argv$rasterize_nmin)) {
    ix_raster <- which( !is.na(env$Xa_n) & env$Xa_n<argv$rasterize_nmin)
    # update the grid and the vectors
    if (length(ix_raster) > 0) {
      env$Xa[ix_raster]<-NA
      env$Xa_sd[ix_raster]<-NA
      r[] <- env$Xa
      y_env$yo$value_a <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
      if (env$cv_mode | env$cv_mode_random)
        y_env$yov$value_a <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
      r[] <- env$Xa_sd
      y_env$yo$value_a_sd <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
      if (env$cv_mode | env$cv_mode_random)
        y_env$yov$value_a_sd <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    }
    cat( paste("    number of boxes with observations within=", length(which(!is.na(env$Xa_n))),"\n"))
    cat( paste("    number of masked boxes with observations within=", length(ix_raster),"\n"))
    cat( paste(" -> number of unmasked boxes=",length(which(!is.na(env$Xa))),"\n"))
  }

  # additional output based on alternative aggregation method(s)
  if (!any(is.na(argv$rasterize_q))) {
    env$Xa_q <- array( data=NA, dim=c(length(env$Xa_n),length(argv$rasterize_q)))
    y_env$yo$value_a_q <- array( data=NA, dim=c( y_env$yo$n, length(argv$rasterize_q)))
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a_q <- array( data=NA, dim=c( y_env$yov$n, length(argv$rasterize_q)))
    for (i in 1:length(argv$rasterize_q)) {
      r <- rasterize( x=cbind(y_env$yo$x,y_env$yo$y), y=env$rmaster, 
                      field=y_env$yo$value, 
                      fun=function(x,na.rm){
                       quantile(na.omit(x),prob=argv$rasterize_q[i],type=4)} )
      env$Xa_q[,i] <- getValues(r)
      if ( length(ix_raster) > 0) env$Xa_q[ix_raster,i]<-NA
      r[] <- env$Xa_q[,i]
      y_env$yo$value_a_q[,i] <- extract( r, cbind( y_env$yo$x, y_env$yo$y), method="simple")
      if (env$cv_mode | env$cv_mode_random)
        y_env$yov$value_a_q[,i] <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    }
  }

}
