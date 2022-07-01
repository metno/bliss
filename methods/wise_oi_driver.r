#+
wise_oi_driver <- function( argv, y_env, fg_env, env) {

  cat( "-- wise oi --\n")

  # number of observations
  env$p_dim <- length(y_env$yo$x)

  # initialization
  Xa <- array( data=NA, dim=c( env$ngrid, env$a_dim))
  if (env$cv_mode | env$cv_mode_random) 
    yov_value_a <- array( data=NA, dim=c( y_env$yov$n, env$a_dim)) 
  yo_value_a <- array( data=NA, dim=c( y_env$yo$n, env$a_dim)) 
  r <- env$rmaster

  if (length( ixb <- which( !is.na( env$Xa[,1]))) == 0) boom()
  env$m_dim <- length( ixb)
  xy <- xyFromCell( r, 1:ncell(r))
  envtmp$x <- xy[ixb,1]
  envtmp$y <- xy[ixb,2]
  envtmp$eps2 <- rep( argv$wise_oi_eps2, env$m_dim)
  envtmp$nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
                     query = cbind(envtmp$x,envtmp$y), 
                     k = argv$wise_oi_pmax, searchtype = "radius", 
                     radius = (7*argv$wise_oi_dh))
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))

  # loop over ensemble members
  for (e in 1:env$k_dim) {
    envtmp$xb <- env$Xa[ixb,e]
    r[] <- env$Xa[,e]
    envtmp$yb <- extract( r, cbind(y_env$yo$x, y_env$yo$y))

    # run oi gridpoint by gridpoint (idi the first time only)
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( oi_basic_gridpoint_by_gridpoint,
                          1:env$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          pmax=argv$wise_oi_pmax,
                          corr=argv$wise_oi_corrfun,
                          dh=argv$wise_oi_dh,
                          idi=(e==1)))
      
    # no-multicores
    } else {
      res <- t( mapply( oi_basic_gridpoint_by_gridpoint,
                        1:env$m_dim,
                        SIMPLIFY=T,
                        pmax=argv$wise_oi_pmax,
                        corr=argv$wise_oi_corrfun,
                        dh=argv$wise_oi_dh,
                        idi=(e==1)))
    }

    # store analysis
    if (e==1) { env$Xa[ixb,e] <- res[,1] } else { env$Xa[ixb,e] <- res[1,] }
    # Safe checks
    if (!is.na(argv$wise_oi_range[1])) 
      env$Xa[env$Xa[,e]<argv$wise_oi_range[1],e] <- argv$wise_oi_range[1]  
    if (!is.na(argv$oi_range[2])) 
      env$Xa[env$Xa[,e]>argv$wise_oi_range[2],e] <- argv$wise_oi_range[2]  
    
    # extract point values (crossvalidation)
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) { 
      y_env$yov$value_a[,e] <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
      if (!is.na(argv$wise_oi_range[1])) 
        y_env$yov$value_a[y_env$yov$value_a[,e]<argv$wise_oi_range[1],e] <- argv$wise_oi_range[1]  
      if (!is.na(argv$wise_oi_range[2])) 
        y_env$yov$value_a[y_env$yov$value_a[,e]>argv$wise_oi_range[2],e] <- argv$wise_oi_range[2]  
    }
    
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y))
    if (!is.na(argv$wise_oi_range[1])) 
      y_env$yo$value_a[y_env$yo$value_a[,e]<argv$wise_oi_range[1],e] <- argv$wise_oi_range[1]  
    if (!is.na(argv$wise_oi_range[2])) 
      y_env$yo$value_a[y_env$yo$value_a[,e]>argv$wise_oi_range[2],e] <- argv$wise_oi_range[2] 
    
    # IDI (just do it once)
    if ( e == 1 ) {
      r[] <- NA
      r[ixb] <- res[,2]
      env$Xidi <- getValues(r)
      if (env$cv_mode | env$cv_mode_random) 
        y_env$yov$idi <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
      y_env$yo$idi <- extract( r, cbind(  y_env$yo$x, y_env$yo$y))
    }
     
    cat("\n")
  } # end loop over ensemble members
}
