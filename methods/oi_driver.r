#+
oi_driver <- function( argv, y_env, fg_env, env) {

  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # number of observations
  env$p_dim <- length(y_env$yo$x)

  # initialization
  env$Xa <- array( data=NA, dim=c( env$ngrid, env$a_dim))
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$a_dim)) 
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$a_dim)) 
  r <- env$rmaster

  # loop over ensemble members
  for (e in 1:env$k_dim) {
    i <- fg_env$ixs[e]
    rb <- subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i])
    if ( e == 1) {
      aux <- getValues(rb)
      if (length( ixb <- which( !is.na( aux))) == 0) boom()
      env$m_dim <- length( ixb)
      xy <- xyFromCell( rb, 1:ncell(rb))
      envtmp$x <- xy[ixb,1]
      envtmp$y <- xy[ixb,2]
      envtmp$eps2 <- rep( argv$eps2, env$m_dim)
      envtmp$nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
                         query = cbind(envtmp$x,envtmp$y), 
                         k = argv$oi_pmax, searchtype = "radius", 
                         radius = (7*argv$oi_dh))
      cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
      cat( paste( "number of observations, p dim >", env$p_dim, "\n"))
    }
    envtmp$xb <- getValues(rb)[ixb]
    envtmp$yb <- extract( rb, cbind(y_env$yo$x, y_env$yo$y))
    # run oi gridpoint by gridpoint (idi only the first time)
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( oi_basic_gridpoint_by_gridpoint,
                          1:env$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          pmax=argv$oi_pmax,
                          corr=argv$oi_corrfun,
                          dh=argv$oi_dh,
                          idi=(e==1)))
      
    # no-multicores
    } else {
      res <- t( mapply( oi_basic_gridpoint_by_gridpoint,
                        1:env$m_dim,
                        SIMPLIFY=T,
                        pmax=argv$oi_pmax,
                        corr=argv$oi_corrfun,
                        dh=argv$oi_dh,
                        idi=(e==1)))
    }

    # store analysis
    if (e==1) { env$Xa[ixb,e] <- res[,1] } else { env$Xa[ixb,e] <- res[1,] }
    # Safe checks
    if (!is.na(argv$oi_range[1])) 
      env$Xa[env$Xa[,e]<argv$oi_range[1],e] <- argv$oi_range[1]  
    if (!is.na(argv$oi_range[2])) 
      env$Xa[env$Xa[,e]>argv$oi_range[2],e] <- argv$oi_range[2]  
    
    # extract point values (crossvalidation)
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) { 
      y_env$yov$value_a[,e] <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
      if (!is.na(argv$oi_range[1])) 
        y_env$yov$value_a[y_env$yov$value_a[,e]<argv$oi_range[1],e] <- argv$oi_range[1]  
      if (!is.na(argv$oi_range[2])) 
        y_env$yov$value_a[y_env$yov$value_a[,e]>argv$oi_range[2],e] <- argv$oi_range[2]  
    }
    
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y))
    if (!is.na(argv$oi_range[1])) 
      y_env$yo$value_a[y_env$yo$value_a[,e]<argv$oi_range[1],e] <- argv$oi_range[1]  
    if (!is.na(argv$oi_range[2])) 
      y_env$yo$value_a[y_env$yo$value_a[,e]>argv$oi_range[2],e] <- argv$oi_range[2] 
    
    # IDI (just do it once)
    if ( e == 1 ) {
      r[] <- NA
      r[ixb] <- res[,2]
      if (env$cv_mode | env$cv_mode_random) 
        y_env$yov$idi <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
      y_env$yo$idi <- extract( r, cbind(  y_env$yo$x, y_env$yo$y))
    }
     
    cat("\n")
  } # end loop over ensemble members
}
