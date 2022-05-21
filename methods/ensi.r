#+
ensi <- function( argv, y_env, fg_env, env) {
  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)
  # ---~--------------
  # define constants
  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length( getValues(env$rmaster))
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))

  env$Xa <- array( data=NA, dim=c( env$ngrid, env$a_dim))
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$a_dim)) 
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$a_dim)) 
  r <- env$rmaster

  for (e in 1:env$k_dim) {
    i <- fg_env$ixs[e]
    if ( e == 1) {
      aux <- getValues( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]))
      if (length( ixb <- which( !is.na( aux))) == 0) boom()
      env$m_dim <- length( ixb)
      xy <- xyFromCell(fg_env$fg[[fg_env$ixf[i]]]$r_main,1:ncell(fg_env$fg[[fg_env$ixf[i]]]$r_main))
      envtmp$x <- xy[ixb,1]
      envtmp$y <- xy[ixb,2]
      envtmp$eps2 <- rep( 0.1, env$m_dim)
      envtmp$nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
                         query = cbind(envtmp$x,envtmp$y), 
                         k = argv$pmax, searchtype = "radius", 
                         radius = 70000)
    }
    envtmp$xb <- getValues(fg_env$fg[[fg_env$ixf[i]]]$r_main)[ixb]
    envtmp$yb <- extract( fg_env$fg[[fg_env$ixf[i]]]$r_main, cbind( y_env$yo$x, y_env$yo$y))
    if (!is.na(argv$cores)) {
      env$Xa[ixb,e] <- t( mcmapply( oi_basic_gridpoint_by_gridpoint,
                          1:env$m_dim,
                          mc.cores=argv$cores,
                          SIMPLIFY=T,
                          pmax=argv$pmax,
                          corr=argv$corrfun))
    # no-multicores
    } else {
      env$Xa[ixb,e] <- t( mapply( oi_basic_gridpoint_by_gridpoint,
                        1:env$m_dim,
                        SIMPLIFY=T,
                        pmax=argv$pmax,
                        corr=argv$corrfun))[,1]
    }
    
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,a] <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
    y_env$yo$value_a[,a]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y))
  }
}
