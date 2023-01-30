#+
sc_driver <- function( argv, y_env, fg_env, env) {
#------------------------------------------------------------------------------
  t0a <- Sys.time()
  cat( "-- OI --\n")

  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # number of observations
  env$p_dim <- y_env$yo$n

  # initialization
  env$Xa   <- array( data=NA, dim=c( env$ngrid, env$k_dim))
  env$Xidi <- array( data=NA, dim=c( env$ngrid, 1))
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, env$k_dim)) 
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, env$k_dim)) 
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
      envtmp$m_dim <- length( ixb)
      envtmp$Eb <- array( data=NA, dim=c(env$m_dim,env$k_dim))
      envtmp$HE <- array( data=NA, dim=c(env$p_dim,env$k_dim))
      envtmp$x <- xy[ixb,1]
      envtmp$y <- xy[ixb,2]
      envtmp$eps2 <- rep( argv$eps2, env$m_dim)
      envtmp$nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
                         query = cbind(envtmp$x,envtmp$y), 
                         k = argv$pmax, searchtype = "radius", 
                         radius = (7*argv$dh))
      cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
      cat( paste( "number of observations, p dim >", env$p_dim, "\n"))
    }
    envtmp$Eb[,e] <- getValues(rb)[ixb]
    envtmp$HE[,e] <- extract( rb, cbind(y_env$yo$x, y_env$yo$y))
    # perturb the observations
    obs <- y_env$yo$value
    if (argv$obs_perturb) {
      if (!is.na(y_env$rain)) {
        ixo <- which( y_env$yo$value >= y_env$rain)
      } else {
        ixo <- 1:y_env$yo$n
      }
      if ( length(ixo) > 0)
        obs[ixo] <- obs[ixo] * runif( length(ixo), min=argv$obs_perturb_rmin, max=argv$obs_perturb_rmax)
      rm(ixo) 
    }
    envtmp$D[,e] <- obs - extract( rb, cbind(y_env$yo$x, y_env$yo$y))
    rm(obs)
  }
  rm(rb)
#  envtmp$D <- y_env$yo$value - envtmp$HE
  envtmp$HE <- NULL
  envtmp$obs_x <- y_env$yo$x
  envtmp$obs_y <- y_env$yo$y
  # run oi gridpoint by gridpoint (idi only the first time)
  if (!is.na(argv$cores)) {
    res <- t( mcmapply( successivecorrections_gridpoint_by_gridpoint,
                        1:env$m_dim,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        MoreArgs = list( corr=argv$corrfun,
                                         dh=argv$dh,
                                         nSCloops=argv$nSCloops,
                                         idi=T)))
      
  # no-multicores
  } else {
    res <- t( mapply( successivecorrections_gridpoint_by_gridpoint,
                      1:env$m_dim,
                      SIMPLIFY=T,
                      MoreArgs = list( corr=argv$corrfun,
                                       dh=argv$dh,
                                       nSCloops=argv$nSCloops,
                                       idi=T)))
  }
  env$Xa[ixb,1:env$k_dim] <- res[,1:env$k_dim]
  env$Xidi[ixb,1] <- res[,env$k_dim+1]
  rm(res)
  envtmp$Eb <- NULL
  if (any(is.na(env$Xa))) env$Xa[is.na(env$Xa)] <- 0
  # Safe checks
  if (!is.na(argv$analysis_range[1])) 
    env$Xa[env$Xa<argv$analysis_range[1]] <- argv$analysis_range[1]  
  if (!is.na(argv$analysis_range[2])) 
    env$Xa[env$Xa>argv$analysis_range[2]] <- argv$analysis_range[2] 

  # Save results
  r[] <- NA
  for (e in 1:env$k_dim) {
    r[] <- env$Xa[,e]
    if (env$cv_mode | env$cv_mode_random) 
      y_env$yov$value_a[,e] <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
    # extract point values (analysis)
    y_env$yo$value_a[,e]  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y), method="simple")
  }
  r[] <- env$Xidi[,1]
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$idi <- extract( r, cbind( y_env$yov$x, y_env$yov$y), method="simple")
  # extract point values (analysis)
  y_env$yo$idi  <- extract( r, cbind(  y_env$yo$x, y_env$yo$y), method="simple")
  
  t1a <- Sys.time()
  cat( paste( "OI total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))
}
