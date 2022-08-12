#+ Merge different observation data sources
corenks_mergeobs <- function( argv, y_env, u_env, env) {
#
# output
#  env$mergeobs$value - merged obs values (cover the same region as the "background")
#  env$mergeobs$idi - integral data influence (=1 dense obs region; =0 obs void region; set to NAs where number of nearby obs is less than pmax)
#  env$mergeobs$obs_errvar - estimated "observation" error variance at gridpoints (covers the same region as idi)
#  env$mergeobs$mergedobs_errvar - estimated merged observation error variance at gridpoints (covers the same region as idi)
# similar output for y_env
#  y_env$yo(v)$mergeobs_value 
#  y_env$yo(v)$mergeobs_idi 
#  y_env$yo(v)$mergeobs_obs_errvar 
#  y_env$yo(v)$mergeobs_mergedobs_errvar 
#  
#..............................................................................

  t0a <- Sys.time()

  cat( "-- corenks mergeobs --\n")

  # number of observations
  env$p_dim <- length(y_env$yo$x)

  # initialization
  if (env$cv_mode | env$cv_mode_random) 
    yov_value_a <- array( data=NA, dim=c( y_env$yov$n, 1)) 
  yo_value_a <- array( data=NA, dim=c( y_env$yo$n, 1)) 
  r <- env$rmaster

  m_dim <- ncell(r)
  xy <- xyFromCell( r, 1:ncell(r))
  envtmp$m_dim <- m_dim
  envtmp$x <- xy[,1]
  envtmp$y <- xy[,2]
  envtmp$eps2 <- rep( argv$corenks_mergeobs_eps2, m_dim)
  envtmp$nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
                     query = cbind(envtmp$x,envtmp$y), 
                     k = argv$corenks_mergeobs_pmax, searchtype = "radius", 
                     radius = (7*argv$corenks_mergeobs_dh))
  cat( paste( " number of grid points, m dim >", m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))

  envtmp$xb <- getValues( u_env$uo[[1]]$r_main)
  envtmp$yb <- extract( u_env$uo[[1]]$r_main, cbind(y_env$yo$x, y_env$yo$y))

  # run oi gridpoint by gridpoint (idi the first time only)
  if (!is.na(argv$cores)) {
    res <- t( mcmapply( oi_basic_gridpoint_by_gridpoint,
                        1:m_dim,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        pmax=argv$corenks_mergeobs_pmax,
                        corr=argv$corenks_mergeobs_corrfun,
                        dh=argv$corenks_mergeobs_dh,
                        idi=T,
                        uncertainty=T))
    
  # no-multicores
  } else {
    res <- t( mapply( oi_basic_gridpoint_by_gridpoint,
                      1:m_dim,
                      SIMPLIFY=T,
                      pmax=argv$corenks_mergeobs_pmax,
                      corr=argv$corenks_mergeobs_corrfun,
                      dh=argv$corenks_mergeobs_dh,
                      idi=T,
                      uncertainty=T))
  }
  # res: 1 xa, 2 xidi, 3 o_errvar, 4 xa_errvar

  # Safe checks
  if (!is.na(argv$corenks_mergeobs_range[1])) 
    res[,1][res[,1]<argv$corenks_mergeobs_range[1]] <- argv$corenks_mergeobs_range[1]  
  if (!is.na(argv$corenks_mergeobs_range[2])) 
    res[,1][res[,1]>argv$corenks_mergeobs_range[2]] <- argv$corenks_mergeobs_range[2]  

  # output
  env$mergeobs <- list()
  env$mergeobs$value <- res[,1]
  env$mergeobs$idi    <- res[,2]
  env$mergeobs$obs_errvar <- res[,3]
  env$mergeobs$mergedobs_errvar <- res[,4]
  
  # extract point values (crossvalidation)
  if (env$cv_mode | env$cv_mode_random) { 
    r[] <- env$mergeobs$value
    y_env$yov$mergeobs_value <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
    r[] <- env$mergeobs$idi 
    y_env$yov$mergeobs_idi <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
    r[] <- env$mergeobs$obs_errvar 
    y_env$yov$mergeobs_obs_errvar <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
    r[] <- env$mergeobs$mergedobs_errvar 
    y_env$yov$mergeobs_mergedobs_errvar <- extract( r, cbind( y_env$yov$x, y_env$yov$y))
  }
    
  # extract point values (analysis)
  r[] <- env$mergeobs$value
  y_env$yo$mergeobs_value <- extract( r, cbind( y_env$yo$x, y_env$yo$y))
  r[] <- env$mergeobs$idi 
  y_env$yo$mergeobs_idi <- extract( r, cbind( y_env$yo$x, y_env$yo$y))
  r[] <- env$mergeobs$obs_errvar 
  y_env$yo$mergeobs_obs_errvar <- extract( r, cbind( y_env$yo$x, y_env$yo$y))
  r[] <- env$mergeobs$mergedobs_errvar 
  y_env$yo$mergeobs_mergedobs_errvar <- extract( r, cbind( y_env$yo$x, y_env$yo$y))

  t1a <- Sys.time()
  cat( paste( "\n", "mergeobs .", "dim =", m_dim, ".", "range =", 
              round( min( env$mergeobs$value, na.rm=T), 2), 
              round( max( env$mergeobs$value, na.rm=T), 2), "\n"))
    cat( paste( "total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))
    cat( "+---------------------------------+\n")
 
  cat("\n")

}
