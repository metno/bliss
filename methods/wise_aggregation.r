#+
wise_aggregation <- function( argv, y_env, env, 
                              plot=F,dir_plot=NA) {
#
#------------------------------------------------------------------------------

  t0a <- Sys.time()

  cat( "-- wise aggregation --\n")

  # set dyadic domain
  xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
  xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
  ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
  ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2
  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)
  n <- ceiling( log( max(nx,ny), 2))
  rdyad <- raster( extent( xmn, xmx, ymn, ymx), ncol=2**n, nrow=2**n, 
                   crs=argv$grid_master.proj4)
  rdyad[] <- 0
  sqrt_m_dim <- sqrt( env$m_dim)

  # initialization
  y_env$yo$value_fgsel_pp <- array( data=NA, dim=c( y_env$yo$n, env$k_dim)) 
  # distance used in aggregation to look up for neighbours is set to 
  #  the resolution of the master grid
  agg_dh <- max(res(env$rmaster))
  # set the max number of dyadic gridpoints to consider (close to a point)
  agg_pmax <- ceiling( 4 * agg_dh**2 / min(res(rdyad))**2)

  # coordinates of notNA points on the dyadic grid
  if (length( ixdyad <- which( !is.na( env$Xbpp_dyad[,1]))) == 0) boom()
  xydyad <- xyFromCell( rdyad, 1:ncell(rdyad))
  xdyad <- xydyad[ixdyad,1]
  ydyad <- xydyad[ixdyad,2]

  # get the list of points on the dyadic domain closest to observation locations
  nn2_yo <- nn2( cbind( xdyad, ydyad),
                 query = cbind( y_env$yo$x, y_env$yo$y), 
                 k = agg_pmax, searchtype = "radius", 
                 radius = agg_dh)

  # cv_mode 
  if (env$cv_mode | env$cv_mode_random) { 
    y_env$yov$value_fgsel_pp <- array( data=NA, dim=c( y_env$yov$n, env$k_dim)) 
    # get the list of points on the dyadic domain closest to cv-observation locations
    nn2_yov <- nn2( cbind( xdyad, ydyad),
                           query = cbind( y_env$yov$x, y_env$yov$y),
                           k = agg_pmax, searchtype = "radius", 
                           radius = agg_dh)
 
  # spatial analysis on the grid
  } else {
    env$Xbpp <- array( data=NA, dim=c( env$ngrid, env$k_dim))
    # get the list of points on the dyadic domain closest to points on the master grid
    if (length( ixm <- which( !is.na( env$mask))) == 0) boom()
    xy <- xyFromCell( env$rmaster, 1:ncell(env$rmaster))
    x <- xy[ixm,1]
    y <- xy[ixm,2]
    nn2_x <- nn2( cbind( xdyad, ydyad),
                  query = cbind( x, y), 
                  k = agg_pmax, searchtype = "radius", 
                  radius = agg_dh)
  }

  # loop over ensembles
  for (e in 1:env$k_dim) {
    if (e == 1) cat ("loop over ensembles ")
    cat(".")
    t0 <- Sys.time()

    # select the ensemble member
    xbpp_dyad <- env$Xbpp_dyad[ixdyad,e]

    # parallel
    if (!is.na(argv$cores)) {
      y_env$yo$value_fgsel_pp[,e] <- t( mcmapply( wise_aggregation_helper,
                                                  1:y_env$yo$n,
                                                  mc.cores=argv$cores,
                                                  SIMPLIFY=T,
                                                  MoreArgs = list(
                                                    q_prob=argv$wise_agg_qprob,
                                                    nn2=nn2_yo,
                                                    values=xbpp_dyad)))
    # no-multicores
    } else {
      y_env$yo$value_fgsel_pp[,e] <- t( mapply( wise_aggregation_helper,
                                                1:y_env$yo$n,
                                                SIMPLIFY=T,
                                                MoreArgs = list(
                                                  q_prob=argv$wise_agg_qprob,
                                                  nn2=nn2_yo,
                                                  values=xbpp_dyad)))
    }

    # cv-mode
    if (env$cv_mode | env$cv_mode_random) {
      # parallel
      if (!is.na(argv$cores)) {
        y_env$yov$value_fgsel_pp[,e] <- t( mcmapply( wise_aggregation_helper,
                                                     1:y_env$yov$n,
                                                     mc.cores=argv$cores,
                                                     SIMPLIFY=T,
                                                     MoreArgs = list(
                                                       q_prob=argv$wise_agg_qprob,
                                                       nn2=nn2_yov,
                                                       values=xbpp_dyad)))
      # no-multicores
      } else {
        y_env$yov$value_fgsel_pp[,e] <- t( mapply( wise_aggregation_helper,
                                                   1:y_env$yov$n,
                                                   SIMPLIFY=T,
                                                   MoreArgs = list(
                                                     q_prob=argv$wise_agg_qprob,
                                                     nn2=nn2_yov,
                                                     values=xbpp_dyad)))
      }

    # spatial analysis on the grid
    } else {
      # parallel
      if (!is.na(argv$cores)) {
        env$Xbpp[,e] <- t( mcmapply( wise_aggregation_helper,
                                     1:env$ngrid,
                                     mc.cores=argv$cores,
                                     SIMPLIFY=T,
                                     MoreArgs = list(
                                       q_prob=argv$wise_agg_qprob,
                                       nn2=nn2_x,
                                       values=xbpp_dyad)))
      # no-multicores
      } else {
        env$Xbpp[,e] <- t( mapply( wise_aggregation_helper,
                                   1:env$ngrid,
                                   SIMPLIFY=T,
                                   MoreArgs = list(
                                     q_prob=argv$wise_agg_qprob,
                                     nn2=nn2_x,
                                     values=xbpp_dyad)))
      }
    } # end if cv-mode or spatial analysis

  } # end loop over ensembles

  t1a <- Sys.time()
  cat( paste( "\n total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))
} # END FUNCTION

#+
wise_aggregation_helper <- function(i, q_prob, nn2, values) {
  if ( (p <- length( aux <- which(nn2$nn.idx[i,]!=0))) == 0) {
    res <- NA
  } else {
    # observations available
    ix <- nn2$nn.idx[i,aux]
    val <- values[ix]
    if ( length( ixgt0 <- which( val >=  y_env$rain)) < floor(p/2)) {
      res <- 0
    } else {
      res <- quantile( val[ixgt0], probs=q_prob, type=4, na.rm=T)
    } 
  }
  res
}
