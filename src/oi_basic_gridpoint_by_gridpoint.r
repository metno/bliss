#+ optimal interpolation 
oi_basic_gridpoint_by_gridpoint<-function( i,
                                           corr="soar",
                                           pmax,
                                           isolation_factor=7,
                                           dh=10000,
                                           idi=T,
                                           uncertainty=F) {
# How to call this:
#  assume that: r is a raster with the background; 
#               y_env$yo$x/y/value are the observations
#  m_dim <- ncell(r)
#  xy <- xyFromCell( r, 1:ncell(r))
#  envtmp$m_dim <- m_dim
#  envtmp$x <- xy[,1]
#  envtmp$y <- xy[,2]
#  envtmp$eps2 <- rep( 0.5, m_dim)
#  envtmp$nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
#                     query = cbind(envtmp$x,envtmp$y), 
#                     k = argv$mergeobs_pmax, searchtype = "radius", 
#                     radius = (7*argv$mergeobs_dh))
#
# returned values: analysis, IDI, observation error var, analysis error var
#------------------------------------------------------------------------------
  xa <- NA; xidi <- NA; o_errvar <- NA; xa_errvar <- NA
  
  # print some stuff now and then
  if( i %% (round(envtmp$m_dim/10)) == 0) cat(".")

  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    xa <- envtmp$xb[i]
  } else {

    # selected observations
    ixa  <- envtmp$nn2$nn.idx[i,aux]

    # define vectors
    dist <- envtmp$nn2$nn.dists[i,aux]
    x <- y_env$yo$x[ixa]
    y <- y_env$yo$y[ixa]
    di <- y_env$yo$value[ixa] - envtmp$yb[ixa]
    eps2 <- envtmp$eps2[i]

    # compute correlations
    rloc <- corr1d( dist, dh, corr) 
    S    <- corr2d( cbind(x,y), dh, corr)

    # OI equations
    SRinv <- chol2inv(chol( (S+diag(x=eps2,p)) ))
    SRinv_di <- crossprod(SRinv,di)       
    xa <- envtmp$xb[i] + sum( rloc * as.vector(SRinv_di))
    if (idi) xidi <- sum(rloc*as.vector(rowSums(SRinv)))
    if (uncertainty) {
      o_errvar  <- mean( di * ( di - crossprod(S,SRinv_di)))
      xa_errvar <- ( o_errvar/ eps2) * 
                   ( 1 - sum( as.vector( crossprod( rloc, SRinv)) * rloc))
    }
  }

  # Exit
  return( c( xa, xidi, o_errvar, xa_errvar))
}

