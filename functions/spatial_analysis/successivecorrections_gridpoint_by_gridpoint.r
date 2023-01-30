#+ ensemble optimal interpolation 
successivecorrections_gridpoint_by_gridpoint<-function( i,
                                             corr = "soar",
                                             dh = 10000,
                                             nSCloops = 100,
                                             idi = F,
                                             uncertainty = F) {
# returned values: analysis, observation error var, analysis error var
#------------------------------------------------------------------------------
  xidi <- NA; #o_errvar <- NA; xa_errvar <- NA

  if( i%%(round(envtmp$m_dim/10)) == 0) cat(".")

  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    Ea <- envtmp$Eb[i,]
  } else {

    # observations available
    ixa  <- envtmp$nn2$nn.idx[i,aux]

    # define vectors
    dist <- envtmp$nn2$nn.dists[i,aux]
    x <- envtmp$obs_x[ixa]
    y <- envtmp$obs_y[ixa]
    dist2 <- outer(x,x,FUN="-")**2. + outer(y,y,FUN="-")**2.
    di <- array( data=envtmp$D[ixa,], dim=c(p,length(envtmp$Eb[i,])))
    eps2 <- envtmp$eps2[i]

    # compute correlations
    rloc <- corr1d( values=dist, dh, corr) 
    S    <- corr2d( values=dist2, dh, corr)

    # SC equations
    SR <- S + diag(x=eps2,p)
    Minv <- 1 / rowSums(abs(SR))
    I_SR_Minv <- diag(p) - t( t(SR) * Minv )
    di_loop <- di
    for (sc in 1:nSCloops)  di_loop <-crossprod( t(I_SR_Minv), di_loop)
    di_loop <- di + di_loop
    Ea <- envtmp$Eb[i,] + crossprod( rloc * (Minv*di_loop))
    if (idi) {
      di_loop <- rep(1,p)
      for (sc in 1:nSCloops) 
        di_loop <-as.vector(crossprod(t(I_SR_Minv),as.vector(di_loop)))
      di_loop <- rep(1,p) + di_loop
      xidi <- sum( rloc * (Minv*di_loop))
    }
  }

  #
  return( c( Ea, xidi))
}

