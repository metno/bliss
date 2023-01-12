#+ ensemble optimal interpolation 
enoi_basic_gridpoint_by_gridpoint<-function( i,
                                             corr = "soar",
                                             dh = 10000,
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
    di <- array( data=envtmp$D[ixa,], dim=c(p,length(envtmp$Eb[i,])))
    eps2 <- envtmp$eps2[i]

    # compute correlations
    rloc <- corr1d( dist, dh, corr) 
    S    <- corr2d( cbind(x,y), dh, corr)

    # OI equations
    SRinv <- chol2inv(chol( (S+diag(x=eps2,p)) ))
    SRinv_di <- crossprod(SRinv,di)       
    Ea <- envtmp$Eb[i,] + crossprod( rloc, SRinv_di)
    if (idi) xidi <- sum( rloc * as.vector(rowSums(SRinv)))
  }

  #
  return( c( Ea, xidi))
}

