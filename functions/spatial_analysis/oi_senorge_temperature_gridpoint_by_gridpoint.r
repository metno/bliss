#+ ensemble optimal interpolation 
oi_senorge_temperature_gridpoint_by_gridpoint <- function( i,
                                                           corr = c("soar","gaussian","linear"),
#                                                           dh = 10000,
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
    di <- array( data=envtmp$D[ixa,], dim=c(p,1))

    dist_1d    <- envtmp$nn2$nn.dists[i,aux]
    diffz_1d   <- envtmp$obs_z[i] - envtmp$obs_z[ixa]
    difflaf_1d <- envtmp$obs_laf[i] - envtmp$obs_laf[ixa]

    dist_2d    <- envtmp$dist2[ixa,ixa]
    diffz_2d   <- envtmp$dist2_z[ixa,ixa]
    difflaf_2d <- envtmp$dist2_laf[ixa,ixa]

    # compute correlations
    rloc <- corr1d( values=dist_1d, values2=diffz_1d, values3=difflaf_1d, par=envtmp$par[i,], label=corr) 
    S    <- corr2d( values=dist_2d, values2=diffz_2d, values3=difflaf_2d, par=envtmp$par[i,], label=corr)

    # OI equations
    SRinv <- chol2inv(chol( (S+diag(x=envtmp$eps2[i],p)) ))
    SRinv_di <- crossprod(SRinv,di)       
    Ea <- envtmp$Eb[i,] + crossprod( rloc, SRinv_di)
    if (idi) xidi <- sum( rloc * as.vector(rowSums(SRinv)))
  }

  #
  return( c( Ea, xidi))
}

