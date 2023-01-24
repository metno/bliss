#+ ensemble optimal interpolation 
enoi_basic_xyzlaf_gridpoint_by_gridpoint<-function( i,
                                                    corr = c("soar","gaussian","linear"),
#                                                    dh = 10000,
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

    arr_for_corr1d <- array( data=NA, dim=c(p,3))
    arr_for_corr1d[,1] <- envtmp$nn2$nn.dists[i,aux]
    arr_for_corr1d[,2] <-   envtmp$obs_z[i] - envtmp$obs_z[ixa]
    arr_for_corr1d[,3] <- envtmp$obs_laf[i] - envtmp$obs_laf[ixa] 

    arr_for_corr2d <- array( data=NA, dim=c(p,p,3))
    arr_for_corr2d[,,1] <- outer(  envtmp$obs_x[ixa],  envtmp$obs_x[ixa],FUN="-")**2. + 
                           outer(  envtmp$obs_y[ixa],  envtmp$obs_y[ixa],FUN="-")**2
    arr_for_corr2d[,,2] <- outer(  envtmp$obs_z[ixa],  envtmp$obs_z[ixa],FUN="-")**2.
    arr_for_corr2d[,,3] <- outer(envtmp$obs_laf[ixa],envtmp$obs_laf[ixa],FUN="-")**2.
#    x <- envtmp$obs_x[ixa]
#    y <- envtmp$obs_y[ixa]
#    y <- envtmp$obs_z[ixa]
    di <- array( data=envtmp$D[ixa,], dim=c(p,length(envtmp$Eb[i,])))

    # compute correlations
    rloc <- corr1d( arr_for_corr1d, envtmp$par[i,], corr) 
    S    <- corr2d( arr_for_corr2d, envtmp$par[i,], corr)

    # OI equations
    SRinv <- chol2inv(chol( (S+diag(x=envtmp$eps2[i],p)) ))
    SRinv_di <- crossprod(SRinv,di)       
    Ea <- envtmp$Eb[i,] + crossprod( rloc, SRinv_di)
    if (idi) xidi <- sum( rloc * as.vector(rowSums(SRinv)))
  }

  #
  return( c( Ea, xidi))
}

