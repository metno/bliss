#+ ensemble optimal interpolation (Faster alternative) 
oi_senorge_precipitation_gridpoint_by_gridpoint<-function( i,
                                                           corr = "soar",
                                                           dh = 10000,
                                                           safecheck=F,
                                                           idi = F,
                                                           uncertainty = F) {
# returned values: analysis, IDI
# Faster alternative to use when the number of observations is such that one
# can keep the matrices in the observation space into memory
# taken from enoi_basicFaster_gridpoint_by_gridpoint.r
#------------------------------------------------------------------------------
  xidi <- NA; #o_errvar <- NA; xa_errvar <- NA

#  if( i>100 & i%%(round(envtmp$m_dim/10)) == 0) cat(".")


  # define vectors
  dist <- sqrt((envtmp$x[i]-envtmp$obs_x)**2 + (envtmp$y[i]-envtmp$obs_y)**2)

  # compute correlations
  rloc <- corr1d( dist, dh, corr) 


  # OI equations
  Ea <- envtmp$Eb[i,] + crossprod( rloc, envtmp$SRinv_di)
  if (safecheck) {
    res <- Ea - envtmp$Eb[i,]
    if (sum(rloc)>0) {
      res_alt <- as.numeric( crossprod( rloc/sum(rloc), envtmp$di))
      if ( abs(res) > abs(res_alt) & sign(res)==sign(res_alt)) {
        Ea <- envtmp$Eb[i,] + res_alt 
      }
    }
  }
  if (idi) xidi <- sum( rloc * as.vector(rowSums(envtmp$SRinv)))

  #
  return( c( Ea, xidi))
}

