#+ 
enoi_Evensen2003_gridpoint_by_gridpoint <- function( i,
                                                     corr = "soar",
                                                     dh = 10000,
                                                     dh_loc = 10000,
                                                     alpha = 0.5,
                                                     k_dim_corr=10,
                                                     idi = F) {
# returned values: analysis, observation error var, analysis error var
# note: the call ‘t(x) %*% y’ (‘crossprod’) or ‘x %*% t(y)’ (‘tcrossprod’)
#------------------------------------------------------------------------------
  xidi <- NA

  if( i%%(round(envtmp$m_dim/10)) == 0) cat(".")

  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    Ea <- envtmp$Eb[i,]
  } else {

    # available observations
    ixa  <- envtmp$nn2$nn.idx[i,aux]

    di <- array( data=envtmp$D[ixa,], dim=c(p,envtmp$k_dim))
    if (envtmp$k_dim>k_dim_corr) {
      ixb <- order( colMeans( abs(di)))[1:k_dim_corr]
    } else {
      ixb <- 1:envtmp$k_dim
      k_dim_corr <- envtmp$k_dim
    }

    # define vectors
    dist <- envtmp$nn2$nn.dists[i,aux]
    x <- envtmp$obs_x[ixa]
    y <- envtmp$obs_y[ixa]
#    di <- array( data=envtmp$D[ixa,], dim=c(p,envtmp$k_dim))
    Zi <- array( data=envtmp$Z[i,ixb],   dim=c(1,k_dim_corr))
    Yi <- array( data=envtmp$Y[ixa,ixb], dim=c(p,k_dim_corr))
    eps2 <- envtmp$eps2[i]

    # localization
    loc1d <- corr1d( values=dist, par=dh_loc, label=corr)
    loc2d <- corr2d( values=outer(x,x,FUN="-")**2.+outer(y,y,FUN="-")**2., par=dh_loc, label=corr)

    # static correlations
    rloc <- corr1d( values=dist, par=dh, label=corr)
    Cyy_s <- corr2d( values=outer(x,x,FUN="-")**2.+outer(y,y,FUN="-")**2., par=dh, label=corr)

    # combine static and dynamic correlations
    Cxy <- alpha * rloc + (1-alpha) * loc1d * tcrossprod( Zi, Yi)
    Cyy <- alpha * Cyy_s + (1-alpha) * loc2d * tcrossprod( Yi, Yi)

    #
    CyyCdd_inv <- chol2inv( chol( (Cyy+diag(x=eps2,p)) ))
    CyyCdd_inv_di <- crossprod( CyyCdd_inv, di)       
    Ea <- envtmp$Eb[i,] + crossprod( t(Cxy), CyyCdd_inv_di)

    if (idi) xidi <- sum( Cxy * as.vector(rowSums(CyyCdd_inv)))
  }
  return( c( Ea, xidi))
}

