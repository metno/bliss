#+  Change-Of-Resolution Ensemble Rauch-Tung-Striebel smoother DownSweep (coarse-to-fine)
corensi_down_gridpoint_by_gridpoint<-function( i,
                                               idi = F) {
# returned values: analysis, observation error var, analysis error var
# note: the call ‘t(x) %*% y’ (‘crossprod’) or ‘x %*% t(y)’ (‘tcrossprod’)
#------------------------------------------------------------------------------
  xidi <- NA; #o_errvar <- NA; xa_errvar <- NA

  if( i%%(round(envtmp$m_dim/10)) == 0) cat(".")

  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    Ea <- envtmp$Ea[i,]
  } else {

    # available observations
    ixa  <- envtmp$nn2$nn.idx[i,aux]

    # define vectors
    dist <- envtmp$nn2$nn.dists[i,aux]
    x <- envtmp$obs_x[ixa]
    y <- envtmp$obs_y[ixa]
    di <- array( data=envtmp$D[ixa,], dim=c(p,envtmp$k_dim))
    Xa_j  <- array( data=envtmp$Xa_j[i,],   dim=c(1,envtmp$k_dim))
    X_j1 <- array( data=envtmp$X_j1[ixa,], dim=c(p,envtmp$k_dim))
    eps2 <- envtmp$eps2[i]

    Cxx_jj1  <- tcrossprod( Xa_j, X_j1)     
    Cxx_j1j1 <- tcrossprod( X_j1, X_j1)     
    #
    CxxCqq_inv <- chol2inv(chol( ( Cxx_j1j1+diag(x=eps2,p)) ))
    CxxCqq_inv_di <- crossprod(CxxCqq_inv,di)       
    Ea <- envtmp$Ea[i,] + crossprod( t(Cxx_jj1), CxxCqq_inv_di)
    if (idi) xidi <- sum( Cxx_jj1 * as.vector(rowSums(CxxCqq_inv)))
  }
  return( c( Ea, xidi))
}

