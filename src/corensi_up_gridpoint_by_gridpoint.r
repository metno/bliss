#+ Change-Of-Resolution Ensemble Rauch-Tung-Striebel smoother UpSweep (fine-to-coarse) 
corensi_up_gridpoint_by_gridpoint<-function( i,
                                             corr = "soar",
                                             dh = 10000,
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
    Ea <- envtmp$E[i,]
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

    # static correlations
    if (corr=="gaussian") {
      rloc <- exp( -0.5* (dist*dist) / (dh*dh) )
    } else if (corr=="soar")  {
      rloc <- (1+dist/dh)*exp(-dist/dh)
    } else if (corr=="powerlaw")  {
      rloc <- 1 / (1 + 0.5*(dist*dist)/(dh*dh))
    } else if (corr=="toar")  {
      rloc <- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
    }

    if (corr=="gaussian") {
      Cyy_s<-exp(-0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
    } else if (corr=="soar")  {
      distnorm<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2) / dh 
      Cyy_s<-(1+distnorm)*exp(-distnorm)
      rm(distnorm)
    } else if (corr=="powerlaw")  {
      Cyy_s<-1 / (1 + 0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
    } else if (corr=="toar")  {
      dist<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)
      Cyy_s<- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
      rm(dist)
    }

    # combine static and dynamic correlations
    Cxy <- alpha * rloc + (1-alpha) * rloc * tcrossprod( Zi, Yi)
    Cyy <- alpha * Cyy_s + (1-alpha) * Cyy_s * tcrossprod( Yi, Yi)

    #
    CyyCdd_inv <- chol2inv( chol( (Cyy+diag(x=eps2,p)) ))
    CyyCdd_inv_di <- crossprod( CyyCdd_inv, di)       
    Ea <- envtmp$E[i,] + crossprod( t(Cxy), CyyCdd_inv_di)

    if (idi) xidi <- sum( Cxy * as.vector(rowSums(CyyCdd_inv)))
  }
  return( c( Ea, xidi))
}

