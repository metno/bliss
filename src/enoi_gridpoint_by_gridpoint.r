#+ ensemble optimal interpolation 
enoi_gridpoint_by_gridpoint<-function( i,
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
    if ( any( class(envtmp$D) == "matrix")) {
      di <- envtmp$D[ixa,]
    } else {
      di <- envtmp$D[ixa]
    }
    eps2 <- envtmp$eps2[i]

    # correlations
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
      S<-exp(-0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
    } else if (corr=="soar")  {
      distnorm<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2) / dh 
      S<-(1+distnorm)*exp(-distnorm)
      rm(distnorm)
    } else if (corr=="powerlaw")  {
      S<-1 / (1 + 0.5*(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)/(dh*dh))
    } else if (corr=="toar")  {
      dist<-sqrt(outer(y,y,FUN="-")**2. + outer(x,x,FUN="-")**2)
      S<- (1 + dist/dh + (dist*dist)/(3*dh*dh)) * exp(-dist/dh)
      rm(dist)
    }
    #
    SRinv <- chol2inv(chol( (S+diag(x=eps2,p)) ))
    SRinv_di <- crossprod(SRinv,di)       
    if ( any( class(envtmp$D) == "matrix")) {
      Ea <- envtmp$Eb[i,] + crossprod( rloc, SRinv_di)
    } else {
      Ea <- envtmp$Eb[i,] + sum( rloc * as.vector(SRinv_di))
    }
    if (idi) xidi <- sum( rloc * as.vector(rowSums(SRinv)))
#    if (uncertainty) {
#      o_errvar  <- mean( di * ( di - crossprod(S,SRinv_di)))
#      xa_errvar <- ( o_errvar/ eps2) * 
#                   ( 1 - sum( as.vector( crossprod( rloc, SRinv)) * rloc))
#    }
  }
  return( c( Ea, xidi))
}

