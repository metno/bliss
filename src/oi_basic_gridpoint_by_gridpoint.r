#+ optimal interpolation 
oi_basic_gridpoint_by_gridpoint<-function( i,
                                           corr="soar",
                                           pmax,
                                           isolation_factor=7,
                                           dh=10000,
                                           idi=T,
                                           uncertainty=F) {
# returned values: analysis, IDI, observation error var, analysis error var
#------------------------------------------------------------------------------
  xa <- NA; xidi <- NA; o_errvar <- NA; xa_errvar <- NA

  if( i%%(round(envtmp$m_dim/10)) == 0) cat(".")

  # select the observations to use
  if ( (p <- length( aux <- which(envtmp$nn2$nn.idx[i,]!=0))) == 0) {

    # no observations, analysis=background
    xa <- envtmp$xb[i]
  } else {

    # observations available
    ixa  <- envtmp$nn2$nn.idx[i,aux]

    # define vectors
    dist <- envtmp$nn2$nn.dists[i,aux]
    x <- y_env$yo$x[ixa]
    y <- y_env$yo$y[ixa]
    yo <- y_env$yo$value[ixa]
    yb <- envtmp$yb[ixa]
    di <- yo - yb
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
    xa <- envtmp$xb[i] + sum( rloc * as.vector(SRinv_di))
    if (idi) xidi <- sum(rloc*as.vector(rowSums(SRinv)))
    if (uncertainty) {
      o_errvar  <- mean( di * ( di - crossprod(S,SRinv_di)))
      xa_errvar <- ( o_errvar/ eps2) * 
                   ( 1 - sum( as.vector( crossprod( rloc, SRinv)) * rloc))
    }
  }
  return( c( xa, xidi, o_errvar, xa_errvar))
}

