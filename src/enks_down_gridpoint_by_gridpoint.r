#+ ensemble optimal interpolation 
enks_down_gridpoint_by_gridpoint<-function( i,
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
    Ea <- envtmp$Ea[i,]
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
#print("===============")
#print("envtmp$Xa_j[i,]")
#print(envtmp$Xa_j[i,])
#print("envtmp$Xa_j1[ixa,]")
#print(envtmp$Xa_j1[ixa,])
    rloc <- tcrossprod( envtmp$Xa_j[i,], envtmp$Xa_j1[ixa,])     
    S    <- tcrossprod( envtmp$Xa_j1[ixa,], envtmp$Xa_j1[ixa,])     
#print("rloc")
#print(rloc)
#print("S")
#print(S)
#print("di")
#print(di)

    #
    SRinv <- chol2inv(chol( (S+diag(x=eps2,p)) ))
#print("SRinv")
#print(SRinv)
    SRinv_di <- crossprod(SRinv,di)       
    if ( any( class(envtmp$D) == "matrix")) {
      Ea <- envtmp$Ea[i,] + crossprod( t(rloc), SRinv_di)
    } else {
      Ea <- envtmp$Ea[i,] + sum( rloc * as.vector(SRinv_di))
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

