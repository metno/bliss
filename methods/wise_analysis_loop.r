#+
wise_analysis_loop <- function( argv, y_env, fg_env, env,
                                supob_nobs=50, supob_radius=1500, supob_q=0.99,
                                max_it=100, opttol=0.02,
                                En2_adj_fun="Gaussian",
                                En2_adj_min=0,
                                rescale_min_obs=100, rescale_min_cells=100,
                                plot=F, dir_plot=NA) {
#
#------------------------------------------------------------------------------

  t0a <- Sys.time()

  cat( "-- wise analysis --\n")

  # number of background ensemble members
  if ( ( nfg <- length( fg_env$ixs)) == 0) return( FALSE)

  # set dyadic domain
  xmn <- as.numeric(argv$grid_master.x1) - as.numeric(argv$grid_master.resx)/2
  xmx <- as.numeric(argv$grid_master.xn) + as.numeric(argv$grid_master.resx)/2
  ymn <- as.numeric(argv$grid_master.y1) - as.numeric(argv$grid_master.resy)/2
  ymx <- as.numeric(argv$grid_master.yn) + as.numeric(argv$grid_master.resy)/2

  nx <- ncol( env$rmaster)
  ny <- nrow( env$rmaster)

  n <- ceiling( log( max(nx,ny), 2))
  rdyad <- raster( extent( xmn, xmx, ymn, ymx), 
                   ncol=2**n, nrow=2**n, crs=argv$grid_master.proj4)
  rdyad[] <- 0
  rfxb <- rdyad
  rfyb <- rdyad
  rfin <- rdyad

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

  # ---~--------------
  # observations on the dyadic grid
  rfobs <- resample( env$mergeobs$r, rdyad, method="bilinear", na.rm=T)
  rfobs[rfobs<y_env$rain] <- 0

  r <- env$rmaster
  r[] <- env$mergeobs$idi
  rfidi <- resample( r, rdyad, method="bilinear", na.rm=T)
  rfidi[rfidi<0] <- 0
  rfidi[rfidi>1] <- 1
  r[] <- env$mergeobs$rall
  rfobsall <- resample( r, rdyad, method="bilinear", na.rm=T)
  rfobsall[rfobsall<y_env$rain] <- 0

  vobs <- getValues(rfobs)
  vobs_val <- !is.na(getValues(rfobs))
  vobs_nan <- is.na(getValues(rfobs))

  c_xy <- which( !is.na( getValues(rfobs)))
  cneg_xy <- which( is.na( getValues(rfobs)))
  env$rfobs <- rfobs
  env$rfidi <- rfidi
  env$rfobsall <- rfobsall

  rfobs[cneg_xy] <- 0
  rfidi[is.na(rfidi)] <- 0
  rfobsall[is.na(rfobsall)] <- 0
  vidi <- getValues(rfidi)
  vobsall <- getValues(rfobsall)

  rfobsidi <- rfidi
#  rfobsidi[cneg_xy] <- 0
#  rfobsidi[c_xy] <- 1

  rfnoidi <- 1 - rfidi
#  rfnoidi[cneg_xy] <- 1
#  rfnoidi[c_xy] <- 0

  rfobsidiwet <- rfidi
  rfobsidiwet[] <- 0
  rfobsidiwet[which(vobsall >= y_env$rain & vidi >= 0.1)] <- 1
#  rfobsidiwet[which(vobs  < y_env$rain & vobs_val)] <- 0
#  rfobsidiwet[which(vobs >= y_env$rain & vobs_val)] <- 1
#  rfobsidiwet[cneg_xy] <- 0

  rfobsididry <- rfobs
  rfobsididry[] <- 0
  rfobsididry[which(vobsall < y_env$rain & vidi >= 0.1)] <- 1
#  rfobsididry[which(vobs  < y_env$rain & vobs_val)] <- 1
#  rfobsididry[which(vobs >= y_env$rain & vobs_val)] <- 0
#  rfobsididry[cneg_xy] <- 0
  
  # ---~--------------
  # Initialization of the wavelet structures
  #  dwt_out. dwt.2d class.
  #    for the i-th level (i=1,...,env$n_levs_mx) 
  #      dwt_out[[3*(i-1)+1]] LHi - wavelet coefficients
  #      dwt_out[[3*(i-1)+2]] HLi - wavelet coefficients
  #      dwt_out[[3*(i-1)+3]] HHi - wavelet coefficients
  #  resolution of the dyadic tree. coefficients = 2**i; base = 2**(i-1)
  #  resolution is the number of original grid points (in each direction) within the box a dydadic tree at level i
  # then i=1 is the finest resolution and i=n_levs_mx the coarser
  dwt_out <- dwt.2d( as.matrix(rdyad), wf=env$wf, J=env$n_levs_mx, boundary=env$boundary)
  for (i in 1:length(dwt_out)) dwt_out[[i]][] <- 0

  # ---~--------------
  # define constants
  # total number of coefficients used. Dimension of the state variable
  env$n_dim <- length( as.vector( unlist( dwt_out)))
  # total number of coefficients used. Dimension of the state variable
  env$m_dim <- length( getValues(rdyad))
  sqrt_m_dim <- sqrt( env$m_dim)
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( "state variable, wavelet coefficients, n dim >", env$n_dim, "\n"))
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))
  cat( paste( "number of observations, superob >", length(c_xy), "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
 
  # transformation operator
  dwtobs  <- dwt.2d( as.matrix(rfobs),    wf=env$wf, J=n, boundary=env$boundary)
  dwtidi  <- dwt.2d( as.matrix(rfobsidi), wf=env$wf, J=n, boundary=env$boundary)
  dwtnoidi   <- dwt.2d( as.matrix(rfnoidi), wf=env$wf, J=n, boundary=env$boundary)
  dwtidiwet  <- dwt.2d( as.matrix(rfobsidiwet), wf=env$wf, J=n, boundary=env$boundary)
  dwtididry  <- dwt.2d( as.matrix(rfobsididry), wf=env$wf, J=n, boundary=env$boundary)
  # multires representations
  mrobs <- list()
  mrnoidi <- list()
  mridi <- list()
  mridiwet <- list()
  mrididry <- list()
  # unlist the result and compute squared energies
  vo    <- vector( mode="numeric", length=env$n_dim); vo[]<-NA
  nn <- 0
  for (l in 1:n) {
print(l)
    ij <- ijFromLev( n, l, F)
    vo[ij[1,1]:ij[3,2]] <- c( dwtobs[[3*l-2]], dwtobs[[3*l-1]], dwtobs[[3*l]])
    mrobs[[l]] <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
#    mrnoidi[[l]] <- dwt.2d( as.matrix(rfnoidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    mridi[[l]] <- dwt.2d( as.matrix(rfobsidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    print(dim(mridi[[l]])[2])
    if ( dim(mridi[[l]])[1] > 2)
      mridi[[l]] <- gauss2dsmooth(x=mridi[[l]],lambda=1,nx=dim(mridi[[l]])[1],ny=dim(mridi[[l]])[2])
    mrnoidi[[l]] <- 1 - mridi[[l]]
    mridiwet[[l]] <- dwt.2d( as.matrix(rfobsidiwet), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    if ( dim(mridiwet[[l]])[1] > 2)
      mridiwet[[l]] <- gauss2dsmooth(x=mridiwet[[l]],lambda=1,nx=dim(mridiwet[[l]])[1],ny=dim(mridiwet[[l]])[2])
#    mrididry[[l]] <- dwt.2d( as.matrix(rfobsididry), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]/2**l
    if ( any( mridi[[l]] > mrnoidi[[l]])) nn <- l
  }
  ij <- ijFromLev( n, n, T)
  vo[ij[1]:ij[2]] <- dwtobs[[3*n+1]]

  if (nn == 0) return(NULL)
  env$nn_dim <- max( ijFromLev( n, (nn+1), F))
  env$nn1_dim <- sqrt( diff(ijFromLev( n, nn, T)) + 1 )
  cat("\n")
#
#moveDown <- function(matin) {
#  dimout <- 2*dim(matin)[1]
#  matout <- array( data=NA, dim=c(dimout,dimout))
#  for (i in 1:dimout) {
#    matout[2*i,]   <- rep( matin[i,], each=2)
#    matout[2*i-1,] <- rep( matin[i,], each=2)
#  }
#}
  #--------------------------------------------------------    
  # Analysis
  cat("Analysis on the transformed space \n")

  Xb_dyad <- array( data=NA, dim=c( env$m_dim, env$k_dim))
  Xb_dyad_original <- array( data=NA, dim=c( env$m_dim, env$k_dim))

  #--------------------------------------------------------    
  # MAIN LOOP
  for (loop in 1:max_it) {

    t0 <- Sys.time()
    cat( paste(" loop", loop, " "))

    #--------------------------------------------------------    
    # set the background

    Ub     <- array( data=NA, dim=c( env$nn_dim, env$k_dim))
    Vb     <- array( data=NA, dim=c( env$nn_dim, env$k_dim))
    vo_Vb  <- array( data=NA, dim=c( env$nn_dim, env$k_dim))

    if (plot & loop == 1) fxb <- vector()
    
    for (e in 1:env$k_dim) {

      # first iteration - background from ensemble 
      if ( loop == 1) {
        i <- fg_env$ixs[e]
        # interpolate onto the dyadic grid
        # background at grid  points
        rfxb <- resample( subset( fg_env$fg[[fg_env$ixf[i]]]$r_main, subset=fg_env$ixe[i]), rdyad, method="bilinear")
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad_original[,e]  <- getValues(rfxb)
        Xb_dyad[,e]  <- getValues(rfxb)
        env$Xa_dyad  <- Xb_dyad
        env$Xa_dyad[] <- NA
      # second iteration onwards - background from previous iteration 
      } else {
        # -- background at grid  points --
        rfxb[] <- env$Xa_dyad[,e]
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad[,e]  <- getValues(rfxb)
      }

      # background at observation points
      rfyb <- rdyad; rfyb[c_xy] <- rfxb[c_xy]
      # innovation 
      rfin <- rdyad; rfin[c_xy] <- rfobs[c_xy] - rfxb[c_xy]

      # transformation operator
      dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=(nn+1), boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:(nn+1)) {
        ij <- ijFromLev( n, l, F)
        Ub[ij[1,1]:ij[3,2],e] <- c( dwtxb[[3*l-2]], dwtxb[[3*l-1]], dwtxb[[3*l]])
        Vb[ij[1,1]:ij[3,2],e] <- c( dwtyb[[3*l-2]], dwtyb[[3*l-1]], dwtyb[[3*l]])
        vo_Vb[ij[1,1]:ij[3,2],e] <- c( dwtin[[3*l-2]], dwtin[[3*l-1]], dwtin[[3*l]])
      }
    } # end loop over ensemble members

    # set the background: end
    #--------------------------------------------------------    

    #--------------------------------------------------------    
    # Background error covariance matrices
    Ab <- Ub - rowMeans(Ub)
    Db <- Vb - rowMeans(Vb)
    var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) 
    var_d1_yy <- 1/(env$k_dim-1) * rowSums( vo_Vb *vo_Vb) 
    # scaling coefficient used in spatial analysis
    coeff <- var_b1_xy / var_d1_yy
    coeff[!is.finite(coeff)] <- 0
    rm( Ab, Db)

#    #--------------------------------------------------------    
#    # sort of root mean squared error of the decomposed innovation
#    env$costf[loop] <- mean( sqrt( var_d1_yy))
#    cat( paste("rmse", round(env$costf[loop],5)), "\n")

    #--------------------------------------------------------    
    # Rescale the different components
    
    Ua  <- array( data=NA, dim=c(     env$nn_dim, env$k_dim))
    for (e in 1:env$k_dim) {
      cat(".")
      Ua[,e] <- Ub[,e] + coeff * vo_Vb[,e]
      rfxb[] <- Xb_dyad[,e]
      r <- rdyad
      for (l in (nn+1):1) {
        dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=l, boundary=env$boundary)
        for (i in 1:(length(dwt)-1)) dwt[[i]][] <- 0
        if (l < (nn+1)) dwt[[3*l+1]][] <- mrxa
        ij <- ijFromLev( n, l, F)
        dwt[[3*l-2]][] <- Ua[ij[1,1]:ij[1,2],e]
        dwt[[3*l-1]][] <- Ua[ij[2,1]:ij[2,2],e]
        dwt[[3*l]][]   <- Ua[ij[3,1]:ij[3,2],e]
        xbadj <- idwt.2d( dwt) 
        r[] <- array(data=as.matrix(xbadj),dim=c(sqrt_m_dim,sqrt_m_dim))
        if (l > 1) {
          dwtxbadj <- dwt.2d( as.matrix(r), wf=env$wf, J=(l-1), boundary=env$boundary)
          mrxbadj <- dwtxbadj[[3*(l-1)+1]]
#          mrxbadj[mrididry[[l-1]]>mridiwet[[l-1]]] <- 0
          mrxbadj <- mrxbadj * mridiwet[[l-1]]
          mrxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*(l-1)+1]]
          mrxa <- mrnoidi[[l-1]] * mrxb + mridi[[l-1]] * mrxbadj
          if (!is.na(y_env$rain)) mrxa[mrxa<y_env$rain] <- 0
        } else {
          xbadj <- getValues(r)
#          xbadj[getValues(rfobsididry)>getValues(rfobsidiwet)] <- 0
          xbadj <- getValues(rfobsidiwet) * xbadj
          env$Xa_dyad[,e] <- getValues(rfnoidi) * Xb_dyad[,e] + getValues(rfobsidi) * xbadj
          if (!is.na(y_env$rain)) env$Xa_dyad[,e][env$Xa_dyad[,e]<y_env$rain] <- 0 
        }
      }  # end loop over scales
    }  # end loop over ensembles

    #
    #--------------------------------------------------------    
    # strict condition for convergence
    r[] <-  rowMeans(Xb_dyad)
    aux <- sqrt( mean( (extract(r,cbind(y_env$yo$x,y_env$yo$y)) - y_env$yo$value))**2)
    r[] <-  rowMeans(env$Xa_dyad)
    env$costf[loop] <- sqrt( mean( (extract(r,cbind(y_env$yo$x,y_env$yo$y)) - y_env$yo$value))**2)
    cat(paste("\n","rmse before after Delta% =",round(aux,5),round(env$costf[loop],5),
                                             round(100*(aux-env$costf[loop])/aux,1),"\n"))
    # less strict condition for convergence
    t1 <- Sys.time()
#    cat( paste( "time=", round(t1-t0,1), attr(t1-t0,"unit"), "costf - Delta En2 index =", round(env$costf[loop],5), "\n"))
    # break out of the main loop early if variations  
#    if ( env$costf[loop] < opttol) break
if (loop==20) {
save(file="tmp.rdata",coeff,env,mrobs,mrnoidi,mridi,mridiwet,mrididry,mrxa,n,nn,rdyad,Xb_dyad,Ua,Ub,vo_Vb,y_env,mrxb,mrxbadj,Xb_dyad_original,rfobs,rfnoidi,rfobsidi,rfobsidiwet,rfobsididry)
#i<-1; s<-rdyad; s[]<-env$Xa_dyad[,i]; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#i<-1; s<-rdyad; s[]<-Xb_dyad_original[,i]; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#points(y_env$yo$x,y_env$yo$y,cex=0.25,pch=21,bg="beige",col="beige")
# i<-1; s[]<-Xb_dyad_original[,i]; yb<-extract(s,cbind(y_env$yo$x,y_env$yo$y)); plot(y_env$yo$value,yb)
# i<-1; s[]<-env$Xa_dyad[,i]; ya<-extract(s,cbind(y_env$yo$x,y_env$yo$y)); plot(y_env$yo$value,ya)
#s<-env$rfobs; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
# 
q()
}
  } # end main loop
q()

#save(file="tmp.rdata",env,Xb_dyad,sqrt_m_dim,rfxb,c_xy)
  cat( paste( "Add the background data where no observations are available\n"))
  for (e in 1:env$k_dim) {
    cat(".")
    rfxb[] <- Xb_dyad[,e]
    bclump <- clump(rfxb)
    oclump <- bclump[c_xy]
    fro <- as.data.frame( table( oclump), stringsAsFactors=F) 
    if ( length( ixo <- which( !is.na(fro[,1]) & fro[,2] > 100)) > 0) {
      oclump_ok <- as.integer(fro[ixo,1])
    } else {
      oclump_ok <- integer(0)
    }
    Xb_dyad[which(getValues(bclump) %in% oclump_ok),e] <- 0
    rfxb[] <- Xb_dyad[,e]
    env$Xa_dyad[,e] <- env$Xa_dyad[,e] + getValues(t(rfxb))
  }
  cat( paste( "\n", "End of analysis loop time=", round(t1-t0,1), attr(t1-t0,"unit"), "\n"))
}
