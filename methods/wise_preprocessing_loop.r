#+
wise_preprocessing_loop <- function( argv, y_env, fg_env, env,
                                     max_it=100, plot=F, dir_plot=NA) {
# Initialization of the wavelet structures
#  dwt. dwt.2d class.
#    for the i-th level (i=1,...,env$n_levs_mx) 
#      dwt_out[[3*(i-1)+1]] LHi - wavelet coefficients
#      dwt_out[[3*(i-1)+2]] HLi - wavelet coefficients
#      dwt_out[[3*(i-1)+3]] HHi - wavelet coefficients
#  resolution of the dyadic tree. coefficients = 2**i; base = 2**(i-1)
#  resolution is the number of original grid points (in each direction) within the box a dydadic tree at level i
# then i=1 is the finest resolution and i=n_levs_mx the coarser
#
#------------------------------------------------------------------------------

  t0a <- Sys.time()

  cat( "-- wise preprocessing loop --\n")

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

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

  # ---~--------------
  # observations on the dyadic grid
  # rfobs covers obs-dense areas only
  rfobs <- resample( env$mergeobs$r, rdyad, method="bilinear", na.rm=T)
  rfobs[rfobs<y_env$rain] <- 0

  # rfidi covers wider areas than rfobs (=1 obs-dense; =0 obs-void)
  r <- env$rmaster
  r[] <- env$mergeobs$idi
  rfidi <- resample( r, rdyad, method="bilinear", na.rm=T)
  rfidi[rfidi<0] <- 0
  rfidi[rfidi>1] <- 1

  # rfobsall covers wider areas than rfobs
  r[] <- env$mergeobs$rall
  rfobsall <- resample( r, rdyad, method="bilinear", na.rm=T)
  rfobsall[rfobsall<y_env$rain] <- 0

  # helper
  c_xy <- which( !is.na( getValues(rfobs)))
  # save for future use
  env$rfobs <- rfobs
  env$rfidi <- rfidi
  env$rfobsall <- rfobsall

  # paddle with 0s
  rfobs[is.na(rfobs)] <- 0
  rfidi[is.na(rfidi)] <- 0
  rfobsall[is.na(rfobsall)] <- 0

  vidi <- getValues(rfidi)
  vobsall <- getValues(rfobsall)

  # =1 obs-void; =0 obs_dense
  rfnoidi <- 1 - rfidi

  # wet gridpoints in observed regions
  rfidiwet <- rfidi
  rfidiwet[] <- 0
  rfidiwet[which(vobsall >= y_env$rain & vidi >= 0.1)] <- 1

  # dry gridpoints in observed regions
  rfididry <- rfobs
  rfididry[] <- 0
  rfididry[which(vobsall < y_env$rain & vidi >= 0.1)] <- 1
  
  # ---~--------------
  # define constants
  # number of grid points on the dyadic grid
  env$m_dim <- length( getValues(rdyad))
  sqrt_m_dim <- sqrt( env$m_dim)
  # max number of wavelet coefficients is the same as the number of grid points 
  env$n_dim <- env$m_dim
  # number of observations
  env$p_dim <- length(y_env$yo$x)
  cat( paste( " number of grid points, m dim >", env$m_dim, "\n"))
  cat( paste( "number of observations, p dim >", env$p_dim, "\n"))
  cat( paste( "number of observations, superob >", length(c_xy), "\n"))

  # ---~--------------
  # -- Observations --

  cat("Transform observations and compute error covariance matrix ")
 
  # multires representations
  mrobs <- list()
  mrnoidi <- list()
  mridi <- list()
  mridiwet <- list()
  mrididry <- list()
  mrnorain<- list()
  mrw<- list()
  mrwobs<- list()
  mrwbkg<- list()
  mrwwet<- list()
  mrwdry<- list()
  # aux
  rf1<-rdyad
  rf1[]<-1
  lambda <- 1
  env$wf<-"haar"
  
  # unlist the result and compute squared energies
  nn <- 0
  for (l in 1:n) {
    aux1 <- dwt.2d( as.matrix(rf1), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]
    mridi[[l]] <- dwt.2d( as.matrix(rfidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridi[[l]])[1] > 2)
      mridi[[l]] <- gauss2dsmooth( x=mridi[[l]], lambda=lambda, nx=dim(mridi[[l]])[1], ny=dim(mridi[[l]])[2])
    if ( any( mridi[[l]] > 1)) mridi[[l]][ mridi[[l]]>1] <- 1
    if ( any( mridi[[l]] < 0)) mridi[[l]][ mridi[[l]]<0] <- 0

    mrnoidi[[l]] <- dwt.2d( as.matrix(rfnoidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mrnoidi[[l]])[1] > 2)
      mrnoidi[[l]] <- gauss2dsmooth( x=mrnoidi[[l]], lambda=lambda, nx=dim(mrnoidi[[l]])[1], ny=dim(mrnoidi[[l]])[2])
    if ( any( mrnoidi[[l]] > 1)) mrnoidi[[l]][ mrnoidi[[l]]>1] <- 1
    if ( any( mrnoidi[[l]] < 0)) mrnoidi[[l]][ mrnoidi[[l]]<0] <- 0

    if ( any( mridi[[l]] > mrnoidi[[l]])) nn <- l
  }

#  lambda <- c( seq( 2**n/2, 1, length=nn), 1)
#  lambda <- c( seq( 2**nn/2, 1, length=nn), 1)
  lambda <- c( nn/(nn-1) * (2**nn/4-1) * 1/1:nn + (nn-2**nn/4)/(nn-1), 1)
  for (l in 1:(nn+1)) {
    aux1 <- dwt.2d( as.matrix(rf1), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]
    mrobs[[l]] <- dwt.2d( as.matrix(rfobs), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]]
#    if ( any( mrobs[[l]] < mrnorain[[l]])) mrobs[[l]][mrobs[[l]]<mrnorain[[l]]] <- 0
    if ( any( mrobs[[l]] < 0)) mrobs[[l]][mrobs[[l]]<0] <- 0

    mridiwet[[l]] <- dwt.2d( as.matrix(rfidiwet), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridiwet[[l]])[1] > 2)
      mridiwet[[l]] <- gauss2dsmooth(x=mridiwet[[l]],lambda=lambda[l],nx=dim(mridiwet[[l]])[1],ny=dim(mridiwet[[l]])[2])
    if ( any( mridiwet[[l]] > 1)) mridiwet[[l]][mridiwet[[l]]>1] <- 1
    if ( any( mridiwet[[l]] < 0)) mridiwet[[l]][mridiwet[[l]]<0] <- 0
    mrididry[[l]] <- dwt.2d( as.matrix(rfididry), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mrididry[[l]])[1] > 2)
      mrididry[[l]] <- gauss2dsmooth(x=mrididry[[l]],lambda=lambda[l],nx=dim(mrididry[[l]])[1],ny=dim(mrididry[[l]])[2])
    if ( any( mrididry[[l]] > 1)) mrididry[[l]][mrididry[[l]]>1] <- 1
    if ( any( mrididry[[l]] < 0)) mrididry[[l]][mrididry[[l]]<0] <- 0
#    res <- set_weights(cbind(as.vector(mridiwet[[l]]),as.vector(mrididry[[l]])))
#    mrwwet[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    mrwdry[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    res <- set_weights(cbind(as.vector(mrididry[[l]]),as.vector(mridiwet[[l]])))
    mrwwet[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    mrwdry[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    if ( any( mrwdry[[l]] > mrwwet[[l]])) mrobs[[l]][mrwdry[[l]] > mrwwet[[l]]] <- 0

    mridi[[l]] <- dwt.2d( as.matrix(rfidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mridi[[l]])[1] > 2)
      mridi[[l]] <- gauss2dsmooth( x=mridi[[l]], lambda=lambda[l], nx=dim(mridi[[l]])[1], ny=dim(mridi[[l]])[2])
    if ( any( mridi[[l]] > 1)) mridi[[l]][mridi[[l]]>1] <- 1
    if ( any( mridi[[l]] < 0)) mridi[[l]][mridi[[l]]<0] <- 0

    mrnoidi[[l]] <- dwt.2d( as.matrix(rfnoidi), wf=env$wf, J=l, boundary=env$boundary)[[3*l+1]] / aux1
    if ( dim(mrnoidi[[l]])[1] > 2)
      mrnoidi[[l]] <- gauss2dsmooth( x=mrnoidi[[l]], lambda=lambda[l], nx=dim(mrnoidi[[l]])[1], ny=dim(mrnoidi[[l]])[2])
    if ( any( mrnoidi[[l]] > 1)) mrnoidi[[l]][ mrnoidi[[l]]>1] <- 1
    if ( any( mrnoidi[[l]] < 0)) mrnoidi[[l]][ mrnoidi[[l]]<0] <- 0

    res <- set_weights(cbind(as.vector(mridi[[l]]),as.vector(mrnoidi[[l]])))
    mrwobs[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
    mrwbkg[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    res <- set_weights(cbind(as.vector(mrnoidi[[l]]),as.vector(mridi[[l]])))
#    mrwobs[[l]] <- array(data=as.matrix(res[,2]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
#    mrwbkg[[l]] <- array(data=as.matrix(res[,1]),dim=c(sqrt_m_dim/2**l,sqrt_m_dim/2**l))
  }
  
save(file="tmp.rdata",lambda,mridi,mrnoidi,mridiwet,mrididry,mrobs,mrnorain,rdyad,mrw,mrwobs,mrwbkg,mrwwet,mrwdry)
print(nn)
q()
#i<-1;r<-aggregate(rdyad,fact=2**i);r[]<-mrobs[[i]];image(r,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#i<-1;r<-aggregate(rdyad,fact=2**i);r[]<-mrwobs[[i]];contour(r,levels=c(0,0.4,0.8),add=T)


  if (nn == 0) return(NULL)
  env$nn_dim <- max( ijFromLev( n, (nn+1), F))
  cat("\n")

  #--------------------------------------------------------    
  # Pre-processing
  cat("Pre-processing on the transformed space \n")

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
        Xbpp_dyad  <- Xb_dyad
        Xbpp_dyad[] <- NA
      # second iteration onwards - background from previous iteration 
      } else {
        # -- background at grid  points --
        rfxb[] <- Xbpp_dyad[,e]
        rfxb[rfxb<y_env$rain] <- 0
        Xb_dyad[,e]  <- getValues(rfxb)
      }

      # background at observation points
      rfyb <- rdyad; rfyb[c_xy] <- rfxb[c_xy]
      # innovation 
      rfin <- rdyad; rfin[c_xy] <- rfobs[c_xy] - rfxb[c_xy]
print(rfxb)
print(rfyb)
print(rfobs)
print(rfin)
      # transformation operator
      dwtxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtyb <- dwt.2d( as.matrix(rfyb), wf=env$wf, J=(nn+1), boundary=env$boundary)
      dwtin <- dwt.2d( as.matrix(rfin), wf=env$wf, J=(nn+1), boundary=env$boundary)

      # unlist the result and compute squared energies
      jj <- 0; ii <- 0
      for (l in 1:(nn+1)) {
        ij <- ijFromLev( n, l, F)
print("---------")
print(l)
print(ij)
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
    Ao_b <- vo_Vb - rowMeans(vo_Vb)
    var_b1_xy <- 1/(env$k_dim-1) * rowSums( Ab*Db) 
#    var_d1_yy <- 1/(env$k_dim-1) * rowSums( vo_Vb *vo_Vb) 
    var_d1_yy <- 1/(env$k_dim-1) * rowSums( Ao_b * Ao_b) 
    # scaling coefficient used in pre-processing
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
# revise the weighting function as in set_weights
      for (l in (nn+1):1) {
        dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=l, boundary=env$boundary)
        for (i in 1:(length(dwt)-1)) dwt[[i]][] <- 0
        if (l < (nn+1)) dwt[[length(dwt)]][] <- mrxa
        ij <- ijFromLev( n, l, F)
        dwt[[3*l-2]][] <- Ua[ij[1,1]:ij[1,2],e]
        dwt[[3*l-1]][] <- Ua[ij[2,1]:ij[2,2],e]
        dwt[[3*l]][]   <- Ua[ij[3,1]:ij[3,2],e]
        xbadj <- idwt.2d( dwt) 
        if (!is.na(y_env$rain)) xbadj[xbadj<y_env$rain] <- 0 
        r[] <- array(data=as.matrix(xbadj),dim=c(sqrt_m_dim,sqrt_m_dim))
        if (l > 1) {
          mrxbadj <- dwt.2d( as.matrix(r), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*(l-1)+1]]
          mrxb <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(l-1), boundary=env$boundary)[[3*(l-1)+1]]
          mrxbadj <- mrwbkg[[l-1]] * mrxb + mrwobs[[l-1]] * mrxbadj
#          if (!is.na(y_env$rain)) { mrxb[mrxb<y_env$rain] <- 0; mrxbadj[mrxbadj<y_env$rain] <- 0 }
          mrxa <- mrwbkg[[l-1]] * mrxbadj + mrwobs[[l-1]] * mrobs
        } else {
          xbadj <- getValues(r)
          xbadj <- getValues(rfidiwet)/(getValues(rfidiwet)+getValues(rfididry)) * rfobs +  getValues(rfididry)/(getValues(rfidiwet)+getValues(rfididry)) * xbadj
          xbadj[getValues(rfididry)>getValues(rfidiwet)] <- 0
###          xbadj <- getValues(rfidiwet) * xbadj
          if (!is.na(y_env$rain)) { xbadj[xbadj<y_env$rain] <- 0 }
          Xbpp_dyad[,e] <- getValues(rfnoidi)/(getValues(rfnoidi)+getValues(rfidi)) * Xb_dyad[,e] + getValues(rfidi)/(getValues(rfnoidi)+getValues(rfidi)) * xbadj
          if (!is.na(y_env$rain)) Xbpp_dyad[,e][Xbpp_dyad[,e]<y_env$rain] <- 0 
        }

###        #if (e==1) {
###        if (loop==3) {
###          if (l > 1) {
###            dwt <- dwt.2d( as.matrix(rfxb), wf=env$wf, J=(l-1), boundary=env$boundary)
###            for (i in 1:(length(dwt)-1)) dwt[[i]][] <- 0
###            dwt[[length(dwt)]][] <- mrxa
###          #  s<-aggregate(rdyad,fact=2**(l-1))
###            s<-rdyad
###            sb<-rdyad
###            sbadj<-rdyad
###            sobs<-rdyad
###            s[]<-array(data=as.matrix(idwt.2d( dwt)),dim=c(sqrt_m_dim,sqrt_m_dim))
###            dwt[[length(dwt)]][] <- mrxbadj
###            sbadj[]<-array(data=as.matrix(idwt.2d( dwt)),dim=c(sqrt_m_dim,sqrt_m_dim))
###            dwt[[length(dwt)]][] <- mrobs[[l-1]]
###            sobs[]<-array(data=as.matrix(idwt.2d( dwt)),dim=c(sqrt_m_dim,sqrt_m_dim))
###            dwt[[length(dwt)]][] <- mrxb
###            sb[]<-array(data=as.matrix(idwt.2d( dwt)),dim=c(sqrt_m_dim,sqrt_m_dim))
###          } else {
###            s<-rdyad
###            sb<-rdyad
###            sbadj<-rdyad
###            s[]<-Xbpp_dyad[,e]
###            sbadj[]<-xbadj
###            sb[]<-Xb_dyad[,e]
###            sobs<-env$rfobs
###          }
###          ffout <- paste0( "pngs/fig_",
###                           formatC(loop,width=2,flag="0"),"_",
###                           formatC(e,width=2,flag="0"),"_",
###                           formatC(l,width=2,flag="0"), ".png")
###          ffout1<-paste0(ffout,".1")
###          ffout2<-paste0(ffout,".2")
###          ffout3<-paste0(ffout,".3")
###          ffout4<-paste0(ffout,".4")
###          png(file=ffout1,width=1200,height=1200)
###          image(sobs,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
###          image(s,breaks=c(-1000,0.1,1000),col=c("white","black"),add=T)
###          image(sobs,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))),add=T)
###          aux<-dev.off()
###          png(file=ffout2,width=1200,height=1200)
###          image(sbadj,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
###          aux<-dev.off()
###          png(file=ffout3,width=1200,height=1200)
###          image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
###          aux<-dev.off()
###          png(file=ffout4,width=1200,height=1200)
###          image(sb,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
###          aux<-dev.off()
###          system(paste("convert +append",ffout1,ffout3,"top.png"))
###          system(paste("convert +append",ffout4,ffout2,"bot.png"))
###          system(paste("convert -append top.png bot.png",ffout))
###          system(paste("rm top.png bot.png",ffout1,ffout2,ffout3,ffout4))
###          print(paste("written file",ffout))
###        }

      }  # end loop over scales
    }  # end loop over ensembles

    #
    #--------------------------------------------------------    
    # strict condition for convergence
    r[] <- rowMeans(Xb_dyad)
    vxb <- extract( r, cbind( y_env$yo$x, y_env$yo$y))
    r[] <- rowMeans(Xbpp_dyad)
    vxbpp <- extract( r, cbind( y_env$yo$x, y_env$yo$y))

    ets_thr <- y_env$rain
    if (any(y_env$yo$value>y_env$rain)) ets_thr <- median(y_env$yo$value[y_env$yo$value>y_env$rain])

    ets_b <- wise_ets( vxb, y_env$yo$value, ets_thr)
    ets_bpp <- wise_ets( vxbpp, y_env$yo$value, ets_thr)
    env$costf[loop] <- ets_bpp - ets_b
    cat(paste("\n","ets before after Delta% =",round(ets_b,3),round(ets_bpp,3),
                                             round(100*(ets_b-ets_bpp)/ets_b,1),"\n"))
    # break out of the main loop 
    if ( env$costf[loop] < 0) {
      cat( paste( "break out of the loop and keep the adjustments of loop", (loop-1), "\n"))
      Xbpp_dyad <- Xb_dyad
      break
    }
    t1 <- Sys.time()
#    if (loop==20) {
#      save(file="tmp.rdata",coeff,env,mrobs,mrnoidi,mridi,mridiwet,mrididry,mrxa,n,nn,rdyad,Xb_dyad,Ua,Ub,vo_Vb,y_env,mrxb,mrxbadj,Xb_dyad_original,rfobs,rfnoidi,rfidi,rfidiwet,rfididry)
#i<-1; s<-rdyad; s[]<-Xbpp_dyad[,i]; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#i<-1; s<-rdyad; s[]<-Xb_dyad_original[,i]; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
#points(y_env$yo$x,y_env$yo$y,cex=0.25,pch=21,bg="beige",col="beige")
# i<-1; s[]<-Xb_dyad_original[,i]; yb<-extract(s,cbind(y_env$yo$x,y_env$yo$y)); plot(y_env$yo$value,yb)
# i<-1; s[]<-Xbpp_dyad[,i]; ya<-extract(s,cbind(y_env$yo$x,y_env$yo$y)); plot(y_env$yo$value,ya)
#s<-env$rfobs; image(s,breaks=c(0,0.1,1,2,4,8,16,32,64,128),col=c("gray",rev(rainbow(8))))
# 
#      q()
#    }
  } # end main loop

  # save the adjusted background fields in "env"
  env$Xbpp_dyad <- Xbpp_dyad

} # END FUNCTION

wise_ets <- function(pred,ref,thr) {
  flag <- !is.na(pred) & !is.na(ref)
  hits <- length( which( flag & ref >= thr & pred >= thr))
  fals <- length( which( flag & ref <  thr & pred >= thr))
  miss <- length( which( flag & ref >= thr & pred <  thr))
  corn <- length( which( flag & ref <  thr & pred <  thr))
  hits_random <- (hits+miss) * (hits+fals) / (hits+fals+miss+corn)
  den <- (hits+fals+miss-hits_random)
  if (den == 0) return(NA)
  return( (hits-hits_random)/den )
}

#+ set IDI-based weights for the blending of multiple fields
set_weights <- function( idis, 
                         small=0.0000001) {
#------------------------------------------------------------------------------
# We have n fields to blend toghether in one field. Different fields are the
# results of spatial analysis with the same statistical interpolation (SI)
# method, applied with different settings.
# SI produces a field where each value is an average over a spatial support, 
# and the sizes of those spatial supports vary across the domain.
# Some settings are such that -on average- the spatial support is
# smaller than for others. 
# For each point in each field, an indicator -named IDI (0-1)- is available 
# that summarizes the amount of information derived from nearby observations
# which is used in the SI:
#  - 0= no info from nearby obs= spatial support large, compared to SI setup
#  - 1= a lot of info from nearby obs= spatial support comparable to SI setup
#
# The blending is done such that more weith is given to the field with smaller
# spatial support (higher resolution) and weight has been given to the other 
# (coarser) fields only if they provides additional information.
# The relative content of information in a field wrt another field is 
# computed through a function inspired by the Shannon's measure of information
# content (Tarantola, sec. 1.2.5).
# The fields are ordered from the one with smaller spatial supports to the one
# with larger spatial supports. Then, relative content of information in one
# field wrt its immediate neighbour with smaller spatial supports is computed.
# The first field, the one with smaller spatial support, is assigned the 
# completmentary of the relative content of information in the second field.
#
# The weights are the normalized relative contents of information.
#
# - Inputs -
# idis, matrix (npoints x nscales) IDI, finer scale [,1] to coarser scale [,n]
# small, constant
# - Output -
# w, matrix (npoints x nscales) weights for the blending
#------------------------------------------------------------------------------
  npoints <- dim(idis)[1]
  nscales <- dim(idis)[2]
  if ( nscales < 2) return(NULL)
  # initialization
  idis[idis<=0] <- small
  idis[idis>=1] <- 1-small
  #
  w <- array( data=NA, dim=dim(idis))
  # relative content of information in idi(i) wrt idi(i-1), coarser wrt finer
  for (i in 2:nscales) 
    w[,i] <- idis[,i] * log( idis[,i] / idis[,(i-1)])
  # adjust small things
  w[w<small] <- 0
  w[w>(1-small)] <- 1
  # weight for the finer resolution
  w[,1] <- 1 - w[,2]
  # normalize weights to 1
  w <- w / rowSums( w)
  w
}
#w <- set_weights( cbind( getValues(i200), getValues(i400), 
#                         getValues(i600), getValues(i800)))
#r <- w[,1] * r200 + w[,2] * r400 + w[,3] * r600 + w[,4] * r800


