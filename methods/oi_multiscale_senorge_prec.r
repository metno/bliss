#+
oi_multiscale_senorge_prec <- function( argv, y_env, fg_env, env) {
#------------------------------------------------------------------------------
  t0a <- Sys.time()

  cat( "-- OI multiscale developed for seNorge --\n")
  # If there is at least one observation of precipitation (not dry everywhere)
  if (y_env$yo$nwet > 0) {

    # number of background ensemble members
    if ( ( nfg <- length( fg_env$fg)) != 1) return( FALSE)
    rrf <- subset( fg_env$fg[[1]]$r_main, subset=1)

    # number of observations
    res_master <- mean( res(env$rmaster))

    # define sequence of nearest station values
    if (y_env$yo$n < 100) {
      kseq <- rev( unique( c(2,3,4,seq(5,y_env$yo$n,by=5), y_env$yo$n)))
    } else {
      kseq <- rev( unique( c(2,3,4,seq(5,100,by=5),seq(100,y_env$yo$n,by=200), y_env$yo$n)))
    }
    nn2 <- nn2( cbind(y_env$yo$x,y_env$yo$y), 
                       query = cbind(y_env$yo$x,y_env$yo$y), 
                       k = y_env$yo$n, searchtype = "radius", 
                       radius = 100000000)
    vecd_tmp <- unique( sort( c( (res_master*round(colMeans(nn2$nn.dists)[kseq]/res_master)), (2:100)*res_master), decreasing=T))
    vecf_tmp <- pmin( res_master*min(c(env$nx,env$ny))/3, pmax( 1*res_master, round(vecd_tmp/2)))
    vecd <- unique( sort( c( vecd_tmp[which(!duplicated(vecf_tmp,fromLast=T) & vecd_tmp >= (2*res_master))],
                             vecd_tmp[which(!duplicated(vecf_tmp,fromLast=F) & vecd_tmp >= (2*res_master))]),
                             decreasing=T))
    # aggregation factor is half the horizontal decorrelation length
    vecf <- round( vecd/(2*1000), 0)
    rm( vecf_tmp, vecd_tmp)
    # same weight to the obs and the background from previous iteration
    vece<-rep(argv$eps2,length(vecf))
    nl<-length(vecd)
    vece[]<-1
    print("+---------------------------------------------------------------+")
    print("vecd vecf vece")
    print(" km   -    -")
    print(cbind(vecd,vecf,vece))
    print("+---------------------------------------------------------------+")

    # Distances between all stations
    Disth <- ( outer(y_env$yo$x,y_env$yo$x,FUN="-")**2.+
               outer(y_env$yo$x,y_env$yo$y,FUN="-")**2. )**0.5
  
    yrf <- extract( rrf, cbind(y_env$yo$x, y_env$yo$y))
    # if rescaling factor exists, multi-scale OI operates on relative anomalies 
    if (argv$use_relativeAnomalies) {
      if (any(yrf==0)) yrf[which(yrf==0)] <- 1
      yo_relan <- y_env$yo$value / yrf
    } else {
      yo_relan <- y_env$yo$value
    }
    zero <- 0
    if ( argv$transf == "Box-Cox") {
      yo_relan <- boxcox( yo_relan, argv$transf.boxcox_lambda)
      zero <- boxcox( 0, argv$transf.boxcox_lambda)
    } else if (argv$transf!="none") {
      boom("transformation not defined")
    }

    # multi-scale OI
    for (l in 1:nl) {
      cat( paste0( "scale # ",formatC(l,width=3,flag="0"),
                   " of ",nl,
                   " (",formatC(vecd[l],width=8,flag="0",format="d"),"m ",
                   "fact=",formatC(vecf[l],width=3,flag="0"),")\n"))
      D <- exp(-0.5*(Disth/vecd[l])**2.)
      # ovarc is observation error variance correction factor (read_obs)
      diag(D) <- diag(D) + vece[l] * y_env$yo$ovarc
      envtmp$InvD <- chol2inv(chol(D))
      if (l==nl | vecf[l]==1) {
        r <- env$rmaster
      } else {
        r <- aggregate(env$rmaster, fact=vecf[l], expand=T, na.rm=T)
      }
      zvalues.l <- getValues(r)
      storage.mode(zvalues.l) <- "numeric"
      xy.l <- xyFromCell(r,1:ncell(r))
#      x.l <- sort(unique(xy.l[,1]))
#      y.l <- sort(unique(xy.l[,2]),decreasing=T)
      mask.l  <- which(!is.na(zvalues.l))
#      zgrid.l <- zvalues.l[mask.l]
      xgrid.l <- xy.l[mask.l,1]
      ygrid.l <- xy.l[mask.l,2]
      rm(xy.l,zvalues.l)
      if (l==1) {
        yb <- rep( mean(yo_relan), length=y_env$yo$n)
        xb <- rep( mean(yo_relan), length=length(xgrid.l))
      } else {
        if ("scale" %in% argv$off_x.variables) 
          xl_tmp <- getValues( resample( rl, r, method="ngb"))[mask.l]
        rb <- resample( ra, r, method="bilinear")
        xb <- getValues(rb)[mask.l]
        count <- 0
        while ( any( is.na(xb))) {
          count <- count+1
          buffer_length <- round(vecd[l]/(10-(count-1)))
          if (!is.finite(buffer_length)) break 
          ib <- which(is.na(xb))
          aux <- extract( rb, cbind(xgrid.l[ib],ygrid.l[ib]),
                          na.rm=T,
                          buffer=buffer_length)
          for (ll in 1:length(aux)) xb[ib[ll]] <- mean(aux[[ll]],na.rm=T)
          rb[mask.l] <- xb
          rm(aux,ib)
        }
        yb <- extract( rb, cbind(y_env$yo$x,y_env$yo$y), method="bilinear")
        rm(rb)
      }
      if (any(is.na(xb))) print("xb is NA")
      if (any(is.na(yb))) print("yb is NA")
      xa.l <- OI_RR_fast( yo    = yo_relan,
                          yb    = yb,
                          xb    = xb,
                          xgrid = xgrid.l,
                          ygrid = ygrid.l,
                          VecX  = y_env$yo$x,
                          VecY  = y_env$yo$y,
                          Dh    = vecd[l],
                          zero  = zero) 
      ra <- r
      ra[] <- NA
      ra[mask.l] <- xa.l
      if ("scale" %in% argv$off_x.variables) {
        rl <- ra
        rl[mask.l] <- vecd[l]
        if (l>1) {
          ixl <- which( abs(xa.l-xb) < 0.005)
          if (length(ixl)>0) rl[mask.l[ixl]] <- xl_tmp[ixl]
        }
      }
      print(range(yo_relan,na.rm=T))
      print(range(yb,na.rm=T))
      print(range(xb,na.rm=T))
      print(range(xa.l,na.rm=T))
    } # end of multi-scale OI

    # compute IDI
    if ( !is.na(argv$dh_idi)) {
      D <- exp(-0.5*(Disth/argv$dh_idi)**2.)
      # ovarc is observation error variance correction factor (read_obs)
      diag(D) <- diag(D) + 0.1
      envtmp$InvD <- chol2inv(chol(D))
      xidi <- OI_RR_fast( yo    = rep(1,y_env$yo$n),
                          yb    = rep(0,y_env$yo$n),
                          xb    = rep(0,length(xgrid.l)),
                          xgrid = xgrid.l,
                          ygrid = ygrid.l,
                          VecX  = y_env$yo$x,
                          VecY  = y_env$yo$y,
                          Dh    = argv$dh_idi)
    } else {
      xidi <- rep( 0, length(xgrid.l))
    }

    # back to precipitation values (from relative anomalies)
    if ( argv$transf == "Box-Cox") {
      xa <- tboxcox( xa.l, argv$transf.boxcox_lambda)
    } else if ( argv$transf != "none") {
      xa <- xa.l
    }
    xa_rel <- xa
    if (argv$use_relativeAnomalies) {
      xa <- xa * getValues(rrf)[mask.l]
    }

    # refine analysis
    print("refine prec/no-prec borders and remove wet regions with no obs in them")
    for (frac_rr in c(100,50,10)) {
      if (frac_rr==100) xa[which(!is.na(xa) & xa<(y_env$rain/frac_rr))]<-0
      ra[mask.l]<-xa
      xa_aux<-getValues(ra)
      xa_aux[which(!is.na(xa_aux) & xa_aux<(y_env$rain/frac_rr))]<-NA
      ra[]<-xa_aux
      rclump<-clump(ra)
      oclump<-extract(rclump,cbind(y_env$yo$x[y_env$yo$ixwet],y_env$yo$y[y_env$yo$ixwet]))
      fr<-freq(rclump)
      # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
      ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
      xa[which(getValues(rclump)[mask.l] %in% fr[ix,1])]<-0
      rm(xa_aux,rclump,oclump,fr,ix)
      ra[mask.l]<-xa
    }

  # If dry everywhere
  } else {
    ra <- env$rmaster
    ra[]     <- NA
    ra[mask] <- 0
    rl <- ra
  } # END IF dry everywhere
save(file="tmp.rdata",argv, y_env, fg_env, env,ra,xa_rel,rl,xidi)
  # Initialization
  env$Xa <- array( data=NA, dim=c( env$ngrid, 1))
  env$Xa_rel <- array( data=NA, dim=c( env$ngrid, 1))
  env$Xidi <- array( data=NA, dim=c( env$ngrid, 1))
  env$Xscale <- array( data=NA, dim=c( env$ngrid, 1))
  y_env$yo$value_a <- array( data=NA, dim=c( y_env$yo$n, 1))
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a <- array( data=NA, dim=c( y_env$yov$n, 1))

  # Output  
  env$Xa[,1] <- getValues(ra)[mask.l]
  env$Xa_rel[,1] <- xa_rel
  env$Xscale[,1] <- getValues(rl)[mask.l]
  env$Xidi[,1] <- xidi
  y_env$yo$value_a[,1] <- extract( ra, cbind( y_env$yo$x, y_env$yo$y), method="bilinear")
  if (env$cv_mode | env$cv_mode_random) 
    y_env$yov$value_a[,1] <- extract( ra, cbind( y_env$yov$x, y_env$yov$y), method="simple")

  # Bye-bye
  t1a <- Sys.time()
  cat( paste( "OI multiscale developed for seNorge, total time", round(t1a-t0a,1), attr(t1a-t0a,"unit"), "\n"))

}
