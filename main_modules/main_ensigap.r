# EnSI with Gaussian Anamorphosis for Precipitation
  # NOTE: apparently I feel better with changing all the variable names from time to time...
  # prepare stuff 
  small_const<-10**(-10)
  # backup of the original data
  Xb_bak <- xb
  yo_bak <- yo
  #  init grid stuff (aix index to not NAs background points)
  if ( ( ngrid <- length( aix)) == 0 ) { boom( "found no valid Xb grid points!") }
  grid_x <- xgrid[aix]
  grid_y <- ygrid[aix]
  # xb is a matrix (aix,nens) and I like to call it Xb now, xb become the ensemble mean
  Xb <- xb
  #  observation stuff
  # obsop_precip returns yb, ix (index over x/ycoord where all the ensemble members are finite and not NAs)
  res <- obsop_precip( r=rmaster, xcoord=VecX, ycoord=VecY)
  if ( ( nobs <- length( res$ix)) == 0 ) { boom( "found no valid Yb grid points!") }
  Yb     <- res$yb # this has dimension length(ix_obs)
  Yb_bak <- Yb
  ix_obs <- res$ix # used also for creating the output
  obs_x  <- VecX[ix_obs]
  obs_y  <- VecY[ix_obs]
  yo     <- yo[ix_obs]
  rm(res)
  # filter out unplausible values
  Xb[Xb<0] <- 0
  Yb[Yb<0] <- 0
  yo[yo<0] <- 0
  # -~- data transformation, Gaussian anamorphosis -~-
  # determine the gamma shape and rate
  #  best fitting of the empirical cdf
  #  gaussian anamorphosis for each ensemble member separately
  # Want to check wheter the estimation works?
  # plot(quantile(ecdf(yo_bak),probs=seq(0,1,by=0.001)),seq(0,1,by=0.001),cex=0.5)
  # points(qgamma(seq(0,1,by=0.001),rate=0.1,shape=0.1),seq(0,1,by=0.001),cex=0.5,col="red")
  # points(qgamma(seq(0,1,by=0.001),rate=gammapar_obs[2],shape=gammapar_obs[1]),seq(0,1,by=0.001),cex=0.5,col="gold")
  if ( !argv$ensip.no_data_transf ) {
    if ( !is.na( argv$ensip.shape) & !is.na( argv$ensip.rate)) {
      shape <- argv$ensip.shape
      rate  <- argv$ensip.rate
    } else {
      gammapar_back <- array( data=NA, dim=c(2,nens))
      for (e in 1:nens) {
        res <- gamma_get_shape_rate_from_dataset_constrOptim( Xb[,e], small_const=small_const)
        gammapar_back[1,e] <- res$shape; gammapar_back[2,e] <- res$rate
      }
      ix <- which( gammapar_back[1,]>0 & gammapar_back[2,]>0)
      shape <- mean( gammapar_back[1,ix])
      rate  <- mean( gammapar_back[2,ix])
      rm( res, ix, gammapar_back)
    }
    yo <- gamma_anamorphosis( yo, shape=shape, rate=rate, small_const=small_const)
    for (e in 1:nens) {
      Xb[,e] <- gamma_anamorphosis( Xb[,e], shape=shape, rate=rate, small_const=small_const)
      Yb[,e] <- gamma_anamorphosis( Yb[,e], shape=shape, rate=rate, small_const=small_const)
    }
    xb <- apply( Xb, MAR=1, FUN=mean)
    yb <- apply( Yb, MAR=1, FUN=mean)
  }
  # -~- best estimate of the truth from the background -~-
  backg_best <- suppressWarnings( as.integer( argv$ensip.backg_best))
  if ( is.na( backg_best)) {
    if ( argv$ensip.backg_best == "mean" ) {
      yb_best <- yb
      xb_best <- xb
      cat("background is the ensemble mean\n")
    } else if ( argv$ensip.backg_best == "best" ) {
      cost <- vector( mode="numeric", length=nens)
      for (e in 1:nens) cost[e] <- sum( (Yb[,e]-yo)**2)
      e_best  <- which.min(cost)
      yb_best <- Yb[,e_best]
      xb_best <- Xb[,e_best]
      cat(paste("background is the ensemble member",argv$iff_fg.e[e_best],"\n"))
    }
  } else {
    e_best <- which( argv$iff_fg.e == backg_best )
    yb_best <- Yb[,e_best]
    xb_best <- Xb[,e_best]
    cat(paste("background is the ensemble member",backg_best,"\n"))
  }
  # -~- spatial analysis, Hybrid Ensemble Optimal Interpolation -~-
  # set the henoi parameters
  henoi_Dh     <- vector(mode="numeric",length=ngrid)
  henoi_Dh_loc <- vector(mode="numeric",length=ngrid)
  henoi_eps2   <- vector(mode="numeric",length=ngrid)
  henoi_alpha  <- vector(mode="numeric",length=ngrid)
  var_o_coeff  <- vector(mode="numeric",length=nobs)
  henoi_Dh_loc[] <- argv$ensip.henoi_Dh_loc
  henoi_Dh[]     <- argv$ensip.henoi_Dh
  henoi_eps2[]   <- argv$ensip.henoi_eps2
  henoi_alpha[]  <- argv$ensip.henoi_alpha
  var_o_coeff[]  <- argv$ensip.var_o_coeff
  if (!argv$ensip.henoi_par_notadaptive) {
    rmaster_agg <- aggregate( rmaster, fact=argv$ensip.henoi_reflen_aggfact)
    xy_agg <- xyFromCell( rmaster_agg, 1:ncell(rmaster_agg))
    grid_x_bak <- grid_x; grid_y_bak <- grid_y; ngrid_bak <- ngrid
    grid_x <- xy_agg[,1]; grid_y <- xy_agg[,2]; ngrid <- length(grid_x)
    if (argv$verbose) {t0<-Sys.time(); cat("closest obs ... ")}
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( distance_closest_k_obs,
                          1:ngrid,
                          SIMPLIFY = T,
                          mc.cores = argv$cores,
                          k = argv$ensip.henoi_reflen_k))
    } else {
      res <- t( mapply( distance_closest_k_obs,
                        1:ngrid,
                        SIMPLIFY = T,
                        k = argv$ensip.henoi_reflen_k))
    }
    if (argv$verbose) {
      t1<-Sys.time()
      cat(paste("time",round(t1-t0,1),attr(t1-t0,"unit"),"\n"))
    }
    grid_x<-grid_x_bak; grid_y<-grid_y_bak; ngrid<-ngrid_bak
    rm(grid_x_bak,grid_y_bak)
    rmaster_agg[] <- res[1,]
    rmaster_agg[] <- pmax( pmin( getValues( rmaster_agg), argv$ensip.henoi_reflen_max), 
                      argv$ensip.henoi_reflen_min)
    r <- resample( focal( rmaster_agg, 
           w=matrix( 1, argv$ensip.henoi_reflen_aggfact, argv$ensip.henoi_reflen_aggfact),
           fun=mean, na.rm=T, pad=T),
           rmaster, method="bilinear")
    henoi_Dh[] <- pmax( pmin( getValues(r)[aix], argv$ensip.henoi_reflen_max), 
                                                 argv$ensip.henoi_reflen_min)
    rm( r, rmaster_agg, res, xy_agg)
  }
  # henoi on the grid
  if (!is.na(argv$off_x)) {
    Af     <- Xb - xb
    HAf    <- Yb - yb
    Pfdiag <- apply( Af, MAR=1, FUN=var)
    if (argv$verbose) { t0<-Sys.time(); cat("x henoi ... ") }
    if (!is.na(argv$cores)) {
      res<-t( mcmapply( henoi,
                        1:ngrid,
                        SIMPLIFY = T,
                        mc.cores = argv$cores,
                        pmax     = argv$ensip.pmax,
                        rloc_min = argv$ensip.rloc_min,
                        statcov_backg   = argv$ensip.henoi_statcov_backg,
                        nens     = nens))
    } else {
      res<-t( mapply( henoi,
                      1:ngrid,
                      SIMPLIFY = T,
                      pmax     = argv$ensip.pmax,
                      rloc_min = argv$ensip.rloc_min,
                      statcov_backg   = argv$ensip.henoi_statcov_backg,
                      nens     = nens))
    }
    if (argv$verbose) {
      t1<-Sys.time()
      cat(paste("time",round(t1-t0,1),attr(t1-t0,"unit"),"\n"))
    }
    xa_henoi_mean <- res[,1]
    xa_henoi_var  <- res[,2]
    xidi_res      <- res[,3]
    xb_henoi_varu <- res[,4]
    rm( res, Pfdiag, Af, HAf)
    # data back-transformation (default = 0 everywhere)
    xa_xpv <- vector( mode="numeric", length=ngrid); xa_xpv[] <- 0
    xa_var <- vector( mode="numeric", length=ngrid); xa_var[] <- NA
    # pdf parameters. gamma, 1-shape, 2-rate / gauss, 1-mean, 2-var
    xa_pdf_par <- array( data=0, dim=c(ngrid,2)) 
    # Gaussian anamorphosis
    if ( !argv$ensip.no_data_transf ) {
      if (argv$verbose) { t0<-Sys.time(); cat("inv gamma ... ") }
      if ( length(shape) != ngrid) { shape<-rep(shape[1],ngrid); rate<-rep(rate[1],ngrid) }
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( inv_gamma_anamorphosis_constrOptim, 
                            1:ngrid,
                            SIMPLIFY = T,
                            mc.cores = argv$cores,
                            rr_inf = argv$rrinf))
      } else {
        res <- t( mapply( inv_gamma_anamorphosis_constrOptim,
                          1:ngrid,
                          SIMPLIFY = T,
                          rr_inf = argv$rrinf))
      }
      xa_pdf_par[,1] <- res[,1]
      xa_pdf_par[,2] <- res[,2]
      xa_xpv <- res[,3] 
      xa_var <- res[,4]
      rm( res)
      if (argv$verbose) {
        t1<-Sys.time()
        cat(paste("time",round(t1-t0,1),attr(t1-t0,"unit"),"\n"))
      }
    # no back-transformation
    } else {
      xa_henoi_mean[xa_henoi_mean<0] <- 0
      xa_henoi_var[xa_henoi_var<0]   <- 0
      # the next steps are consistent with "inv_gamma_anamorphosis_constrOptim"
      # default: set everything to henoi results
      xa_xpv <- xa_henoi_mean
      xa_var <- xa_henoi_var
      xa_pdf_par[,1] <- xa_henoi_mean 
      xa_pdf_par[,2] <- xa_henoi_var
      # exception 1: set to no-prec where expected val and uncertainty are small
      if ( length( ix <- which( xa_henoi_mean < argv$rrinf & sqrt(xa_henoi_var) < (argv$rrinf/2) ))>0) {
        xa_xpv[ix] <- 0
        xa_var[ix] <- NA
        xa_pdf_par[ix,1] <- NA
        xa_pdf_par[ix,2] <- NA
      }
      # exception 2: do not trust variance where is most likely raining but uncertainty is too small
      if ( length( ix <- which( xa_henoi_mean >=argv$rrinf & sqrt(xa_henoi_var) < (argv$rrinf/2) ))>0) {
        xa_xpv[ix] <- xa_henoi_mean
        xa_var[ix] <- NA
        xa_pdf_par[ix,1] <- NA
        xa_pdf_par[ix,2] <- NA
      }
    }
    save.image("tmp1.rdata")
    rm( xa_henoi_mean, xa_henoi_var)
  } # end henoi on the grid
  # -~- CV statistics -~-
  if (cv_mode) {
    # init
    grid_x_avbak <- grid_x
    grid_y_avbak <- grid_y
    ngrid_avbak  <- ngrid
    # set background at CV locations (Xb is transformed)
    res <- obsop_precip( r=rmaster, xcoord=VecX_cv, ycoord=VecY_cv, inf=-1000)
    if ( ( ngrid <- length( res$ix)) == 0 ) { boom( "found no valid Xb at cv points!") }
    ix_cv <- res$ix
    Xb <- res$yb
    xb <- rowMeans( Xb)
    grid_x <- VecX_cv[ix_cv]
    grid_y <- VecY_cv[ix_cv]
    rm(res)
    # set henoi parameters in a consistent way as for the grid
    henoi_Dh_avbak     <- henoi_Dh 
    henoi_Dh_loc_avbak <- henoi_Dh_loc
    henoi_eps2_avbak   <- henoi_eps2
    henoi_alpha_avbak   <- henoi_alpha
    r<-rmaster; r[]<-NA
    r[aix]       <- henoi_Dh
    henoi_Dh     <- extract( r, cbind( grid_x, grid_y))
    r[aix]       <- henoi_Dh_loc
    henoi_Dh_loc <- extract( r, cbind( grid_x, grid_y))
    r[aix]       <- henoi_eps2
    henoi_eps2   <- extract( r, cbind( grid_x, grid_y))
    r[aix]       <- henoi_alpha
    henoi_alpha   <- extract( r, cbind( grid_x, grid_y))
    rm(r)
    # spatial analysis
    HAf <- Yb - yb
    Af <- HAf
    Pfdiag <- apply(  Af, MAR=1, FUN=var)
    if (argv$verbose) { t0<-Sys.time(); cat("cv henoi ... ") }
    if (!is.na(argv$cores)) {
      res <- t( mcmapply( henoi, 
                          1:ngrid,
                          SIMPLIFY = T,
                          mc.cores = argv$cores,
                          pmax     = argv$ensip.pmax,
                          rloc_min = argv$ensip.rloc_min,
                          statcov_backg   = argv$ensip.henoi_statcov_backg,
                          nens     = nens))
    } else {
      res <- t( mapply( henoi,
                        1:ngrid,
                        SIMPLIFY = T,
                        pmax     = argv$ensip.pmax,
                        rloc_min = argv$ensip.rloc_min,
                        statcov_backg   = argv$ensip.henoi_statcov_backg,
                        nens     = nens))
    }
    if (argv$verbose) {
      t1<-Sys.time()
      cat(paste("time",round(t1-t0,1),attr(t1-t0,"unit"),"\n"))
    }
    yav_henoi_mean <- res[,1]
    yav_henoi_var  <- res[,2]
    yidiv          <- res[,3]
    yav_henoi_varu <- res[,4]
    rm( res, Pfdiag, Af, HAf)
    # data back-transformation (default = 0 everywhere)
    yav_xpv <- vector( mode="numeric", length=ngrid); yav_xpv[] <- 0
    yav_var <- vector( mode="numeric", length=ngrid); yav_var[] <- NA
    # pdf parameters. gamma, 1-shape, 2-rate / gauss, 1-mean, 2-var
    yav_pdf_par <- array( data=0, dim=c(ngrid,2)) 
    # Gaussian anamorphosis
    if ( !argv$ensip.no_data_transf ) {
      if (argv$verbose) { t0<-Sys.time(); cat("inv gamma ... ") }
      if ( length(shape) != ngrid) { shape<-rep(shape[1],ngrid); rate<-rep(rate[1],ngrid) }
      if (!is.na(argv$cores)) {
        res <- t( mcmapply( inv_gamma_anamorphosis_constrOptim, 
                            1:ngrid,
                            SIMPLIFY = T,
                            mc.cores = argv$cores,
                            rr_inf = argv$rrinf,
                            cv_mode  = T))
      } else {
        res <- t( mapply( inv_gamma_anamorphosis_constrOptim,
                          1:ngrid,
                          SIMPLIFY = T,
                          rr_inf = argv$rrinf,
                          cv_mode  = T))
      }
      yav_pdf_par[,1] <- res[,1]
      yav_pdf_par[,2] <- res[,2]
      yav_xpv <- res[,3]
      yav_var <- res[,4]
      rm( res)
      if (argv$verbose) {
        t1<-Sys.time()
        cat(paste("time",round(t1-t0,1),attr(t1-t0,"unit"),"\n"))
      }
    # no back-transformation
    } else {
      yav_henoi_mean[yav_henoi_mean<0] <- 0
      yav_henoi_var[yav_henoi_var<0]   <- 0
      # the next steps are consistent with "inv_gamma_anamorphosis_constrOptim"
      # default: set everything to henoi results
      yav_xpv <- yav_henoi_mean
      yav_var <- yav_henoi_var
      yav_pdf_par[,1] <- yav_henoi_mean 
      yav_pdf_par[,2] <- yav_henoi_var
      # exception 1: set to no-prec where expected val and uncertainty are small
      if ( length( ix <- which( yav_henoi_mean < argv$rrinf & sqrt(yav_henoi_var) < argv$rrinf ))>0) {
        yav_xpv[ix] <- 0
        yav_var[ix] <- NA
        yav_pdf_par[ix,1] <- NA
        yav_pdf_par[ix,2] <- NA
      }
      # exception 2: do not trust variance where is most likely raining but uncertainty is too small
      if ( length( ix <- which( yav_henoi_mean >=argv$rrinf & sqrt(yav_henoi_var) < argv$rrinf ))>0) {
        yav_xpv[ix] <- yav_henoi_mean
        yav_var[ix] <- NA
        yav_pdf_par[ix,1] <- NA
        yav_pdf_par[ix,2] <- NA
      }
    }
    rm( yav_henoi_mean, yav_henoi_var)
    # prepare for output
    ya_cv <- vector( mode="numeric", length=ncv); ya_cv[] <- NA
    yb_cv       <- ya_cv # NAs
    ya_cv_var   <- ya_cv # NAs
    yb_cv       <- ya_cv # NAs
    yidi_cv     <- ya_cv # NAs
    ya_cv_henoi_varu <- ya_cv # NAs
    ya_cv_pdf_par <- array( data=NA, dim=c(ncv,2))
    Yb_cv         <- array( data=NA, dim=c(ncv,nens))
    ya_cv[ix_cv]     <- yav_xpv
    ya_cv_var[ix_cv] <- yav_var
    yidi_cv[ix_cv]   <- yidiv
    ya_cv_henoi_varu[ix_cv]    <- yav_henoi_varu
    ya_cv_pdf_par[ix_cv,] <- yav_pdf_par
    Xb <- Xb_bak 
    Yb_cv <- obsop_precip( r=rmaster, xcoord=VecX_cv, ycoord=VecY_cv)$yb
    yb_cv <- rowMeans( Yb_cv)
    save.image("tmp2.rdata")
    rm(ix_cv, yav_xpv, yav_var, yidiv, yav_henoi_varu, yav_pdf_par)
    # restore old grid variables 
    grid_x <- grid_x_avbak
    grid_y <- grid_y_avbak
    ngrid  <- ngrid_avbak
    yav_henoi_Dh     <- henoi_Dh
    yav_henoi_Dh_loc <- henoi_Dh_loc
    yav_henoi_eps2   <- henoi_eps2
    yav_henoi_alpha  <- henoi_alpha
    henoi_Dh     <- henoi_Dh_avbak
    henoi_Dh_loc <- henoi_Dh_loc_avbak
    henoi_eps2   <- henoi_eps2_avbak
    henoi_alpha   <- henoi_alpha_avbak
    rm(list = ls()[grep("_avbak$", ls())])
  } # end of cv_mode
  # -~- prepare for output at observation point -~-
  # init
  ya_pdf_par <- array( data=NA, dim=c(n0,2))
  yb <- vector( mode="numeric", length=n0); yb[]<-NA
  ya       <- yb #NAs
  ya_var   <- yb #NAs
  yidi     <- yb #NAs
  ya_henoi_varu <- yb #NAs
  # observations
  yo <- yo_bak
  # background ensemble mean at observation point
  yb[ix_obs] <- rowMeans( Yb_bak)
  rm(ix_obs)
  # fill in only if the fields have been created
  if ( !is.na( argv$off_x)) {
    # analysis at observation points
    r <- rmaster
    r[aix] <- xa_xpv
    ya     <- extract( r, cbind( VecX, VecY))
    r[aix] <- xa_var
    ya_var <- extract( r, cbind( VecX, VecY))
    # idi at observation points
    r[aix] <- xidi_res
    yidi   <- extract( r, cbind( VecX, VecY))
    # unexplained variance at observation points
    r[aix]   <- xb_henoi_varu
    ya_henoi_varu <- extract( r, cbind( VecX, VecY))
    #
    r[aix]         <- xa_pdf_par[,1]
    ya_pdf_par[,1] <- extract(r,cbind(VecX,VecY))
    r[aix]         <- xa_pdf_par[,2]
    ya_pdf_par[,2] <- extract(r,cbind(VecX,VecY))
    rm(r)
  }
  # -~- prepare for gridded output -~-
  if ( !is.na( argv$off_x) ) {
    xa <- xa_xpv
    xa_errsd <- xa_var; xa_errsd[] <- NA
    if ( length( ix <- which( xa_var < 0 ))) {
      xa_errsd[ix] <- xa_var[ix]
      xa_errsd[-ix] <- sqrt( xa_var[-ix])
      cat(paste("@@ Warning: negative variance found for ",length(ix),"points.\n"))
      cat("For those points, the negative variance is reported as output.\n")
    } else {
      xa_errsd <- sqrt( xa_var)
    }
    xb   <- Xb_bak
    xidi <- xidi_res
    if ( length( ix <- which( xidi < 0))) 
      cat(paste("@@ Warning: negative IDI found for ",length(ix),"points.\n"))
    xdh  <- henoi_Dh_loc
    rm(list = ls()[grep("_bak$", ls())])
    rm( xa_xpv, xa_var, henoi_Dh_loc, xidi_res)
  }
