  if (argv$verbose) {
    print("OI_firstguess")
  }
  # set eps2 (station-dependent and value-dependent)
  eps2<-vector(mode="numeric",length=n0)
  ixdef<-which(is.na(argv$oifg.eps2_prId))
  for (r in 1:(length(argv$oifg.eps2_r)-1)) {
    ixr<-which(yo>=argv$oifg.eps2_r[r] & 
               yo<argv$oifg.eps2_r[(r+1)])
    eps2[ixr]<-argv$oifg.eps2[r,ixdef]
    for (i in 1:length(argv$oifg.eps2_prId)) {
      if (is.na(argv$oifg.eps2_prId[i])) next
      ixr<-which(yo>=argv$oifg.eps2_r[r] & 
                 yo<argv$oifg.eps2_r[(r+1)] & 
                 prId==argv$oifg.eps2_prId[i]) 
      eps2[ixr]<-argv$oifg.eps2[r,i]
    }
  }
# something like this will create problems in the S matrix inversion
#  dh<-vector(mode="numeric",length=n0)
#  ixdef<-which(is.na(argv$oifg.Dh_prId))
#  for (r in 1:(length(argv$oifg.Dh_r)-1)) {
#    ixr<-which(yo>=argv$oifg.Dh_r[r] & 
#               yo<argv$oifg.Dh_r[(r+1)])
#    dh[ixr]<-argv$oifg.Dh[r,ixdef]
#    for (i in 1:length(argv$oifg.Dh_prId)) {
#      if (is.na(argv$oifg.Dh_prId[i])) next
#      ixr<-which(yo>=argv$oifg.Dh_r[r] & 
#                 yo<argv$oifg.Dh_r[(r+1)] & 
#                 prId==argv$oifg.Dh_prId[i]) 
#      dh[ixr]<-argv$oifg.Dh[r,i]
#    }
#  }
  # dh (m), one for all gridpoints
  dh<-argv$oifg.Dh*1000
  dh2<-dh*dh
  # background at station locations
  yb<-extract(rfg, cbind(VecX,VecY), method="bilinear")
  # if CVmode, then fix the background
  if (cv_elab) 
    yb_cv<-extract(rfg, cbind(VecX_cv,VecY_cv), method="bilinear")
  rm(rfg)
  # -- elaboration on the grid, x_elab --
  if (x_elab) {
    t00<-Sys.time()
    xgrid_spint<-xgrid[aix]
    ygrid_spint<-ygrid[aix]
    if (argv$transf=="Box-Cox") {
      xb_spint<-xb
      yo_spint<-yo
      yb_spint<-yb
      xb_spint[which(xb_spint<0)]<-0
      yb_spint[which(yb_spint<0)]<-0
      yo_spint[which(yo_spint<0)]<-0
      yo_spint<-boxcox(yo_spint,argv$transf.boxcox_lambda)
      yb_spint<-boxcox(yb_spint,argv$transf.boxcox_lambda)
      xb_spint<-boxcox(xb_spint,argv$transf.boxcox_lambda)
    } else {
      xb_spint<-xb
      yo_spint<-yo
      yb_spint<-yb
    }
    # multicores
    # xa=arr[,1], xa_errvar=arr[,2], o_errvar=arr[,3], xidi=arr[,4], idiv=arr[,5]
    if (!is.na(argv$cores)) {
      arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                      1:ngrid,
                      mc.cores=argv$cores,
                      SIMPLIFY=T,
                      pmax=argv$pmax,
                      corr=argv$corrfun))
    # no-multicores
    } else {
      arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                    1:ngrid,
                    SIMPLIFY=T,
                    pmax=argv$pmax,
                    corr=argv$corrfun))
    }
    xidi<-arr[,4]
    argv$oifg.xta_errvar_smooth<-50000
    if (argv$transf=="Box-Cox") {
      xta<-arr[,1]
      # smooth out xta_errvar 
      # first, fill NAs on
      r<-rmaster;r[]<-NA
      r[aix]<-arr[,2]
      r1<-r
      focalx<-ceiling(argv$oifg.xta_errvar_smooth/xres(rmaster))
      focaly<-ceiling(argv$oifg.xta_errvar_smooth/yres(rmaster))
      if ((focalx%%2)==0) focalx<-focalx+1
      if ((focaly%%2)==0) focaly<-focaly+1
      for (ii in 1:10) {
        r1<-focal(r1,w=matrix(1,focalx,focaly),fun=mean,na.rm=T,NAonly=T)
        if (!any(is.na(getValues(r1)[aix]))) break
      }
      # second, smooth and make a continuous transition with non-prec regions
      for (ii in 1:3) {
        if (ii==2) {focalx<-5;focaly<-5}
        if (ii==3) {focalx<-3;focaly<-3}
        aux<-getValues(r1)[aix]
        aux[is.na(aux)|xta<boxcox(argv$rrinf,argv$transf.boxcox_lambda)]<-0
        r1[aix]<-aux
        r1[is.na(r1)]<-0
        r1<-focal(r1,w=matrix(1,focalx,focaly),fun=mean,na.rm=T)
      }
      xta_errvar<-getValues(r1)[aix]; rm(r,r1,focalx,focaly,aux)
      #
      xtb<-xb_spint
      xbctmp<-cbind(xta,sqrt(abs(xta_errvar)))
      if (!is.na(argv$cores)) {
        res<-t(mcmapply(tboxcox,
                        1:ngrid,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        x_mean_sd_in_xbctmp=T,
                        lambda=argv$transf.boxcox_lambda,
                        threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                        distribution=T,
                        statistics="mean_sd",
                        method="Taylor"))
      } else {
        res<-t(mapply(tboxcox,
                      1:ngrid,
                      SIMPLIFY=T,
                      x_mean_sd_in_xbctmp=T,
                      lambda=argv$transf.boxcox_lambda,
                      threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                      distribution=T,
                      statistics="mean_sd",
                      method="Taylor"))
      }
      rm(xbctmp)
      xa<-res[,1]
      xa_errsd<-res[,2]
      rm(res)
    } else {
      xa<-arr[,1]
      xa_errsd<-sqrt(abs(arr[,2]))
    }
    rm(arr,xgrid_spint,ygrid_spint,xb_spint,yb_spint,yo_spint)
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("x_elab, time=",
                  round(t11-t00,1),attr(t11-t00,"unit")))
    }
  } # end x_elab
  #
  # -- elaboration on station points, y_elab --
  if (y_elab) {
    t00<-Sys.time()
    xgrid_spint<-VecX
    ygrid_spint<-VecY
    if (argv$transf=="Box-Cox") {
      xb_spint<-yb
      yo_spint<-yo
      yb_spint<-yb
      xb_spint[which(xb_spint<0)]<-0
      yb_spint[which(yb_spint<0)]<-0
      yo_spint[which(yo_spint<0)]<-0
      yo_spint<-boxcox(yo_spint,argv$transf.boxcox_lambda)
      yb_spint<-boxcox(yb_spint,argv$transf.boxcox_lambda)
      xb_spint<-boxcox(xb_spint,argv$transf.boxcox_lambda)
    } else {
      xb_spint<-yb
      yo_spint<-yo
      yb_spint<-yb
    }
    # multicores
    # xa=arr[,1], xa_errvar=arr[,2], o_errvar=arr[,3], xidi=arr[,4], idiv=arr[,5]
    if (!is.na(argv$cores)) {
      arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                      1:n0,
                      mc.cores=argv$cores,
                      SIMPLIFY=T,
                      y_elab=T,
                      pmax=argv$pmax,
                      corr=argv$corrfun))
    # no-multicores
    } else {
      arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                    1:n0,
                    SIMPLIFY=T,
                    y_elab=T,
                    pmax=argv$pmax,
                    corr=argv$corrfun))
    }
    yidi<-arr[,4]
    yidiv<-arr[,5]
    if (argv$idiv_instead_of_elev) {
      elev_for_verif_y<-yidiv
    } else {
      elev_for_verif_y<-VecZ
    }
    if (argv$transf=="Box-Cox") {
      yta<-arr[,1]
      yta_errvar<-arr[,2]
      ytb<-yb_spint
      yto<-yo_spint
      ytav<-arr[,6]
      xbctmp<-cbind(yta,sqrt(abs(yta_errvar)))
      if (!is.na(argv$cores)) {
        res<-t(mcmapply(tboxcox,
                        1:n0,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        x_mean_sd_in_xbctmp=T,
                        lambda=argv$transf.boxcox_lambda,
                        threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                        distribution=T,
                        statistics="mean_sd",
                        method="Taylor"))
      } else {
        res<-t(mapply(tboxcox,
                      1:n0,
                      SIMPLIFY=T,
                      x_mean_sd_in_xbctmp=T,
                      lambda=argv$transf.boxcox_lambda,
                      threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                      distribution=T,
                      statistics="mean_sd",
                      method="Taylor"))
      }
      rm(xbctmp)
      ya<-res[,1]
      ya_errsd<-res[,2]
      rm(res)
      # NOTE: if no transformation, yav is the leave-one-out cv.
      #       if transformation, yav is not the leave-one-out cv.
#      ya<-apply(cbind(yta,sqrt(abs(arr[,2]))),
#                      MARGIN=1,
#                      FUN=tboxcox4pdf_apply,
#                          lambda=argv$transf.boxcox_lambda,
#                          brrinf=brrinf)
      yav<-tboxcox(ytav,lambda=argv$transf.boxcox_lambda,distribution=F)
    } else {
      ya<-arr[,1]
      ya_errsd<-sqrt(abs(arr[,2]))
      yav<-arr[,6]
    }
    rm(arr,xgrid_spint,ygrid_spint,xb_spint,yb_spint,yo_spint)
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("cv_elab, time=",
                  round(t11-t00,1),attr(t11-t00,"unit")))
    }
  } # end y_elab
  #
  # -- elaboration on station points, cv_elab --
  if (cv_elab) {
    t00<-Sys.time()
    xgrid_spint<-VecX_cv
    ygrid_spint<-VecY_cv
    if (argv$transf=="Box-Cox") {
      xb_spint<-yb_cv
      yo_spint<-yo
      yb_spint<-yb
      xb_spint[which(xb_spint<0)]<-0
      yb_spint[which(yb_spint<0)]<-0
      yo_spint[which(yo_spint<0)]<-0
      yo_spint<-boxcox(yo_spint,argv$transf.boxcox_lambda)
      yb_spint<-boxcox(yb_spint,argv$transf.boxcox_lambda)
      xb_spint<-boxcox(xb_spint,argv$transf.boxcox_lambda)
    } else {
      xb_spint<-yb_cv
      yo_spint<-yo
      yb_spint<-yb
    }
    # multicores
    # xa=arr[,1], xa_errvar=arr[,2], o_errvar=arr[,3], xidi=arr[,4], idiv=arr[,5]
    if (!is.na(argv$cores)) {
      arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                      1:ncv,
                      mc.cores=argv$cores,
                      SIMPLIFY=T,
                      pmax=argv$pmax,
                      corr=argv$corrfun))
    # no-multicores
    } else {
      arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                    1:ncv,
                    SIMPLIFY=T,
                    pmax=argv$pmax,
                    corr=argv$corrfun))
    }
    yidi_cv<-arr[,4]
    if (argv$idiv_instead_of_elev) {
      elev_for_verif_cv<-yidi_cv
    } else {
      elev_for_verif_cv<-VecZ_cv
    }
    if (argv$transf=="Box-Cox") {
      yta_cv<-arr[,1]
      yta_errvar_cv<-arr[,2]
      ytb_cv<-xb_spint
      yto_cv<-boxcox(yo_cv,argv$transf.boxcox_lambda)
      xbctmp<-cbind(yta_cv,sqrt(abs(yta_errvar_cv)))
      if (!is.na(argv$cores)) {
        res<-t(mcmapply(tboxcox,
                        1:ncv,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        x_mean_sd_in_xbctmp=T,
                        lambda=argv$transf.boxcox_lambda,
                        threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                        distribution=T,
                        statistics="mean_sd",
                        method="Taylor"))
      } else {
        res<-t(mapply(tboxcox,
                      1:ncv,
                      SIMPLIFY=T,
                      x_mean_sd_in_xbctmp=T,
                      lambda=argv$transf.boxcox_lambda,
                      threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                      distribution=T,
                      statistics="mean_sd",
                      method="Taylor"))
      }
      rm(xbctmp)
      ya_cv<-res[,1]
      ya_errsd_cv<-res[,2]
      rm(res)
    } else {
      ya_cv<-arr[,1]
      ya_errsd_cv<-sqrt(abs(arr[,2]))
    }
    rm(arr,xgrid_spint,ygrid_spint,xb_spint,yb_spint,yo_spint)
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("cv_elab, time=",
                  round(t11-t00,1),attr(t11-t00,"unit")))
    }
  } # end cv_elab
  #
  # -- elaboration on station points, loocv_elab --
  if (loocv_elab) {
    t00<-Sys.time()
    xgrid_spint<-VecX
    ygrid_spint<-VecY
    if (argv$transf=="Box-Cox") {
      xb_spint<-yb
      yo_spint<-yo
      yb_spint<-yb
      xb_spint[which(xb_spint<0)]<-0
      yb_spint[which(yb_spint<0)]<-0
      yo_spint[which(yo_spint<0)]<-0
      yo_spint<-boxcox(yo_spint,argv$transf.boxcox_lambda)
      yb_spint<-boxcox(yb_spint,argv$transf.boxcox_lambda)
      xb_spint<-boxcox(xb_spint,argv$transf.boxcox_lambda)
    } else {
      xb_spint<-yb
      yo_spint<-yo
      yb_spint<-yb
    }
    # multicores
    # xa=arr[,1], xa_errvar=arr[,2], o_errvar=arr[,3], xidi=arr[,4], idiv=arr[,5]
    if (!is.na(argv$cores)) {
      arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                      1:n0,
                      mc.cores=argv$cores,
                      SIMPLIFY=T,
                      corr=argv$corrfun,
                      pmax=argv$pmax,
                      loocv=T,
                      y_elab=T))
    # no-multicores
    } else {
      arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                    1:n0,
                    SIMPLIFY=T,
                    corr=argv$corrfun,
                    pmax=argv$pmax,
                    loocv=T,
                    y_elab=T))
    }
    yidi_lcv<-arr[,4]
    yidiv_lcv<-arr[,5]
    if (argv$idiv_instead_of_elev) {
      elev_for_verif_loocv<-yidiv_lcv
    } else {
      elev_for_verif_loocv<-VecZ
    }
    if (argv$transf=="Box-Cox") {
      yta_lcv<-arr[,1]
      yta_errvar_lcv<-arr[,2]
      ytb_lcv<-yb_spint
      yto_lcv<-yo_spint
      xbctmp<-cbind(yta_lcv,sqrt(abs(yta_errvar_lcv)))
      if (!is.na(argv$cores)) {
        res<-t(mcmapply(tboxcox,
                        1:n0,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        x_mean_sd_in_xbctmp=T,
                        lambda=argv$transf.boxcox_lambda,
                        threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                        distribution=T,
                        statistics="mean_sd",
                        method="Taylor"))
      } else {
        res<-t(mapply(tboxcox,
                      1:n0,
                      SIMPLIFY=T,
                      x_mean_sd_in_xbctmp=T,
                      lambda=argv$transf.boxcox_lambda,
                      threshold=boxcox(argv$rrinf,argv$transf.boxcox_lambda),
                      distribution=T,
                      statistics="mean_sd",
                      method="Taylor"))
      }
      rm(xbctmp)
      ya_lcv<-res[,1]
      ya_errsd_lcv<-res[,2]
      rm(res)
    } else {
      ya_lcv<-arr[,1]
      ya_errsd_lcv<-sqrt(abs(arr[,2]))
    }
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("loocv_elab, time=",
                  round(t11-t00,1),attr(t11-t00,"unit")))
    }
    rm(arr,xgrid_spint,ygrid_spint,xb_spint,yb_spint,yo_spint)
  } # end loocv_elab
