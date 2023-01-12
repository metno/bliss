  if (argv$verbose) print("OI_Barnes")
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
    xobs_spint<-VecX
    yobs_spint<-VecY
    zobs_spint<-VecZ
    xb_spint<-xb
    yo_spint<-yo
    yb_spint<-yb
    nobs<-length(yo_spint)
    # multicores
    # xa=arr[,1], dh=arr[,2]
    if (!is.na(argv$cores)) {
      arr<-t(mcmapply(oi_bratseth_gridpoint_by_gridpoint,
                      1:ngrid,
                      mc.cores=argv$cores,
                      SIMPLIFY=T,
                      nSCloops=argv$oibr.nSCloops,
                      delta_hor=argv$oibr.delta_hor, #m
                      dh=argv$oibr.dh,
                      box_o_nearest_halfwidth=argv$oibr.box_o_nearest_halfwidth, #m
                      dz=argv$oibr.dz,
                      lafmin=argv$oibr.lafmin,
                      dh_adaptive=argv$oibr.dh_adaptive,
                      dh_adaptive_min=argv$oibr.dh_adaptive_min,
                      dh_adaptive_max=argv$oibr.dh_adaptive_max,
                      loocv=F,
                      pmax=argv$pmax,
                      corr=argv$corrfun))
    # no-multicores
    } else {
      arr<-t(mapply(oi_bratseth_gridpoint_by_gridpoint,
                    1:ngrid,
                    SIMPLIFY=T,
                    nSCloops=argv$oibr.nSCloops,
                    delta_hor=argv$oibr.delta_hor, #m
                    dh=argv$oibr.dh,
                    box_o_nearest_halfwidth=argv$oibr.box_o_nearest_halfwidth, #m
                    dz=argv$oibr.dz,
                    lafmin=argv$oibr.lafmin,
                    dh_adaptive=argv$oibr.dh_adaptive,
                    dh_adaptive_min=argv$oibr.dh_adaptive_min,
                    dh_adaptive_max=argv$oibr.dh_adaptive_max,
                    loocv=F,
                    pmax=argv$pmax,
                    corr=argv$corrfun))
    }
    xa<-arr[,1]
    xdh<-arr[,2]
    rm(arr,xgrid_spint,ygrid_spint,xb_spint,yb_spint,yo_spint)
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("x_elab, time=",
                  round(t11-t00,1),attr(t11-t00,"unit")))
    }
  } # end x_elab
  #
# END of SCB, Successive Correction Barnes' scheme (with background)  <===
