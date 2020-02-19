  if (argv$verbose) t00<-Sys.time()
  # save original vectors for future use
  VecX_orig<-VecX
  VecY_orig<-VecY
  VecZ_orig<-VecZ
  prId_orig<-prId
  yo_orig<-yo
  xb_orig<-xb
  aix<-which(!is.na(xgrid) & !is.na(ygrid) & !is.na(dem))
  if ((argv$cv_mode|argv$cv_mode_random)) {
    xgrid<-VecX_cv
    ygrid<-VecY_cv
  } else {
    xgrid_orig<-xgrid
    ygrid_orig<-ygrid
    dem_orig<-dem
    # define grid-space
    # note: xb is defined only over "aix" points
    xgrid<-xgrid[aix]
    ygrid<-ygrid[aix]
    dem<-dem[aix]
  }
  # define obs-space, background
  r<-rmaster
  yb0<-array(data=NA,dim=c(n0,nens))
  argv$letkf.obsop<-"obsop_LapseRateConst"
  argv$letkf.obsop.mMinElevDiff<-30
  argv$letkf.obsop.mMinGradient<-(-0.01)
  argv$letkf.obsop.mMaxGradient<-0.01
  argv$letkf.obsop.mSearchRadius<-10
  for (e in 1:nens) {
    if (argv$letkf.obsop=="obsop_LapseRateConst") {
      yb0[,e]<-obsop_LapseRateConst(s=n0,
                                    m=length(aix),
                                    xm=xgrid_orig[aix],
                                    ym=ygrid_orig[aix],
                                    zm=dem_orig[aix],
                                    xs=VecX,
                                    ys=VecY,
                                    zs=VecZ,
                                    ifield=xb[,e],
                                    oval=rep(0,n0),
                                    mMinElevDiff=argv$letkf.obsop.mMinElevDiff,
                                    mMinGradient=argv$letkf.obsop.mMinGradient,
                                    mMaxGradient=argv$letkf.obsop.mMaxGradient,
                                    mSearchRadius=argv$letkf.obsop.mSearchRadius)$oval
    }
  }
  ix<-which(apply(yb0,MAR=1,
                  FUN=function(x){length(which(!is.na(x) & !is.nan(x)))})
            ==nens)
  n1<-length(ix)
  yb<-array(data=NA,dim=c(n1,nens))
  for (e in 1:nens) yb[,e]<-yb0[ix,e]
  rm(yb0)
  # define obs-space, observations
  ix_orig<-ix
  yo<-yo[ix]
  VecX<-VecX[ix]
  VecY<-VecY[ix]
  VecZ<-VecZ[ix]
  prId<-prId[ix] 
  rm(ix)
  # define the inverse of the diagonal observation error matrix
  argv$letkf.eps2_default<-0.5
  argv$letkf.eps2_prId<-NA
  diagRinv<-rep(argv$letkf.eps2_default,n1)
  for (i in 1:length(argv$letkf.eps2_prId)) {
    ix<-which(prId==argv$letkf.eps2_prId[i])
    if (length(ix)>0) diagRinv[ix]<-argv$letkf.eps2[i]
    rm(ix)
  }
  if (any(diagRinv==0 | is.na(diagRinv) | is.nan(diagRinv))) 
    boom("ERROR in the definition of Rinv")
  diagRinv<-1./diagRinv
  # if CVmode, then adjust the background
  if ((argv$cv_mode|argv$cv_mode_random)) {
    xb<-array(data=NA,dim=c(ncv,nens))
    for (e in 1:nens) {
      if (argv$letkf.obsop=="obsop_LapseRateConst") {
        xb[,e]<-obsop_LapseRateConst(s=n0,
                                      m=length(aix),
                                      xm=xgrid_orig[aix],
                                      ym=ygrid_orig[aix],
                                      zm=dem_orig[aix],
                                      xs=VecX_cv,
                                      ys=VecY_cv,
                                      zs=VecZ_cv,
                                      ifield=xb_orig[,e],
                                      oval=rep(0,n0),
                                      mMinElevDiff=argv$letkf.obsop.mMinElevDiff,
                                      mMinGradient=argv$letkf.obsop.mMinGradient,
                                      mMaxGradient=argv$letkf.obsop.mMaxGradient,
                                      mSearchRadius=argv$letkf.obsop.mSearchRadius)$oval
      }
    }
  }
  # ensemble mean
  xbm<-apply(xb,MAR=1,FUN=mean)
  ybm<-apply(yb,MAR=1,FUN=mean)
  # ensemble perturbations
  Xb<-xb-xbm
  Yb<-yb-ybm
  # ensemble variances
  ybvar<-apply(Yb,MAR=1,FUN=sd)**2
  # intialize analysis vector
  xa<-xb
  xa[]<-NA
  # letkf: loop over gridpoints
  argv$letkf.Dh<-10
  argv$letkf.rloc_min<-0.0013
  argv$letkf.pmax<-200
  argv$letkf.sigma2_min<-0.1
  if (argv$verbose) t00<-Sys.time()
  for (i in 1:length(xgrid)) {
    dist<-(((xgrid[i]-VecX)**2+(ygrid[i]-VecY)**2 )**0.5)/1000.
    rloc<-exp(-0.5*(dist**2./argv$letkf.Dh**2))
    sel<-which(rloc>argv$letkf.rloc_min)
    # i-th gridpoint analysis 
    if (length(sel)>0) {
#      if (length(sel)>argv$letkf.pmax) sel<-order(dist,decreasing=T)[1:argv$letkf.pmax]
      if (length(sel)>argv$letkf.pmax) sel<-order(dist,decreasing=F)[1:argv$letkf.pmax]
      sigma2<-max(c(mean(ybvar[sel]),argv$letkf.sigma2_min))
      Yb.i<-Yb[sel,,drop=F]
      d.i<-yo[sel]-ybm[sel]
      C.i<-t(1./sigma2*diagRinv[sel]*rloc[sel]*Yb.i)
      C1.i<-crossprod(t(C.i),Yb.i)
      Cd.i<-crossprod(t(C.i),d.i)
      rm(C.i)
      diag(C1.i)<-diag(C1.i) + (nens-1)
      Pa.i<-chol2inv(chol(C1.i))
      rm(C1.i)
      # Pa*(k-1) is the matrix for which we want to square root:
      a.eig <- eigen(Pa.i*(nens-1),symmetric=T)
      Wa <- tcrossprod(
             tcrossprod( a.eig$vectors, diag(sqrt(a.eig$values))),
             t(solve(a.eig$vectors)) )
      waa<-crossprod(Pa.i, Cd.i )
      rm(Cd.i)
      W<-Wa+as.vector(waa)
      xa[i,]<-xbm[i]+Xb[i,] %*% W
      rm(W,Wa,waa,a.eig)
    # i-th gridpoint is isolated
    } else {
      xa[i,]<-xb[i,]
    }
  } # end: letkf loop over gridpoints
  if (argv$verbose) {
    t11<-Sys.time()
    print(paste("letkf time=",round(t11-t00,1),attr(t11-t00,"unit")))
  }
  if (!(argv$cv_mode|argv$cv_mode_random) & 
      ("idi" %in% argv$off_x.variables | 
       argv$idiv_instead_of_elev)) {
    if (argv$verbose) t00<-Sys.time()
    D<-exp(-0.5*((outer(VecY,VecY,FUN="-")**2.+
                  outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
                 /argv$hyletkf.Dh_oi)**2.)
    diag(D)<-diag(D)+argv$hyletkf.eps2_oi
    InvD<-chol2inv(chol(D))
    rm(D)
    xidi<-OI_RR_fast(yo=rep(1,length(VecX)),
                     yb=rep(0,length(VecX)),
                     xb=rep(0,length(xgrid)),
                     xgrid=xgrid,
                     ygrid=ygrid,
                     VecX=VecX,
                     VecY=VecY,
                     Dh=argv$hyletkf.Dh_oi)
    if (argv$idiv_instead_of_elev) {
      G<-exp(-0.5*((outer(VecY_orig,VecY,FUN="-")**2.+
                    outer(VecX_orig,VecX,FUN="-")**2.)**0.5/1000.
                    /argv$hyletkf.Dh_oi)**2.)
      W<-tcrossprod(G,InvD)
      rm(G,InvD)
      # this is the cross-validation integral data influence ("yidiv")
      elev_for_verif<-rep(1,n0) + 1./(1.-diag(W)) * (rowSums(W)-rep(1,n0))
      rm(W)
    }
    # InvD could be used in the following
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("idi time=",round(t11-t00,1),
                              attr(t11-t00,"unit")))
    }
  }
  if (!argv$idiv_instead_of_elev) elev_for_verif<-VecZ
  # back to the original vectors
  yo<-yo_orig
  prId<-prId_orig
  VecX<-VecX_orig
  VecY<-VecY_orig
  rm(yo_orig,prId_orig,VecX_orig,VecY_orig)
  # prepare point / gridded output
  if (!(argv$cv_mode|argv$cv_mode_random)) {
    ya<-array(data=NA,dim=c(n0,length(argv$iff_fg.e)))
    yb<-array(data=NA,dim=c(n0,length(argv$iff_fg.e)))
    for (j in 1:length(argv$off_x.variables)) {
      if (!exists("r.list")) r.list<-list()
      if (argv$off_x.variables[j]=="analysis") {
        r<-rmaster
        r[]<-NA
        grid<-array(data=NA,dim=c(length(x),length(y),length(argv$iff_fg.e),1))
        for (e in 1:length(argv$iff_fg.e)) {
          if (argv$iff_fg.e[e] %in% valens) {
            i<-which(valens==argv$iff_fg.e[e])
            r[aix]<-xa[,i]
            grid[,,e,1]<-matrix(data=getValues(r),
                                ncol=length(y),
                                nrow=length(x))
            ya[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
            auxx<-which(ya[,e]<0 | is.na(ya[,e]) | is.nan(ya[,e]))
            if (length(auxx)>0) 
              ya[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
            rm(auxx)
          }
        }
        rm(r)
        r.list[[j]]<-grid
        rm(grid)
      } else if (argv$off_x.variables[j]=="background") {
        r<-rmaster
        r[]<-NA
        grid<-array(data=NA,dim=c(length(x),length(y),length(argv$iff_fg.e),1))
        for (e in 1:length(argv$iff_fg.e)) {
          if (argv$iff_fg.e[e] %in% valens) {
            i<-which(valens==argv$iff_fg.e[e])
            r[aix]<-xb[,i]
            grid[,,e,1]<-matrix(data=getValues(r),
                                ncol=length(y),
                                nrow=length(x))
            yb[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
            auxx<-which(yb[,e]<0 | is.na(yb[,e]) | is.nan(yb[,e]))
            if (length(auxx)>0) 
              yb[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
            rm(auxx)
          }
        }
        rm(r)
        r.list[[j]]<-grid
        rm(grid)
      } else if (argv$off_x.variables[j]=="idi") {
        r<-rmaster
        r[]<-NA
        r[aix]<-xidi
        r.list[[j]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      }
    }
  } # end prepare for gridded output
# END of Local Ensemble Transform Kalman Filter
