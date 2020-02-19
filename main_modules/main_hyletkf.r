  argv$hyletkf.Dh<-1000*argv$hyletkf.Dh # km to m
  argv$hyletkf.Dh_oi<-1000*argv$hyletkf.Dh_oi # km to m
  Dh2<-argv$hyletkf.Dh*argv$hyletkf.Dh
  Dh_oi2<-argv$hyletkf.Dh_oi*argv$hyletkf.Dh_oi
# xb, aix,nens
# yo, VecX,VecY, prId
  # save original vectors for future use
  VecX_orig<-VecX
  VecY_orig<-VecY
  prId_orig<-prId
  yo_orig<-yo
  xb_orig<-xb
  if ((argv$cv_mode|argv$cv_mode_random)) {
    xgrid<-VecX_cv
    ygrid<-VecY_cv
  } else {
    xgrid_orig<-xgrid
    ygrid_orig<-ygrid
    # define grid-space
    # note: xb is defined only over "aix" points
    xgrid<-xgrid[aix]
    ygrid<-ygrid[aix]
  }
  # background at observation locations
  r<-rmaster
  yb0<-array(data=NA,dim=c(n0,nens))
  for (e in 1:nens) {
    r[]<-NA
    r[aix]<-xb[,e]
    yb0[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
    auxx<-which(yb0[,e]<0 | is.na(yb0[,e]) | is.nan(yb0[,e]))
    if (length(auxx)>0) yb0[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
    rm(auxx)
  }
  rm(r)
  # index over yb with all the ensemble members finite and not NAs
  ix<-which(apply(yb0,MAR=1,
                  FUN=function(x){length(which(!is.na(x) & !is.nan(x)))})
            ==nens)
  n1<-length(ix)
  yb<-array(data=NA,dim=c(n1,nens))
  for (e in 1:nens) yb[,e]<-yb0[ix,e]
  rm(yb0)
  # background, probability of rain
  yb_pwet<-apply(yb,MAR=1,FUN=function(x){length(which(x>=argv$rrinf))/nens})
  # consider only observations where all yb members are finite and not NAs
  ix_orig<-ix
  yo<-yo[ix]
  VecX<-VecX[ix]
  VecY<-VecY[ix]
  prId<-prId[ix] 
  rm(ix)
  # yo_wet, 1=wet 0=dry
  yo_wet<-yo
  yo_wet[]<-0
  yo_wet[which(yo>argv$rrinf)]<-1
  # define the inverse of the diagonal observation error matrix
  # var(obs err) = eps2 * var(backg err)
  diagRinv<-rep(argv$hyletkf.eps2_prec_default,n1)
  for (i in 1:length(argv$hyletkf.eps2_prId)) {
    ix<-which(prId==argv$hyletkf.eps2_prId[i])
    if (length(ix)>0) diagRinv[ix]<-argv$hyletkf.eps2[i]
    rm(ix)
  }
  diagRinv[which(yo<argv$rrinf)]<-argv$hyletkf.eps2_noprec_default
  if (any(diagRinv==0 | is.na(diagRinv) | is.nan(diagRinv))) 
    boom("ERROR in the definition of Rinv")
  diagRinv<-1./diagRinv
  # if CVmode, then adjust the background
  if ((argv$cv_mode|argv$cv_mode_random)) {
    xb<-array(data=NA,dim=c(ncv,nens))
    r<-rmaster
    for (e in 1:nens) {
      r[]<-NA
      r[aix]<-xb_orig[,e]
      xb[,e]<-extract(r,cbind(VecX_cv,VecY_cv),method="bilinear")
      auxx<-which(xb[,e]<0 | is.na(xb[,e]) | is.nan(xb[,e]))
      if (length(auxx)>0) 
        xb[auxx,e]<-extract(r,cbind(VecX_cv[auxx],VecY_cv[auxx]))
      rm(auxx)
    }
    rm(r)
  }
  # Gaussian anamorphosis  
  if (argv$transf=="Box-Cox") {
    xb[which(xb<0)]<-0
    yb[which(yb<0)]<-0
    yo[which(yo<0)]<-0
    yo<-boxcox(yo,argv$transf.boxcox_lambda)
    yb<-boxcox(yb,argv$transf.boxcox_lambda)
    xb<-boxcox(xb,argv$transf.boxcox_lambda)
  }
  # ensemble mean
  if (!is.na(argv$cores)) {
    xbm<-mcmapply(function(x) mean(xb[x,]),1:length(xgrid),mc.cores=argv$cores,SIMPLIFY=T)
    ybm<-mcmapply(function(x) mean(yb[x,]),1:n0,mc.cores=argv$cores,SIMPLIFY=T)
  } else {
    xbm<-mapply(function(x) mean(xb[x,]),1:length(xgrid),SIMPLIFY=T)
    ybm<-mapply(function(x) mean(yb[x,]),1:n0,SIMPLIFY=T)
  }
  # ensemble perturbations
  Xb<-xb-xbm
  Yb<-yb-ybm
  # ensemble variances
  # var(backg err) is derived from the background ensemble
  # a minimum value of var(backg err) is set (over the original values)
  if (!is.na(argv$cores)) {
    ybvar<-mcmapply(function(x) mean(Yb[x,]**2),1:n0,mc.cores=argv$cores,SIMPLIFY=T)
  } else {
    ybvar<-mapply(function(x) mean(Yb[x,]**2),1:n0,SIMPLIFY=T)
  }
  if (argv$transf=="Box-Cox") {
    tybm<-tboxcox(ybm,argv$transf.boxcox_lambda)
    ybvar_ref<-ybvar
    sdmin<-sqrt(argv$hyletkf.sigma2_min)
    ix_tmp<-which(tybm>sdmin)
    if (length(ix_tmp)>0) 
      ybvar_ref[ix_tmp]<-0.5*(boxcox((tybm[ix_tmp]+sdmin),argv$transf.boxcox_lambda)-
                              boxcox((tybm[ix_tmp]-sdmin),argv$transf.boxcox_lambda))
    ix_tmp<-which(tybm<=sdmin)
    if (length(ix_tmp)>0) 
      ybvar_ref[ix_tmp]<-boxcox(sdmin,argv$transf.boxcox_lambda)
    ybvar<-pmax(ybvar,ybvar_ref**2)
    rm(tybm,ybvar_ref,ix_tmp)
  } else  {
  ybvar[ybvar<argv$hyletkf.sigma2_min]<-argv$hyletkf.sigma2_min
  }
  #
  #............................................................................
  # first step of hyletkf
  hyletkf_1<-function(i){
  # typeA vector VecX, VecY,... dimension 1:n0
  # typeB vector rloc, ... dimension 1:n.i
    deltax<-abs(xgrid[i]-VecX)
    deltay<-abs(ygrid[i]-VecY)
    if (!any(deltax<(7*argv$hyletkf.Dh))) return(rep(NA,nens))
    if (!any(deltay<(7*argv$hyletkf.Dh))) return(rep(NA,nens))
    ixa<-which( deltax<(7*argv$hyletkf.Dh) & 
                deltay<(7*argv$hyletkf.Dh) )
    # i-th gridpoint analysis 
    if (length(ixa)>0) {
      # exp (-1/2 * dist**2 / dh**2)
      rloc<-exp(-0.5* (deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]) / Dh2)
#      dist<-sqrt(deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa])
#      rloc<-(1+dist/argv$hyletkf.Dh)*exp(-dist/argv$hyletkf.Dh)
      if (length(ixa)>argv$hyletkf.pmax) {
        ixb<-order(rloc, decreasing=T)[1:argv$hyletkf.pmax]
        rloc<-rloc[ixb]
        ixa<-ixa[ixb]
        rm(ixb)
      }
      sigma2<-mean(ybvar[ixa])
      Yb.i<-Yb[ixa,,drop=F]
      d.i<-yo[ixa]-ybm[ixa]
      C.i<-t(1./sigma2*diagRinv[ixa]*rloc*Yb.i)
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
      xa<-xbm[i]+Xb[i,] %*% W
      rm(W,Wa,waa,a.eig)
    # i-th gridpoint is isolated
    } else {
      xa<-xb[i,]
    }
    return(xa)
  }
  if (argv$verbose) t00<-Sys.time()
  if (!is.na(argv$cores)) {
    xa<-t(mcmapply(hyletkf_1,1:length(xgrid),mc.cores=argv$cores,SIMPLIFY=T))
  } else {
    xa<-t(mapply(hyletkf_1,1:length(xgrid),SIMPLIFY=T))
  }
#  xa<-t(mapply(hyletkf_1,1:length(xgrid),SIMPLIFY=T))
  if (argv$verbose) {
    t11<-Sys.time()
    print(paste("hyletkf step1, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  }
  # analysis at observation locations
  if (!(argv$cv_mode|argv$cv_mode_random)) {
    r<-rmaster
    ya<-array(data=NA,dim=c(length(ix_orig),nens))
    for (e in 1:nens) {
      r[]<-NA
      r[aix]<-xa[,e]
      ya[,e]<-extract(r,cbind(VecX,VecY),method="bilinear")
      auxx<-which(ya[,e]<0 | is.na(ya[,e]) | is.nan(ya[,e]))
      if (length(auxx)>0) ya[auxx,e]<-extract(r,cbind(VecX[auxx],VecY[auxx]))
      rm(auxx)
    }
    rm(r)
  } else {
    xgrid_bak<-xgrid
    ygrid_bak<-ygrid
    Xb_bak<-Xb
    xbm_bak<-xbm
    xgrid<-VecX
    ygrid<-VecY
    Xb<-Yb
    xbm<-ybm
    if (argv$verbose) t00<-Sys.time()
    if (!is.na(argv$cores)) {
      ya<-t(mcmapply(hyletkf_1,1:length(xgrid),mc.cores=argv$cores,SIMPLIFY=T))
    } else {
      ya<-t(mapply(hyletkf_1,1:length(xgrid),SIMPLIFY=T))
    }
    if (argv$verbose) {
      t11<-Sys.time()
      print(paste("ya hyletkf step1, time=",round(t11-t00,1),attr(t11-t00,"unit")))
    }
    xgrid<-xgrid_bak
    ygrid<-ygrid_bak
    Xb<-Xb_bak
    xbm<-xbm_bak
    rm(xgrid_bak,ygrid_bak,Xb_bak,xbm_bak)
  }
  # analysis, probability of rain
  ya_pwet<-apply(ya,MAR=1,FUN=function(x){length(which(x>=argv$rrinf))/nens})
  #
  # ensemble mean
  if (!is.na(argv$cores)) {
    yam<-mcmapply(function(x) mean(ya[x,]),1:n0,mc.cores=argv$cores,SIMPLIFY=T)
  } else {
    yam<-mapply(function(x) mean(ya[x,]),1:n0,SIMPLIFY=T)
  }
  #
  #............................................................................
  # 2nd step of hyletkf
  hyletkf_2<-function(i){
  # typeA vector VecX, VecY,... dimension 1:n0
  # typeB vector rloc, ... dimension 1:n.i
    deltax<-abs(xgrid[i]-VecX)
    deltay<-abs(ygrid[i]-VecY)
    if (!any(deltax<(7*argv$hyletkf.Dh_oi))) return(rep(NA,nens))
    if (!any(deltay<(7*argv$hyletkf.Dh_oi))) return(rep(NA,nens))
    ixa<-which( deltax<(7*argv$hyletkf.Dh_oi) & 
                deltay<(7*argv$hyletkf.Dh_oi) )
    # i-th gridpoint analysis 
    if (length(ixa)>0) {
      # exp (-1/2 * dist**2 / dh**2)
      rloc<-exp(-0.5* (deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa]) / Dh_oi2)
#      dist<-sqrt(deltax[ixa]*deltax[ixa]+deltay[ixa]*deltay[ixa])
#      rloc<-(1+dist/Dh_oi)*exp(-dist/Dh_oi)
      if (length(ixa)>argv$hyletkf.pmax) {
        ixb<-order(rloc, decreasing=T)[1:argv$hyletkf.pmax]
        rloc<-rloc[ixb]
        ixa<-ixa[ixb]
        rm(ixb)
      }
      d.i<-yo[ixa]-yam[ixa]
      Disth2<-outer(VecY[ixa],VecY[ixa],FUN="-")**2.+
              outer(VecX[ixa],VecX[ixa],FUN="-")**2.
      D<-exp(-0.5*Disth2/Dh_oi2)
      rm(Disth2)
      diag(D)<-diag(D)+argv$hyletkf.eps2_oi
      InvD<-chol2inv(chol(D))
      xam<-mean(xa[i,])+as.numeric(crossprod(rloc,crossprod(t(InvD),d.i)))
      xa2<-xa[i,]-mean(xa[i,])+xam
    # i-th gridpoint is isolated
    } else {
      xa2<-xa[i,]
    }
    return(xa2)
  }
  if (!is.na(argv$cores)) {
    xa<-t(mcmapply(hyletkf_2,1:length(xgrid),mc.cores=argv$cores,SIMPLIFY=T))
  } else {
    xa<-t(mapply(hyletkf_2,1:length(xgrid),SIMPLIFY=T))
  }
  if (argv$verbose) {
    t11<-Sys.time()
    print(paste("hyletkf step2, time=",round(t11-t00,1),attr(t11-t00,"unit")))
  }
#  # intialize analysis vector
#  xa<-xb
#  xa[]<-NA
#  # letkf: loop over gridpoints
#  if (argv$verbose) t00<-Sys.time()
#  for (i in 1:length(xgrid)) {
#    dist<-(((xgrid[i]-VecX)**2+(ygrid[i]-VecY)**2 )**0.5)/1000.
#    rloc<-exp(-0.5*(dist**2./argv$hyletkf.Dh**2))
#    sel<-which(rloc>argv$hyletkf.rloc_min)
#    if ((i%%10000)==0) print(paste(i,length(sel)))
#    # i-th gridpoint analysis 
#    if (length(sel)>0) {
#      if (length(sel)>argv$hyletkf.pmax) sel<-order(dist)[1:argv$hyletkf.pmax]
#      sel_wet<-sel[which(yb_pwet[sel]>=0.5)]
#      if (length(sel_wet)>0) {
#        sigma2<-max(c(mean(ybvar[sel]),argv$hyletkf.sigma2_min))
#      } else {
#        sigma2<-argv$hyletkf.sigma2_min
#      }
#      Yb.i<-Yb[sel,,drop=F]
#      d.i<-yo[sel]-ybm[sel]
#      C.i<-t(1./sigma2*diagRinv[sel]*rloc[sel]*Yb.i)
#      C1.i<-crossprod(t(C.i),Yb.i)
#      Cd.i<-crossprod(t(C.i),d.i)
#      rm(C.i)
#      diag(C1.i)<-diag(C1.i) + (nens-1)
#      Pa.i<-chol2inv(chol(C1.i))
#      rm(C1.i)
#      # Pa*(k-1) is the matrix for which we want to square root:
#      a.eig <- eigen(Pa.i*(nens-1),symmetric=T)
#      Wa <- tcrossprod(
#             tcrossprod( a.eig$vectors, diag(sqrt(a.eig$values))),
#             t(solve(a.eig$vectors)) )
#      waa<-crossprod(Pa.i, Cd.i )
#      rm(Cd.i)
#      W<-Wa+as.vector(waa)
#      xa[i,]<-xbm[i]+Xb[i,] %*% W
#      rm(W,Wa,waa,a.eig)
#    # i-th gridpoint is isolated
#    } else {
#      xa[i,]<-xb[i,]
#    }
#  } # end: letkf loop over gridpoints
#  if (argv$verbose) {
#    t11<-Sys.time()
#    print(paste("letkf time=",round(t11-t00,1),attr(t11-t00,"unit")))
#  }
#  # OI
#  sel_oi<-which(yb_pwet<0.5 & yo_wet==1)  
#  if (length(sel_oi)>0) {
#    if (argv$verbose) t00<-Sys.time()
#    noi<-length(sel_oi)
#    Disth<-matrix(ncol=noi,nrow=noi,data=0.)
#    Disth<-(outer(VecY[sel_oi],VecY[sel_oi],FUN="-")**2.+
#            outer(VecX[sel_oi],VecX[sel_oi],FUN="-")**2.)**0.5/1000.
#    D<-exp(-0.5*(Disth/argv$hyletkf.Dh_oi)**2.)
#    diag(D)<-diag(D)+argv$hyletkf.eps2_oi
#    InvD<-chol2inv(chol(D))
#    xam<-OI_RR_fast(yo=yo[sel_oi],
#                    yb=ybm[sel_oi],
#                    xb=rowMeans(xa),
#                    xgrid=xgrid,
#                    ygrid=ygrid,
#                    VecX=VecX[sel_oi],
#                    VecY=VecY[sel_oi],
#                    Dh=argv$hyletkf.Dh_oi) 
#    xa<-xa-rowMeans(xa)+xam
#    rm(xam)
#    if (argv$verbose) {
#      t11<-Sys.time()
#      print(paste("oi time=",round(t11-t00,1),
#                             attr(t11-t00,"unit")))
#    }
#  }
  # Back-transformation
  if (argv$transf=="Box-Cox") {
    xa<-tboxcox(xa,argv$transf.boxcox_lambda)
    xb<-tboxcox(xb,argv$transf.boxcox_lambda)
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
# END of Hybrid Local Ensemble Transform Kalman Filter
