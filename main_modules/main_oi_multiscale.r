  if (nwet>0) {
    # if rescaling factor exists, multi-scale OI operates on relative anomalies 
    if (exists("yrf")) {
      if (any(yrf==0)) yrf[which(yrf==0)]<-1
      yo_relan<-yo/yrf
    } else {
      yo_relan<-yo
    }
    zero<-0
    if (argv$transf=="Box-Cox") {
      yo_relan<-boxcox(yo_relan,argv$transf.boxcox_lambda)
      zero<-boxcox(0,argv$transf.boxcox_lambda)
    } else if (argv$transf!="none") {
      boom("transformation not defined")
    }
    # multi-scale OI
    for (l in 1:nl) {
      if (argv$verbose) 
        print(paste("scale # ",formatC(l,width=2,flag="0"),
                    " of ",nl,
                    " (",formatC(vecd[l],width=4,flag="0"),"km ",
                    "fact=",formatC(vecf[l],width=3,flag="0"),")",sep=""))
      D<-exp(-0.5*(Disth/vecd[l])**2.)
      diag(D)<-diag(D)+vece[l]*ovarc
      InvD<-chol2inv(chol(D))
      if (l==nl | vecf[l]==1) {
        r<-rmaster
      } else {
        r<-aggregate(rmaster, fact=vecf[l], expand=T, na.rm=T)
      }
      zvalues.l<-getValues(r)
      storage.mode(zvalues.l)<-"numeric"
      xy.l<-xyFromCell(r,1:ncell(r))
      x.l<-sort(unique(xy.l[,1]))
      y.l<-sort(unique(xy.l[,2]),decreasing=T)
      mask.l<-which(!is.na(zvalues.l))
      zgrid.l<-zvalues.l[mask.l]
      xgrid.l<-xy.l[mask.l,1]
      ygrid.l<-xy.l[mask.l,2]
      rm(xy.l,zvalues.l)
      if (l==1) {
        yb<-rep(mean(yo_relan),length=n0)
        xb<-rep(mean(yo_relan),length=length(xgrid.l))
      } else {
        if ("scale" %in% argv$off_x.variables) 
          xl_tmp<-getValues(resample(rl,r,method="ngb"))[mask.l]
        rb<-resample(ra,r,method="bilinear")
        xb<-getValues(rb)[mask.l]
        count<-0
        while (any(is.na(xb))) {
          count<-count+1
          buffer_length<-round(vecd[l]/(10-(count-1)),0)*1000
          if (!is.finite(buffer_length)) break 
          ib<-which(is.na(xb))
          aux<-extract(rb,cbind(xgrid.l[ib],ygrid.l[ib]),
                       na.rm=T,
                       buffer=buffer_length)
          for (ll in 1:length(aux)) xb[ib[ll]]<-mean(aux[[ll]],na.rm=T)
          rb[mask.l]<-xb
          rm(aux,ib)
        }
        yb<-extract(rb,cbind(VecX,VecY),method="bilinear")
        rm(rb)
      }
      if (any(is.na(xb))) print("xb is NA")
      if (any(is.na(yb))) print("yb is NA")
      xa.l<-OI_RR_fast(yo=yo_relan,
                       yb=yb,
                       xb=xb,
                       xgrid=xgrid.l,
                       ygrid=ygrid.l,
                       VecX=VecX,
                       VecY=VecY,
                       Dh=vecd[l],
                       zero=zero) 
      ra<-r
      ra[]<-NA
      ra[mask.l]<-xa.l
      if ("scale" %in% argv$off_x.variables) {
        rl<-ra
        rl[mask.l]<-vecd[l]
        if (l>1) {
          ixl<-which(abs(xa.l-xb)<0.005)
          if (length(ixl)>0) rl[mask.l[ixl]]<-xl_tmp[ixl]
        }
      }
    } # end of multi-scale OI
    # back to precipitation values (from relative anomalies)
    if (argv$transf=="Box-Cox") {
      xa<-tboxcox(xa.l,argv$transf.boxcox_lambda)
    } else if (argv$transf!="none") {
      xa<-xa.l
    }
    if (exists("rf")) xa<-xa*rf[mask.l]
    print("refine prec/no-prec borders and remove wet regions with no obs in them")
    for (frac_rr in c(100,50,10)) {
      if (frac_rr==100) xa[which(!is.na(xa) & xa<(argv$rrinf/frac_rr))]<-0
      ra[mask.l]<-xa
      xa_aux<-getValues(ra)
      xa_aux[which(!is.na(xa_aux) & xa_aux<(argv$rrinf/frac_rr))]<-NA
      ra[]<-xa_aux
      rclump<-clump(ra)
      oclump<-extract(rclump,cbind(VecX[ixwet],VecY[ixwet]))
      fr<-freq(rclump)
      # remove clumps of YESprec cells less than (4x4)km^2 or not including wet obs
      ix<-which(!is.na(fr[,1]) & !is.na(fr[,2]) & ( (fr[,2]<=16) | !(fr[,1] %in% oclump)) )
      xa[which(getValues(rclump)[mask.l] %in% fr[ix,1])]<-0
      rm(xa_aux,rclump,oclump,fr,ix)
      ra[mask.l]<-xa
    }
    rm(mask.l,xa,xa.l)
    ya<-extract(ra,cbind(VecX,VecY),method="bilinear")
    yb<-rep(-9999,length(ya))
    yav<-rep(-9999,length(ya))
  } else {
    ra<-rmaster
    ra[]<-NA
    ra[mask]<-0
    rl<-ra
    ya<-rep(0,n0)
    yb<-rep(-9999,n0)
    yav<-rep(-9999,n0)
  }
  # compute IDI if needed
  if (!(argv$cv_mode | argv$cv_mode_random) & 
      ("idi" %in% argv$off_x.variables | 
       argv$idiv_instead_of_elev)) {
    xgrid.l<-xgrid[mask]
    ygrid.l<-ygrid[mask]
    if (argv$verbose) t00<-Sys.time()
    D<-exp(-0.5*((outer(VecY,VecY,FUN="-")**2.+
                  outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
                 /argv$oimult.Dh_idi)**2.)
    diag(D)<-diag(D)+argv$oimult.eps2_idi
    InvD<-chol2inv(chol(D))
    if (!argv$idiv_instead_of_elev) rm(D)
    if ("idi" %in% argv$off_x.variables) 
      xidi<-OI_RR_fast(yo=rep(1,length(VecX)),
                       yb=rep(0,length(VecX)),
                       xb=rep(0,length(xgrid.l)),
                       xgrid=xgrid.l,
                       ygrid=ygrid.l,
                       VecX=VecX,
                       VecY=VecY,
                       Dh=argv$oimult.Dh_idi)
    if (argv$idiv_instead_of_elev) {
      diag(D)<-diag(D)-argv$oimult.eps2_idi
      W<-tcrossprod(D,InvD)
      rm(D,InvD)
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
  } else {
    elev_for_verif<-VecZ
  }
  # prepare gridded output
  if (!(argv$cv_mode|argv$cv_mode_random)) {
    for (i in 1:length(argv$off_x.variables)) {
      if (!exists("r.list")) r.list<-list()
      if (argv$off_x.variables[i]=="analysis") {
        r.list[[i]]<-matrix(data=getValues(ra),
                            ncol=length(y),
                            nrow=length(x))
      } else if (argv$off_x.variables[i]=="background") {
        print("Warning: no background for OI_multiscale")
      } else if (argv$off_x.variables[i]=="idi") {
        r<-rmaster
        r[]<-NA
        r[mask]<-xidi
        r.list[[i]]<-matrix(data=getValues(r),
                            ncol=length(y),
                            nrow=length(x))
        rm(r)
      } else if (argv$off_x.variables[i]=="scale") {
        r.list[[i]]<-matrix(data=getValues(rl),
                            ncol=length(y),
                            nrow=length(x))
        rm(rl)
      } else if (argv$off_x.variables[i]=="rel_an") {
        r.list[[i]]<-matrix(data=getValues(ra)/rf,
                            ncol=length(y),
                            nrow=length(x))
      }
    }
  }
# END of OI multiscale (without background)
