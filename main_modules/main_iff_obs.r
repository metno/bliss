dat<-read.table(file=argv$iff_obs,
                header=T,
                sep=argv$iff_obs.sep,
                stringsAsFactors=F,
                strip.white=T)
# read x_orig,y_orig,value & set x,y 
varidxtmp<-match(c(argv$iff_obs.x,
                   argv$iff_obs.y,
                   argv$iff_obs.value),
                   names(dat))
if (any(is.na(varidxtmp))) {
  print("ERROR in the specification of the variable names")
  print(paste("        x=",argv$iff_obs.x))
  print(paste("        y=",argv$iff_obs.y))
  print(paste("    value=",argv$iff_obs.value))
  print("header of input file:")
  print(argv$argv$iff_obs)
  print(names(dat))
  quit(status=1)
}
data<-data.frame(dat[,varidxtmp])
names(data)<-c("x_orig","y_orig","value")
data$x_orig<-suppressWarnings(as.numeric(data$x_orig))
data$y_orig<-suppressWarnings(as.numeric(data$y_orig))
data$value<-suppressWarnings(as.numeric(data$value))
if (argv$iff_obs.proj4!=argv$grid_master.proj4) {
  xymaster<-spTransform(SpatialPoints(cbind(data$x_orig,data$y_orig),
                                      proj4string=CRS(argv$iff_obs.proj4)) ,
                        CRS(argv$grid_master.proj4))
  data$x<-attr(xymaster,"coords")[,1]
  data$y<-attr(xymaster,"coords")[,2]
  rm(xymaster)
} else {
  data$x<-data$x_orig
  data$y<-data$y_orig
}
# lat-lon coords are required by verif
if (argv$iff_obs.proj4!=proj4.llwgs84) {
  xyll<-spTransform(SpatialPoints(cbind(data$x_orig,data$y_orig),
                                  proj4string=CRS(argv$iff_obs.proj4)) ,
                    CRS(proj4.llwgs84))
  data$lon<-attr(xyll,"coords")[,1]
  data$lat<-attr(xyll,"coords")[,2]
  rm(xyll)
} else {
  data$lon<-data$x_orig
  data$lat<-data$y_orig
}
# read z 
if (argv$iff_obs.z!="none") {
  varidxtmp<-which(argv$iff_obs.z==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     z=",argv$iff_obs.z))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$z<-dat[,varidxtmp]
} else {
  data$z<-rep(0,length(data$x))
}
# read sourceId 
if (argv$iff_obs.sourceId!="none") {
  varidxtmp<-which(argv$iff_obs.sourceId==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     sourceId=",argv$iff_obs.sourceId))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$sourceId<-dat[,varidxtmp]
} else {
  data$sourceId<-1:length(data$x)
}
# read prId 
if (argv$iff_obs.prId!="none") {
  varidxtmp<-which(argv$iff_obs.prId==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     prId=",argv$iff_obs.prId))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$prId<-dat[,varidxtmp]
} else {
  data$prId<-rep(0,length(data$x))
}
# read dqc flag 
if (argv$iff_obs.dqc!="none") {
  varidxtmp<-which(argv$iff_obs.dqc==names(dat))
  if (length(varidxtmp)==0) {
    print("ERROR in the specification of the variable names")
    print(paste("     dqc=",argv$iff_obs.dqc))
    print("header of input file:")
    print(argv$iff_obs)
    print(names(dat))
    quit(status=1)
  }
  data$dqc<-dat[,varidxtmp]
} else {
  data$dqc<-rep(0,length(data$x))
}
# data, dataframe: x,y,z,x_orig,y_orig,value,sourceId,prId,dqc
#
if (file.exists(argv$iff_black)) {
  bstid<-read.csv(argv$iff_black,header=T,stringsAsFactors=F,strip.white=T)
} else {
  bstid<-integer(0)
}
# select only observation within the master grid
flag_in_master<-!is.na(extract(rmaster,cbind(data$x,data$y)))
flag_in_fg<-rep(T,length(data$x))
if (file.exists(argv$iff_fg)) 
  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) 
# on-the-fly dqc, used for testing
#  flag_in_fg<-!is.na(extract(rfg,cbind(data$x,data$y))) &
#              data$value > (extract(rfg,cbind(data$x,data$y))-0.2*extract(rfg,cbind(data$x,data$y)))
#CVmode
if (cv_mode) {
# prId=1 MET-WMO stations
  ixcv<-which( data$dqc==0 & 
               data$prId %in% argv$prId.cv & 
               !is.na(data$value) &
               flag_in_master &
               flag_in_fg   )
  if (length(ixcv)==0) boom("ERROR \"cv_mode\" running without CV-observations ")
  VecX_cv<-data$x[ixcv]
  VecY_cv<-data$y[ixcv]
  VecXorig_cv<-data$x_orig[ixcv]
  VecYorig_cv<-data$y_orig[ixcv]
  VecLat_cv<-data$lat[ixcv]
  VecLon_cv<-data$lon[ixcv]
  VecZ_cv<-data$z[ixcv]
  VecS_cv<-data$sourceId[ixcv]
  yo_cv<-data$value[ixcv]
  if (exists("rrf")) yrf_cv<-extract(rrf,cbind(VecX_cv,VecY_cv),na.rm=T)
  data$value[ixcv]<-NA
  data$dqc[ixcv]<-999
  ncv<-length(ixcv)
  if (!is.na(argv$rrinf)) {
    ixwet_cv<-which(yo_cv>=argv$rrinf)
    ixdry_cv<-which(yo_cv< argv$rrinf)
    nwet_cv<-length(ixwet_cv)
    ndry_cv<-length(ixdry_cv)
  }
} else if (cv_mode_random) {
  stmp<-vector()
  xtmp<-vector()
  ytmp<-vector()
  ncv<-0
  for (i in 1:length(data$x)) {
    if ( is.na(data$value[i]) | data$dqc[i]!=0 | !flag_in_master[i] | !flag_in_fg[i] ) next
    if (ncv==0) {
      ncv<-ncv+1
      stmp[ncv]<-data$sourceId[i]
      xtmp[ncv]<-data$x[i]
      ytmp[ncv]<-data$y[i]
    } else {
      dtmp<-sqrt( (data$x[i]-xtmp[1:ncv])**2 + (data$y[i]-ytmp[1:ncv])**2 )/1000
      if (!any(dtmp<argv$d.cv)) {
        ncv<-ncv+1
        stmp[ncv]<-data$sourceId[i]
        xtmp[ncv]<-data$x[i]
        ytmp[ncv]<-data$y[i]
      }
    }
  }
  ixcv<-which( data$sourceId %in% stmp[1:ncv] )
  if (length(ixcv)==0) boom("ERROR \"cv_mode\" running without CV-observations ")
  VecX_cv<-data$x[ixcv]
  VecY_cv<-data$y[ixcv]
  VecXorig_cv<-data$x_orig[ixcv]
  VecYorig_cv<-data$y_orig[ixcv]
  VecLat_cv<-data$lat[ixcv]
  VecLon_cv<-data$lon[ixcv]
  VecZ_cv<-data$z[ixcv]
  VecS_cv<-data$sourceId[ixcv]
  yo_cv<-data$value[ixcv]
  if (exists("rrf")) yrf_cv<-extract(rrf,cbind(VecX_cv,VecY_cv),na.rm=T)
  data$value[ixcv]<-NA
  data$dqc[ixcv]<-999
  ncv<-length(ixcv)
  if (!is.na(argv$rrinf)) {
    ixwet_cv<-which(yo_cv>=argv$rrinf)
    ixdry_cv<-which(yo_cv< argv$rrinf)
    nwet_cv<-length(ixwet_cv)
    ndry_cv<-length(ixdry_cv)
  }
}
if (any(!is.na(argv$prId.exclude))) {
  ix0<-which(data$dqc==0 & 
             !(data$prId %in% argv$prId.exclude) &
             flag_in_master &
             flag_in_fg &
             !is.na(data$value) &
             !is.nan(data$value) )
} else {
  ix0<-which(data$dqc==0 &
             flag_in_master &
             flag_in_fg  &
             !is.na(data$value) &
             !is.nan(data$value) )
}
n0<-length(ix0)
# definitive station list
if (n0==0) boom("No observations. Stop here.")
VecX<-data$x[ix0]
VecY<-data$y[ix0]
VecXorig<-data$x_orig[ix0]
VecYorig<-data$y_orig[ix0]
VecLat<-data$lat[ix0]
VecLon<-data$lon[ix0]
VecZ<-data$z[ix0]
VecS<-data$sourceId[ix0]
yo<-data$value[ix0]
if (exists("rrf"))  yrf<-extract(rrf,cbind(VecX,VecY),na.rm=T)
if (exists("rlaf")) VecLaf<-extract(rlaf,cbind(VecX,VecY),na.rm=T)
prId<-data$prId[ix0]
ydqc.flag<-rep(0,length=n0)
rm(data)
if (!is.na(argv$rrinf)) {
  ixwet<-which(yo>=argv$rrinf)
  ixdry<-which(yo< argv$rrinf)
  nwet<-length(ixwet)
  ndry<-length(ixdry)
}
# Aggregate observations
if (argv$obspp.agg) {
  for (i in 1:length(argv$obspp.agg_prId)) {
    if (is.na(argv$obspp.agg_prId[i])) {
      ix<-1:n0
    } else {
      ix<-which(prId==argv$obspp.agg_prId[i])
    }
    if (length(ix)==0) next
    if (argv$obspp.agg_fact[i]>1) {
      r<-aggregate(rmaster, fact=argv$obspp.agg_fact[i], expand=T, na.rm=T)
    } else {
      r<-rmaster
    }
    r1<-rasterize(x=cbind(VecX[ix],VecY[ix]),
                  y=r,
                  field=yo[ix],
                  fun=mean,
                  na.rm=T)
    if (length(ix)!=n0) {
      VecX<-VecX[-ix]
      VecY<-VecY[-ix]
      VecXorig<-VecXorig[-ix]
      VecYorig<-VecYorig[-ix]
      VecLat<-VecLat[-ix]
      VecLon<-VecLon[-ix]
      VecZ<-VecZ[-ix]
      VecS<-VecS[-ix]
      yo<-yo[-ix]
      if (exists("rrf"))  yrf<-yrf[-ix]
      if (exists("rlaf")) VecLaf<-VecLaf[-ix]
      prId<-prId[-ix]
    }
    ix1<-which( !is.na(getValues(r1)) )
    if (exists("rrf")) ix1<-ix1[which(!is.na( extract(rrf,
                                               cbind(xyFromCell(r1,ix1)[,1],
                                                     xyFromCell(r1,ix1)[,2]),na.rm=T)))]
    if (length(ix1)==0) next
    v1<-getValues(r1)[ix1]
    x1<-xyFromCell(r1,ix1)[,1]
    y1<-xyFromCell(r1,ix1)[,2]
    if (argv$iff_obs.proj4!=argv$grid_master.proj4) {
      xymaster<-spTransform(SpatialPoints(cbind(x1,y1),
                                          proj4string=CRS(argv$grid_master.proj4)) ,
                            CRS(argv$iff_obs.proj4))
      x1_orig<-attr(xymaster,"coords")[,1]
      y1_orig<-attr(xymaster,"coords")[,2]
      rm(xymaster)
    } else {
      x1_orig<-x1
      y1_orig<-y1
    }
    if (proj4.llwgs84!=argv$grid_master.proj4) {
      xyll<-spTransform(SpatialPoints(cbind(x1,y1),
                                      proj4string=CRS(argv$grid_master.proj4)) ,
                        CRS(proj4.llwgs84))
      x1_ll<-attr(xyll,"coords")[,1]
      y1_ll<-attr(xyll,"coords")[,2]
      rm(xyll)
    } else {
      x1_ll<-x1
      y1_ll<-y1
    }
    VecX<-c(VecX,x1)
    VecY<-c(VecY,y1)
    VecXorig<-c(VecXorig,x1_orig)
    VecYorig<-c(VecYorig,y1_orig)
    VecLat<-c(VecLat,y1_ll)
    VecLon<-c(VecLon,x1_ll)
    VecZ<-c(VecZ,rep(0,length(x1)))
    VecS<-c(VecS,rep(NA,length(x1)))
    yo<-c(yo,v1)
    if (exists("rrf"))  yrf<-c(yrf,extract(rrf,cbind(x1,y1),na.rm=T))
    if (exists("rlaf")) VecLaf<-c(VecLaf,extract(rlaf,cbind(x1,y1),na.rm=T))
    prId<-c(prId,rep(argv$obspp.agg_prId[i],length(x1)))
  } # end loop over argv$obspp.agg_prId
  n0<-length(yo)
  ydqc.flag<-rep(0,length=n0)
  if (!is.na(argv$rrinf)) {
    ixwet<-which(yo>=argv$rrinf)
    ixdry<-which(yo< argv$rrinf)
    nwet<-length(ixwet)
    ndry<-length(ixdry)
  }
  rm(x1,y1,x1_ll,y1_ll,x1_orig,y1_orig,r,r1,ix)
} # end observation aggregation
#
# ad-hoc data quality control
if (argv$obspp.dqcpuddle) {
  if (argv$mode=="OI_firstguess") {
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
    # dh (m), one for all gridpoints
    dh<-argv$oifg.Dh*1000
    dh2<-dh*dh
    # background at station locations
    yb<-extract(rfg, cbind(VecX,VecY), method="bilinear")
    # create a coarser grid
    if (argv$obspp.dqcpuddle_fact>1) {
      r<-aggregate(rfg, fact=argv$obspp.dqcpuddle_fact, expand=T, na.rm=T)
    } else {
      r<-rmaster
    }
    xix<-which(!is.na(getValues(r)))
    # check: at least 10 events are needed on the grid
    if (length(xix)>10) {
      xgrid_spint<-xyFromCell(r,xix)[,1]
      ygrid_spint<-xyFromCell(r,xix)[,2]
      ngrid_spint<-length(ygrid_spint)
      # -- elaboration on the grid, x_elab --
      t00<-Sys.time()
      if (argv$transf=="Box-Cox") {
        xb_spint<-getValues(r)[xix]
        yo_spint<-yo
        yb_spint<-yb
        xb_spint[which(xb_spint<0)]<-0
        yb_spint[which(yb_spint<0)]<-0
        yo_spint[which(yo_spint<0)]<-0
        yo_spint<-boxcox(yo_spint,argv$transf.boxcox_lambda)
        yb_spint<-boxcox(yb_spint,argv$transf.boxcox_lambda)
        xb_spint<-boxcox(xb_spint,argv$transf.boxcox_lambda)
      } else {
        xb_spint<-getValues(r)[xix]
        yo_spint<-yo
        yb_spint<-yb
      }
      # multicores
      # xa=arr[,1], xa_errvar=arr[,2], o_errvar=arr[,3], xidi=arr[,4], idiv=arr[,5]
      if (!is.na(argv$cores)) {
        arr<-t(mcmapply(oi_var_gridpoint_by_gridpoint,
                        1:ngrid_spint,
                        mc.cores=argv$cores,
                        SIMPLIFY=T,
                        pmax=argv$pmax,
                        corr=argv$corrfun))
      # no-multicores
      } else {
        arr<-t(mapply(oi_var_gridpoint_by_gridpoint,
                      1:ngrid_spint,
                      SIMPLIFY=T,
                      pmax=argv$pmax,
                      corr=argv$corrfun))
      }
      if (argv$verbose) {
        t11<-Sys.time()
        print(paste("obspp_dqcpuddle elab, time=",
                    round(t11-t00,1),attr(t11-t00,"unit")))
      }
      # loop over dqc thresholds
      for (i in 1:length(argv$obspp.dqcpuddle_thres)) {
        minarea<-ceiling(1000000*argv$obspp.dqcpuddle_minarea[i]/(xres(r)*yres(r)))
        r[]<-NA
        thr<-argv$obspp.dqcpuddle_thres[i]
        if (argv$transf=="Box-Cox") 
          thr<-boxcox(argv$obspp.dqcpuddle_thres[i],argv$transf.boxcox_lambda)
        # T= event occurred; F=event not occurred
        if (argv$obspp.dqcpuddle_cond[i]=="lt") {
          xcond<-arr[,1]<thr
          ycond<-yo<argv$obspp.dqcpuddle_thres[i]
        } else if (argv$obspp.dqcpuddle_cond[i]=="le") {
          xcond<-arr[,1]<=thr
          ycond<-yo<=argv$obspp.dqcpuddle_thres[i]
        } else if (argv$obspp.dqcpuddle_cond[i]=="gt") {
          xcond<-arr[,1]>thr 
          ycond<-yo>argv$obspp.dqcpuddle_thres[i]
        } else if (argv$obspp.dqcpuddle_cond[i]=="ge") {
          xcond<-arr[,1]>=thr
          ycond<-yo>=argv$obspp.dqcpuddle_thres[i]
        }
        if (!any(xcond)) next
        r[xix[which(xcond)]]<-1
        yix<-which(ycond)
        if (length(yix)==0) next
        rc<-clump(r)
        dc<-getValues(rc)
        f<-freq(rc)
        ytmp<-extract(rc,cbind(VecX[yix],VecY[yix]),na.rm=T,small=T)
        if (!any(!is.na(ytmp))) next
        nobs_in_f<-mapply(function(x){
                           if (is.na(f[x,1])) {
                             length(which(is.na(ytmp)))
                           } else {
                             length(which(!is.na(ytmp) & ytmp==f[x,1]))
                           }},
                           1:length(f[,1]),
                           SIMPLIFY=T)
        # observation is suspect if:
        # a. the region where the event occurred is too small
        # b. the are too few observation in the region where the event occurred
        ixsus<-yix[which(
            is.na(ytmp) |
            f[match(ytmp,f),2]<minarea | # a.
            nobs_in_f[match(ytmp,f)]<argv$obspp.dqcpuddle_nminobs[i])] # b.
        if (length(ixsus)==0) next
        VecX<-VecX[-ixsus]
        VecY<-VecY[-ixsus]
        VecXorig<-VecXorig[-ixsus]
        VecYorig<-VecYorig[-ixsus]
        VecLat<-VecLat[-ixsus]
        VecLon<-VecLon[-ixsus]
        VecZ<-VecZ[-ixsus]
        VecS<-VecS[-ixsus]
        yo<-yo[-ixsus]
        if (exists("rrf"))  yrf<-yrf[-ixsus]
        if (exists("rlaf")) VecLaf<-VecLaf[-ixsus]
        prId<-prId[-ixsus]
      } # loop over dqc thresholds
      n0<-length(yo)
      ydqc.flag<-rep(0,length=n0)
      if (!is.na(argv$rrinf)) {
        ixwet<-which(yo>=argv$rrinf)
        ixdry<-which(yo< argv$rrinf)
        nwet<-length(ixwet)
        ndry<-length(ixdry)
      }
    } # endif, check: at least 10 events are needed on the grid 
  } #end if (argv$mode=="OI_firstguess") {
}
# observation error variance correction factor
ovarc<-rep(1,n0)
if (any(!is.na(argv$ovarc.prId))) {
  for (i in 1:length(argv$ovarc.prId)) {
    if (any(prId==argv$ovarc.prId[i])) 
      ovarc[which(prId==argv$ovarc.prId[i])]<-argv$ovarc[i]
  }
}
if (argv$verbose) { 
  print("+---------------------------------------------------------------+")
  if (!is.na(argv$rrinf)) {
    print(paste("#observations (wet/dry) =",n0,"(",nwet,"/",ndry,")"))
    if (cv_mode | cv_mode_random) {
      print(paste("#cv-observations (wet/dry) =",ncv,"(",nwet_cv,"/",ndry_cv,")"))
    }
  } else {
    print(paste("#observations =",n0))
  }
  print("+...............................................................+")
}
if (!is.na(argv$off_obspp)) {
  dataout<-array(data=NA, dim=c(length(VecLat),8))
  names(dataout)<-c("lat","lon","elev","value","prid","dqc","sct","rep")
  dataout[,1]<-round(VecLat,6)
  dataout[,2]<-round(VecLon,6)
  dataout[,3]<-round(VecZ,0)
  dataout[,4]<-round(yo,3)
  dataout[,5]<-prId
  dataout[,6]<-round(ydqc.flag,0)
  dataout[,7]<-0
  dataout[,8]<-round(eps2[-ixsus],6)
  write.table(file=argv$off_obspp,
              dataout,
              quote=F,
              col.names=c("lat","lon","elev","value","prid","dqc","sct","rep"),
              row.names=F,
              dec=".",
              sep=";")
  print(paste("output saved on file",argv$off_obspp))
  quit(status=0)
}
