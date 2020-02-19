for (file in c(argv$off_y_verif_a,argv$off_y_verif_b,argv$off_y_verif_av,
               argv$off_yt_verif_a,argv$off_yt_verif_b,argv$off_yt_verif_av,
               argv$off_cv_verif_a,argv$off_cv_verif_b,
               argv$off_cvt_verif_a,argv$off_cvt_verif_b,
               argv$off_lcv_verif_a,argv$off_lcv_verif_b,
               argv$off_lcvt_verif_a,argv$off_lcvt_verif_b)) {
  if (is.na(file)) next
  if (file==argv$off_y_verif_a & !is.na(argv$off_y_verif_a)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_y
    obs_vals<-yo
    fcst_vals<-ya
  } else if (file==argv$off_y_verif_b & !is.na(argv$off_y_verif_b)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_y
    obs_vals<-yo
    fcst_vals<-yb
  } else if (file==argv$off_y_verif_av & !is.na(argv$off_y_verif_av)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_y
    obs_vals<-yo
    fcst_vals<-yav
  } else if (file==argv$off_yt_verif_a & !is.na(argv$off_yt_verif_a)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_y
    obs_vals<-yto
    fcst_vals<-yta
  } else if (file==argv$off_yt_verif_b & !is.na(argv$off_yt_verif_b)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_y
    obs_vals<-yto
    fcst_vals<-ytb
  } else if (file==argv$off_yt_verif_av & !is.na(argv$off_yt_verif_av)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_y
    obs_vals<-yto
    fcst_vals<-ytav
  } else if (file==argv$off_cv_verif_a & !is.na(argv$off_cv_verif_a)) {
    loc_vals<-VecS_cv
    lat_vals<-VecLat_cv
    lon_vals<-VecLon_cv
    elev_vals<-elev_for_verif_cv
    obs_vals<-yo_cv
    fcst_vals<-ya_cv
  } else if (file==argv$off_cv_verif_b & !is.na(argv$off_cv_verif_b)) {
    loc_vals<-VecS_cv
    lat_vals<-VecLat_cv
    lon_vals<-VecLon_cv
    elev_vals<-elev_for_verif_cv
    obs_vals<-yo_cv
    fcst_vals<-yb_cv
  } else if (file==argv$off_cvt_verif_a & !is.na(argv$off_cvt_verif_a)) {
    loc_vals<-VecS_cv
    lat_vals<-VecLat_cv
    lon_vals<-VecLon_cv
    elev_vals<-elev_for_verif_cv
    obs_vals<-yto_cv
    fcst_vals<-yta_cv
  } else if (file==argv$off_cvt_verif_b & !is.na(argv$off_cvt_verif_b)) {
    loc_vals<-VecS_cv
    lat_vals<-VecLat_cv
    lon_vals<-VecLon_cv
    elev_vals<-elev_for_verif_cv
    obs_vals<-yto_cv
    fcst_vals<-ytb_cv
  } else if (file==argv$off_lcv_verif_a & !is.na(argv$off_lcv_verif_a)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_lcv
    obs_vals<-yo_lcv
    fcst_vals<-ya_lcv
  } else if (file==argv$off_lcv_verif_b & !is.na(argv$off_lcv_verif_b)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_lcv
    obs_vals<-yo_lcv
    fcst_vals<-yb_lcv
  } else if (file==argv$off_lcvt_verif_a & !is.na(argv$off_lcvt_verif_a)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_lcv
    obs_vals<-yto_lcv
    fcst_vals<-yta_lcv
  } else if (file==argv$off_lcvt_verif_b & !is.na(argv$off_lcvt_verif_a)) {
    loc_vals<-VecS
    lat_vals<-VecLat
    lon_vals<-VecLon
    elev_vals<-elev_for_verif_lcv
    obs_vals<-yto_lcv
    fcst_vals<-ytb_lcv
  }
  nloc<-length(loc_vals)
  if (!exists("tstamp_nc")) {
    tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                      format="%Y%m%d%H%M",tz="GMT")
    dim_t<-list(name="time",
                units="H",
                vals=tstamp_nc,
                is_time=T,
                unlim=T,
                times.unit="S",
                times.format="%Y%m%d%H%M",
                timeref="197001010000",
                timeref.format="%Y%m%d%H%M")
  }
  dim_lt<-list(name="leadtime",units="",vals=0)
  dim_loc<-list(name="location",units="",vals=loc_vals)
  nlt<-1
  nt<-1
  var_lat<-list(name="lat",
                units="",
                dim_names=c("location"),
                vals=array(lat_vals,dim=c(nloc)))
  var_lon<-list(name="lon",
                units="",
                dim_names=c("location"),
                vals=array(lon_vals,dim=c(nloc)))
  var_ele<-list(name="altitude",
                units="",
                dim_names=c("location"),
                vals=array(elev_vals,dim=c(nloc)))
  obs<-array(data=NA,dim=c(nloc,nlt,nt))
  obs[,nlt,nt]<-obs_vals
  var_obs<-list(name="obs",
                units="",
                dim_names=c("location","leadtime","time"),
                vals=obs)
  rm(obs)
  if (!iff_is_ens) { 
    dim_list<-list(dim_t,dim_lt,dim_loc)
    fcst<-array(data=NA,dim=c(nloc,nlt,nt))
    fcst[,nlt,nt]<-fcst_vals
    var_fcst<-list(name="fcst",
                   units="",
                   dim_names=c("location","leadtime","time"),
                   vals=fcst)
    rm(fcst)
    var_list<-list(var_lat,var_lon,var_ele,var_obs,var_fcst)
    rm(var_lat,var_lon,var_ele,var_obs,var_fcst)
  # case of ensemble analysis
  } else {
    dim_ens<-list(name="ensemble_member",units="",vals=argv$iff_fg.e)
    nr<-length(argv$off_vernc.thresholds)
    nq<-length(argv$off_vernc.quantile)
    dim_r<-list(name="threshold",units="",vals=argv$off_vernc.thresholds)
    dim_q<-list(name="quantile",units="",vals=argv$off_vernc.quantile)
    dim_list<-list(dim_t,dim_lt,dim_loc,dim_ens,dim_r,dim_q)
    fcst<-array(data=NA,dim=c(nloc,nlt,nt))
    fcst[,nlt,nt]<-rowMeans(fcst_vals,na.rm=T)
    var_fcst<-list(name="fcst",
                   units="",
                   dim_names=c("location","leadtime","time"),
                   vals=fcst)
    rm(fcst)
    xval<-array(data=NA,dim=c(length(argv$iff_fg.e),nloc,nlt,nt))
    xval[,,nlt,nt]<-fcst_vals
    var_ens<-list(name="ensemble",
                  units="",
                  dim_names=c("ensemble_member","location","leadtime","time"),
                  vals=xval)
    rm(xval)
    pit<-array(data=NA,dim=c(nloc,nlt,nt))
    pit[,nlt,nt]<-apply(cbind(fcst_vals,obs_vals),
                        MARGIN=1,
                        FUN=function(x){
                         if (is.na(x[length(x)])) return(NA)
                         ix<-which(!is.na(x[1:(length(x)-1)]))
                         if (length(ix)==0) return(NA)
                         ecdf(x[ix])(x[length(x)])})
    var_pit<-list(name="pit",
                  units="",
                  dim_names=c("location","leadtime","time"),
                  vals=pit)
    rm(pit)
    cdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    pdf<-array(data=NA,dim=c(nr,nloc,nlt,nt))
    cdf[1:nr,1:nloc,nlt,nt]<-apply(fcst_vals, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        ecdf(x[ix])(argv$off_vernc.thresholds)})
    pdf[1:nr,1:nloc,nlt,nt]<-apply(fcst_vals, MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        approxfun(density(x[ix]),yleft=0,yright=0)(argv$off_vernc.thresholds)})
    var_cdf<-list(name="cdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=cdf)
    rm(cdf)
    var_pdf<-list(name="pdf",
                  units="",
                  dim_names=c("threshold","location","leadtime","time"),
                  vals=pdf)
    rm(pdf)
    qx<-array(data=NA,dim=c(nq,nloc,nlt,nt))
    qx[1:nq,1:nloc,nlt,nt]<-apply(fcst_vals,MARGIN=1,
      FUN=function(x){ix<-which(!is.na(x)); if (length(ix)==0) return(NA)
        as.vector(quantile(x[ix],probs=argv$off_vernc.quantile))})
    var_qx<-list(name="x",
                 units="",
                 dim_names=c("quantile","location","leadtime","time"),
                 vals=qx)
    rm(qx)
    var_list<-list(var_lat,var_lon,var_ele,var_obs,var_fcst,var_ens,
                   var_pit,var_cdf,var_pdf,var_qx)
    rm(var_lat,var_lon,var_ele,var_obs,var_fcst,var_ens,
       var_pit,var_cdf,var_pdf,var_qx)
  }
  gatt_longn<-list(attname="long_name",attval=argv$off_vernc.varname,prec="text")
  gatt_stdn<-list(attname="standard_name",attval=argv$off_vernc.stdvarname,prec="text")
  gatt_units<-list(attname="units",attval=argv$off_vernc.varunits,prec="text")
  res<-write_generic_dotnc(nc.file=file,
                           dims=dim_list,
                           vars=var_list,
                           glob_attrs=list(gatt_longn,gatt_stdn,gatt_units))
  if (is.null(res)) print(paste("ERROR while writing ",file))
  print(paste("output saved on file",file))
}
