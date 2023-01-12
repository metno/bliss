  if (!exists("r.list")) r.list<-list()
  r<-rmaster; r[]<-NA
  for (i in 1:length(argv$off_x.variables)) {
    if (argv$off_x.variables[i]=="analysis") {
      if (!exists("xa")) {xa<-aix;xa[]<-NA}
      xout<-xa
    } else if (argv$off_x.variables[i]=="analysis_errsd") {
      if (!exists("xa_errsd")) {xa_errsd<-aix;xa_errsd[]<-NA}
      xout<-xa_errsd
    } else if (argv$off_x.variables[i]=="background") {
      if (!exists("xb")) {xb<-aix;xb[]<-NA}
      xout<-xb
    } else if (argv$off_x.variables[i]=="rel_an") {
      if (!exists("xrela")) {xrela<-aix;xrela[]<-NA}
      xout<-xrela
    } else if (argv$off_x.variables[i]=="idi") {
      if (!exists("xidi")) {xidi<-aix;xidi[]<-NA}
      xout<-xidi
    } else if (argv$off_x.variables[i]=="dh" | 
               argv$off_x.variables[i]=="scale") {
      if (!exists("xdh")) {xdh<-aix;xdh[]<-NA}
      xout<-xdh
    } else if (argv$off_x.variables[i]=="gamma_shape") {
      if (!exists("a_gamma_shape")) {a_gamma_shape<-aix;a_gamma_shape[]<-NA}
      xout<-xa_pdf_par[,1]
    } else if (argv$off_x.variables[i]=="gamma_rate") {
      if (!exists("a_gamma_rate")) {a_gamma_rate<-aix;a_gamma_rate[]<-NA}
      xout<-xa_pdf_par[,2]
    } else if (argv$off_x.variables[i]=="varu") {
      if (!exists("xb_henoi_varu")) {xb_henoi_varu<-aix;xb_henoi_varu[]<-NA}
      xout<-xb_henoi_varu
    } else if (argv$off_x.variables[i]=="observations") {
      xout<-getValues(rasterize(x=cbind(VecX,VecY),y=rmaster,field=yo,fun=mean,na.rm=T))
    } else if (argv$off_x.variables[i]=="mean_raster") { 
      xout<-xr 
    } else if (argv$off_x.variables[i]=="sd_raster") { 
      xout<-xr_sd 
    } else if (argv$off_x.variables[i]=="n_raster") { 
      xout<-xr_n 
    } else if (substr(argv$off_x.variables[i],1,8)=="q_raster") { 
      if (length( ix_q<-which(as.integer(100*argv$rasterize_q) ==
           as.numeric(substr(argv$off_x.variables[i],10,11))))==0)
        boom(paste("problem with variable ",argv$off_x.variables[i],
                   "associated with quantile",
                   as.numeric(substr(argv$off_x.variables[i],10,11))))
      xout<-xr_q[,ix_q] 
    }
    if (is.array(xout)) {
      if (dim(xout)[1]==length(aix) & length(aix)!=ng) {
        xout1<-xout
        xout<-array(data=NA,dim=c(ng,dim(xout1)[2]))
        for (j in 1:dim(xout1)[2]) xout[aix,j]<-xout1[,j]
        rm(xout1)
      }
      r.list[[i]]<-array( data=NA, dim=c(length(x),length(y),dim(xout)[2],1) )
      for (j in 1:dim(xout)[2]) 
        r.list[[i]][,,j,1]<-matrix(xout[,j], ncol=length(y), nrow=length(x))
    } else {
      if (length(xout)==length(aix) & length(aix)!=ng) {
        xout1<-xout
        xout<-vector(mode="numeric",length=ng)
        xout[aix]<-xout1
        rm(xout1)
      }
      r.list[[i]]<-array(data=NA,dim=c(length(x),length(y),1))
      r.list[[i]][,,1]<-matrix(data=xout, ncol=length(y), nrow=length(x))
    }
  }
  # define time for output
  tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                    format="%Y%m%d%H%M",tz="GMT")
  time_bnds<-array(format(rev(seq(strptime(argv$date_out,argv$date_out_fmt),
                                           length=2,by=argv$time_bnds_string)),
                   format="%Y%m%d%H%M",tz="GMT"),dim=c(1,2))
  out<-write_dotnc(grid.list=r.list,
                   file.name=argv$off_x,
                   grid.type=argv$off_x.grid,
                   x=x,
                   y=y,
                   var.name=argv$off_x.varname,
                   var.longname=argv$off_x.varlongname,
                   var.standardname=argv$off_x.varstandardname,
                   var.version=argv$off_x.varversion,
                   var.unit=argv$off_x.varunit,
                   times=tstamp_nc,
                   times.unit=argv$off_x.timesunit,
                   reference=argv$off_x.reference,
                   proj4.string=argv$grid_master.proj4,
                   lonlat.out=argv$off_x.write_lonlat,
                   round.dig=argv$off_x.diground,
                   summary=argv$off_x.summary,
                   source.string=argv$off_x.sourcestring,
                   title=argv$off_x.title,
                   comment=argv$off_x.comment,
                   atts.var.add=NULL,
                   var.cell_methods=argv$off_x.cell_methods,
                   time_bnds=time_bnds,
                   cf_1.7=T)
  print(paste("output saved on file",argv$off_x))
