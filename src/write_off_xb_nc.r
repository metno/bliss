#
write_off_xb_nc <- function( argv,  y_env, fg_env, u_env, env) {

  if (!exists("r.list")) r.list<-list()

  r <- env$rmaster
  r[]<-NA

  for (i in 1:length(argv$off_xb.variables)) {

    if (argv$off_xb.variables[i] == "background") {
#      if (!exists("xb")) {xb<-aix;xb[]<-NA}
#      xout <- xb
      xout <- array( data=NA, dim=c( env$ngrid, env$k_dim))
      for (e in 1:env$k_dim) {
        ii <- fg_env$ixs[e]
        rfxb <- subset( fg_env$fg[[fg_env$ixf[ii]]]$r_main, subset=fg_env$ixe[ii])
        rfxb[rfxb<y_env$rain] <- 0
        xout[,e] <- getValues(rfxb)[env$mask]
      } 
    } else if (argv$off_xb.variables[i] == "obs_align") {
      xout <- array( data=NA, dim=c( env$ngrid, 1))
      uo <- getValues( u_env$uo[[1]]$r_main)
      uo[uo<u_env$rain] <- 0
      uo[uo>u_env$rain] <- 1
      xout <- uo[env$mask]
    } else if (argv$off_xb.variables[i]=="obs_raster") {
      xout <- array( data=NA, dim=c( env$ngrid, 1))
      r <- env$rmaster
      r[] <- env$mergeobs$value
      xout <- getValues( r)[env$mask]
    } else if (argv$off_xb.variables[i]=="obs_idi_raster") {
      xout <- array( data=NA, dim=c( env$ngrid, 1))
      r <- env$rmaster
      r[] <- env$mergeobs$idi
      xout <- getValues( r)[env$mask]
#    } else if (argv$off_xb.variables[i]=="rel_an") {
#      if (!exists("xrela")) {xrela<-aix;xrela[]<-NA}
#      xout<-xrela
#    } else if (argv$off_xb.variables[i]=="idi") {
#      if (!exists("xidi")) {xidi<-aix;xidi[]<-NA}
#      xout<-xidi
    } else if ( argv$off_xb.variables[i]=="scale") {
      xout <- as.vector(env$Xscale)
#    } else if (argv$off_xb.variables[i]=="gamma_shape") {
#      if (!exists("a_gamma_shape")) {a_gamma_shape<-aix;a_gamma_shape[]<-NA}
#      xout<-xa_pdf_par[,1]
#    } else if (argv$off_xb.variables[i]=="gamma_rate") {
#      if (!exists("a_gamma_rate")) {a_gamma_rate<-aix;a_gamma_rate[]<-NA}
#      xout<-xa_pdf_par[,2]
#    } else if (argv$off_xb.variables[i]=="varu") {
#      if (!exists("xb_henoi_varu")) {xb_henoi_varu<-aix;xb_henoi_varu[]<-NA}
#      xout<-xb_henoi_varu
#    } else if (argv$off_xb.variables[i]=="observations") {
#      xout<-getValues(rasterize(x=cbind(VecX,VecY),y=rmaster,field=yo,fun=mean,na.rm=T))
#    } else if (argv$off_xb.variables[i]=="mean_raster") { 
#      xout<-xr 
#    } else if (argv$off_xb.variables[i]=="sd_raster") { 
#      xout<-xr_sd 
#    } else if (argv$off_xb.variables[i]=="n_raster") { 
#      xout<-xr_n 
#    } else if (substr(argv$off_xb.variables[i],1,8)=="q_raster") { 
#      if (length( ix_q<-which(as.integer(100*argv$rasterize_q) ==
#           as.numeric(substr(argv$off_xb.variables[i],10,11))))==0)
#        boom(paste("problem with variable ",argv$off_xb.variables[i],
#                   "associated with quantile",
#                   as.numeric(substr(argv$off_xb.variables[i],10,11))))
#      xout<-xr_q[,ix_q] 
    }
    if (is.array(xout)) {
      if ( env$ngrid != env$ng) {
        xout1 <- xout
        xout  <- array( data=NA, dim=c(env$ng,dim(xout1)[2]))
        for (j in 1:dim(xout1)[2]) xout[env$mask,j]<-xout1[,j]
        rm(xout1)
      }
      r.list[[i]]<-array( data=NA, dim=c(env$nx,env$ny,dim(xout)[2],1) )
      for (j in 1:dim(xout)[2]) 
        r.list[[i]][,,j,1]<-matrix(xout[,j], ncol=env$ny, nrow=env$nx)
    } else {
      if (env$ngrid != env$ng) {
        xout1 <- xout
        xout  <- vector(mode="numeric",length=env$ng)
        xout[env$mask] <- xout1
        rm(xout1)
      }
      r.list[[i]] <- array( data=NA, dim=c(env$nx,env$ny,1))
      r.list[[i]][,,1]<-matrix(data=xout, ncol=env$ny, nrow=env$nx)
    }
  }
  # define time for output
  tstamp_nc<-format(strptime(argv$date_out,argv$date_out_fmt),
                    format="%Y%m%d%H%M",tz="GMT")
  time_bnds<-array(format(rev(seq(strptime(argv$date_out,argv$date_out_fmt),
                                           length=2,by=argv$time_bnds_string)),
                   format="%Y%m%d%H%M",tz="GMT"),dim=c(1,2))
  out<-write_dotnc(grid.list=r.list,
                   file.name=argv$off_xb,
                   grid.type=argv$off_xb.grid,
                   x=env$x,
                   y=env$y,
                   var.name=argv$off_xb.varname,
                   var.longname=argv$off_xb.varlongname,
                   var.standardname=argv$off_xb.varstandardname,
                   var.version=argv$off_xb.varversion,
                   var.unit=argv$off_xb.varunit,
                   times=tstamp_nc,
                   times.unit=argv$off_xb.timesunit,
                   reference=argv$off_xb.reference,
                   proj4.string=argv$grid_master.proj4,
                   lonlat.out=argv$off_xb.write_lonlat,
                   round.dig=argv$off_xb.diground,
                   summary=argv$off_xb.summary,
                   source.string=argv$off_xb.sourcestring,
                   title=argv$off_xb.title,
                   comment=argv$off_xb.comment,
                   atts.var.add=NULL,
                   var.cell_methods=argv$off_xb.cell_methods,
                   time_bnds=time_bnds,
                   cf_1.7=T)
  print(paste("output saved on file",argv$off_xb))
}
