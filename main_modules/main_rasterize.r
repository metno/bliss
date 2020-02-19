  if (argv$verbose) print("Rasterize")
  r<-rasterize( x=cbind(VecX,VecY), y=rmaster, field=yo, fun=mean, na.rm=T )
  xr<-getValues(r)
  r<-rasterize( x=cbind(VecX,VecY), y=rmaster, field=yo, fun=sd, na.rm=T )
  xr_sd<-getValues(r)
  r<-rasterize( x=cbind(VecX,VecY), y=rmaster, field=yo, 
                fun=function(x,na.rm){length(na.omit(x))} )
  xr_n<-getValues(r)
  if (!is.na(argv$rasterize_nmin)) {
    ix_raster<-which(!is.na(xr_n) & xr_n<argv$rasterize_nmin)
    if (length(ix_raster)==0) {
      rm(ix_raster)
    } else {
      xr[ix_raster]<-NA
      xr_sd[ix_raster]<-NA
    }
    if (argv$verbose) {
      print(paste("    number of boxes with observations within=",
            length(which(!is.na(xr_n)))))
      print(paste("    number of masked boxes with observations within=",
            length(ix_raster)))
      print(paste(" -> number of unmasked boxes=",length(which(!is.na(xr)))))
    }
  }
  if (!any(is.na(argv$rasterize_q))) {
    xr_q<-array(data=NA,dim=c(length(xr_n),length(argv$rasterize_q)))
    for (i in 1:length(argv$rasterize_q)) {
      r<-rasterize( x=cbind(VecX,VecY), y=rmaster, field=yo, 
                    fun=function(x,na.rm){
                     quantile(na.omit(x),prob=argv$rasterize_q[i],type=4)} )
      xr_q[,i]<-getValues(r)
    }
    if (exists("ix_raster")) xr_q[ix_raster,]<-NA 
    aix<-1:length(xr)
  }

