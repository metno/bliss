# superobbing for the ith box (divide it in a n1 x n2 boxes) 
superobs<-function(ixy,         #i-index,x,y
                   res,         #resx,resy
                   n=c(1,1),
                   nmin=2) {  #nx,ny
  ixI<-which(VecI==ixy[1])
  if (n[1]==1 & n[2]==1) {
    nobs<-length(ixI)
    if (nobs>nmin) {
      nout<-1
      x_so<-mean(VecX[ixI])
      y_so<-mean(VecY[ixI])
      z_so<-mean(VecZ[ixI])
      laf_so<-mean(VecLaf[ixI])
      yo_so<-mean(yo[ixI])
      rm_me[ixI]<-T
      assign("rm_me",rm_me,envir=.GlobalEnv)
    } else {
      nout<-0
      x_so<-NULL
      y_so<-NULL
      z_so<-NULL
      laf_so<-NULL
      yo_so<-NULL
    }
  } else {
    r<-raster(xmn=ixy[2]-res[1]/2,xmx=ixy[2]+res[1]/2,
              ymn=ixy[3]-res[2]/2,ymx=ixy[3]+res[2]/2,
              res=c(res[1]/n[1],res[2]/n[2]))
    r[]<-1:ncell(r)
    nobs<-getValues(
           rasterize(cbind(VecX[ixI],VecY[ixI]),
                     r,
                     yo[ixI],
                     fun=function(x,...)length(x)) )
    ii<-extract(r,cbind(VecX[ixI],VecY[ixI]))
    ixnr<-which(!is.na(nobs) & nobs>nmin)
    nout<-length(ixnr)
    if (nout>0) {
      x_so<-vector()
      y_so<-vector()
      z_so<-vector()
      laf_so<-vector()
      yo_so<-vector()
      for (i in 1:nout) {
        x_so[i]<-mean(VecX[ixI][which(ii==ixnr[i])])
        y_so[i]<-mean(VecY[ixI][which(ii==ixnr[i])])
        z_so[i]<-mean(VecZ[ixI][which(ii==ixnr[i])])
        laf_so[i]<-mean(VecLaf[ixI][which(ii==ixnr[i])])
        yo_so[i]<-mean(yo[ixI][which(ii==ixnr[i])])
        rm_me[ixI][which(ii==ixnr[i])]<-T
      }
      assign("rm_me",rm_me,envir=.GlobalEnv)
    } else {
      x_so<-NULL
      y_so<-NULL
      z_so<-NULL
      laf_so<-NULL
      yo_so<-NULL
    }
  }
  return(list(n=nout,x=x_so,y=y_so,z=z_so,laf=laf_so,yo=yo_so))
}

