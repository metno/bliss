#+ header to c function
OI_RR_fast<-function(yo,
                       yb,
                       xb,
                       xgrid,
                       ygrid,
                       zgrid=NULL,
                       VecX,
                       VecY,
                       VecZ=NULL,
                       Dh,
                       Dz=NULL,
                       zero=0) {
#------------------------------------------------------------------------------
  if (is.null(Dz)) {
    Dz<-10000000
    VecZ<-rep(0,length(VecX))
    zgrid<-rep(0,length(xgrid))
  }
  no<-length(yo)
  ng<-length(xb)
  xa<-vector(mode="numeric",length=ng)
  vec<-vector(mode="numeric",length=no)
  d<-yo-yb
  out<-.C("oi_rr_first",no=as.integer(no), 
                        innov=as.double(d),
                        SRinv=as.numeric(InvD),
                        vec=as.double(vec) ) 
  vec[1:no]<-out$vec[1:no]
  rm(out)
  out<-.C("oi_rr_fast",ng=as.integer(ng),
                       no=as.integer(no),
                       xg=as.double(xgrid),
                       yg=as.double(ygrid),
                       zg=as.double(zgrid),
                       xo=as.double(VecX),
                       yo=as.double(VecY),
                       zo=as.double(VecZ),
                       Dh=as.double(Dh),
                       Dz=as.double(Dz),
                       xb=as.double(xb),
                       vec=as.double(vec),
                       xa=as.double(xa),
                       zero=as.double(zero) )
  out$xa[1:ng]
}

